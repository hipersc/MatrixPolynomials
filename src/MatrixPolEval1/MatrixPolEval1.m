function [c_prec,leading_coeff_sign,c_vpa,er_min,all_cvpa,savings,s,p]=MatrixPolEval1(b,precision,type_pol,ndigits,s)
% MatrixPolEval1 computes evaluation formulas that save one matrix product with
% respect to the Paterson-Stockmeyer method [1, 2] for polynomials with
% a non-zero leading coefficient. No savings are possible for polynomials
% with degrees m <= 7, m = 9, and m = 11, and the Paterson-Stockmeyer method 
% is recommended for polynomials with those degrees. The function rejects complex 
% solutions and real solutions whose leading coefficients are negative, since 
% there are similar solutions with positive leading coefficients. The function 
% warns whether the accuracy of the evaluation formulas is likely to be low 
% by using the stability test from [2, Ex. 3.1].
%
% [c_prec,leading_coeff_sign,c_vpa,er_min,all_cvpa,savings,s,p] = MatrixPolEval1(b,precision,type_pol,ndigits,s)
%
% This implementation considers polynomial z1ps(A) (38) from [1].
%
% Input data:
% - b: polynomial coefficients b(1)*I + b(2)*A + b(3)*A^2 + ... + b(m+1)*A^m
%   (default: b = 1./factorial(sym(0:8)); exponential Taylor approximation of
%   degree 8).
% - precision: 'single' or 'double' IEEE precision arithmetic (default: 'double').
% - type_pol: Three possibilities for the evaluation formula y1s(A) are given
%   (default: 1, corresponding to eq. (35) from [1]).
% - ndigits: variable precision arithmetic number of digits (default: 32).
% - s: The polynomial degree must be m = 4s + p with p >= 0. (default: 
%   minimum value of s to save one matrix product with respect to the
%   Paterson-Stockmeyer method).
%   
% Output data:
% - c_prec: coefficients in single or double precision.
% - leading_coeff_sign: +1 or -1. If b_{m+1} < 0, leading_coeff_sign = -1 and the
%   function computes the coefficients for y1s(A) = -Pm(A) to avoid complex 
%   solutions for c_{m+1} (see [1] and [2, p. 237]).
% - c_vpa: coefficients with 'ndigits' variable precision arithmetic.
% - er_min: minimum of the maximum reconstruction relative error of
%   coefficients b.
% - all_cvpa: All the coefficients obtained with 'ndigits' variable precision
%   arithmetic.
% - savings: Indicates the number of matrix products saved compared to the
%   Paterson-Stockmeyer method for the selected values of s and p. This metric
%   helps verify whether actual savings are achieved. Its value is less than 
%   or equal to 1.
% - s, p: Selected values of s and p such that m = 4s + p with p >= 0.
%
% Examples: 
% % Best coefficients for evaluating polynomial I+A+A^2+...+A^m saving one 
% % matrix product with respect to the Paterson-Stockmeyer method for m=8, 10 or m>=12:
% c_prec = MatrixPolEval1(ones(m+1,1)) 
%
% % Full syntax call:
% [c_prec,leading_coeff_sign,c_vpa,er_min,all_cvpa,savings,s,p] = MatrixPolEval1(b,precision,type_pol,ndigits,s)
%
% % Example using default b, double precision, type_pol=1, and 64 digits:
% [c_prec,leading_coeff_sign,c_vpa,er_min,all_cvpa,savings,s,p] = MatrixPolEval1([],'double',1,64)
%
% References:
% [1] J.M. Alonso, J. Sastre, J. Ibáñez, E. Defez, A general framework for 
%     matrix polynomial evaluation with reduced computational cost and 
%     improved stability. (Details that no savings are possible for 
%     m <= 7, 9, and 11, and introduces type 2 and 3 expressions).
% [2] J. Sastre, Efficient evaluation of matrix polynomials, Linear 
%     Algebra Appl., 539, 2018, 229-250.
%
% Authors: José Miguel Alonso, Jorge Sastre
% Revised version: 2025/06/13.
%
% Group of High Performance Scientific Computing (HiPerSC)
% Universitat Politecnica de Valencia (Spain)
% http://hipersc.blogs.upv.es

% Check input data
if nargin<4 || isempty(ndigits)
    ndigits=32;
end
if nargin<3 || isempty(type_pol)
    type_pol=1;
end
if nargin<2 || isempty(precision)
    precision='double';
end 
if nargin<1 || isempty(b)
    % Example: b0, b1,..., b8 coefficients for Taylor approximation to
    % the matrix exponential
    m=8;
    b=1./factorial(sym(0:8));
end

% Get the polynomial order
m=length(b)-1;
if m<=7 || m==9 || m==11
    warning('MatrixPolEval:DegreeNotSupported', ...
            'Matrix product savings are not available for the specified polynomial degree.\nPlease apply the Paterson-Stockmeyer method directly.')
    c_prec=[]; c_vpa=[]; er_min=[]; all_cvpa=[];
    return
end

if b(m+1) == 0
    error('MatrixPolEval:ZeroLeadingCoefficient', ...
          'The leading coefficient of the polynomial must be non-zero.');
elseif b(m+1) < 0
    warning('MatrixPolEval:NegativeLeadingCoefficient', ...
            'Negative polynomial leading coefficient found: computing the coefficients for y1s(A) = - Pm(A), leading_coeff_sign = -1');
    b = -b;
    leading_coeff_sign = -1;
else
    leading_coeff_sign = 1;
end

% Check the polynomial type
if type_pol<1 || type_pol>3
    type_pol=1;
end

% Get parameters s and p from m (remember that m = 4s + p)
if nargin < 5 || isempty(s)
    % Initial estimation for s
    s = floor(sqrt(m));
    p = m - 4*s;
    CostMatrixPolEval1 = s + 1 + ceil(p/s);
    
    % Check the cost for a lower value of s
    s = s - 1;
    p = m - 4*s;
    CostMatrixPolEval1k_s_lower = s + 1 + ceil(p/s);
    
    % Loop to reduce 's' as much as possible while maintaining the same (or lower) cost.
    % This minimizes the memory required to store matrix powers (A^2, ..., A^s).
    while CostMatrixPolEval1 >= CostMatrixPolEval1k_s_lower  
        CostMatrixPolEval1 = CostMatrixPolEval1k_s_lower;
        s = s - 1;
        p = m - 4*s;
        CostMatrixPolEval1k_s_lower = s + 1 + ceil(p/s);
    end
    
    % Restore the optimal s after the while loop breaks
    s = s + 1;
    p = m - 4*s;
else
    % If s is provided as an input argument
    if s < 2
        error('MatrixPolEval:STooSmall', ...
              'The parameter s must be an integer greater than or equal to 2.');
    end
    
    p = m - 4*s;
    if p < 0
        error('MatrixPolEval:InvalidS', ...
              'The condition 4s <= m must be satisfied.');
    end
end

% Final cost of our evaluation method (number of matrix products)
CostMatrixPolEval1 = s + 1 + ceil(p/s);

% Optimal cost for the standard Paterson-Stockmeyer (PS) method
s_ps_ceil = ceil(sqrt(m));
s_ps_floor = floor(sqrt(m));
CostPS = min(ceil(m / s_ps_ceil) + s_ps_ceil, ceil(m / s_ps_floor) + s_ps_floor) - 2;

% Savings in the number of matrix products compared to standard PS
savings = CostPS - CostMatrixPolEval1;

% Set default values according to the chosen precision
switch lower(precision)
    case 'single'
        ndigits = max(16, ndigits);
        maxreal = realmax('single');
        tol = eps('single') / 2;
    case 'double'
        ndigits = max(32, ndigits);
        maxreal = realmax('double');
        tol = eps('double') / 2;
    otherwise
        error('MatrixPolEval:InvalidPrecision', ...
              'Accuracy parameter not valid (choose single or double)');
end

% Set precision to be used for Variable-Precision Arithmetic (VPA)
digits(ndigits);

% Create symbolic coefficient
syms A
c_prec = sym('c',[1 m+1]);

% Evaluate formulas 
y0s=A^s*sum(c_prec(m+2-s:m+1).*A.^(1:s));
fprintf('\nPOLYNOMIAL EVALUATION FORMULAS, polynomial degree m = 4s + p = %i, s = %i, p = %i\n\n',m,s,p)
fprintf('To save one matrix product with respect to the Paterson-Stockmeyer method, follow these steps:\n')
fprintf('    1. Compute and store the matrix powers A^2, A^3,...A^s with s = %i. These will be reused in subsequent computations.\n',s)
fprintf('    2. Evaluate the following formulas using the precomputed matrix powers:\n\n')
fprintf('y0s=A^s*(c(m+1)*A^s + c(m)*A^(s-1) + ... + c(m+1-s+1)*A);\n\n')
if type_pol==1
    y1s=(y0s+sum(c_prec(m+2-2*s:m+1-s).*A.^(1:s)))*(y0s+sum(c_prec(m+3-3*s:m+1-2*s).*A.^(2:s)))+c_prec(m+2-3*s)*y0s+sum(c_prec(m+1-4*s:m+1-3*s).*A.^(0:s));
    if leading_coeff_sign==1 
        fprintf('y1s=(y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
        fprintf('    (y0s + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+2)*A^2) + ...\n')
        fprintf('    c(m+1-3*s+1)*y0s + c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
    else
        fprintf('y1s=-((y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
        fprintf('    (y0s + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+2)*A^2) + ...\n')
        fprintf('    c(m+1-3*s+1)*y0s + c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I);\n\n')
    end
elseif type_pol==2 
    y1s=(y0s+sum(c_prec(m+1-2*s:m+1-s).*A.^(0:s)))*(y0s+sum(c_prec(m+2-3*s:m-2*s).*A.^(2:s)))+sum(c_prec(m+1-4*s:m+1-3*s).*A.^(0:s));
    if leading_coeff_sign==1
        fprintf('y1s=(y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A + c(m+1-2*s)*I) * ...\n')
        fprintf('    (y0s + c(m-2*s)*A^s + c(m-2*s-1)*A^(s-1) + ... + c(m+1-3*s+1)*A^2) + ...\n')
        fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
    else
        fprintf('y1s=-((y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A + c(m+1-2*s)*I) * ...\n')
        fprintf('    (y0s + c(m-2*s)*A^s + c(m-2*s-1)*A^(s-1) + ... + c(m+1-3*s+1)*A^2) + ...\n')
        fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I);\n\n')
    end
elseif type_pol==3
    y1s=(y0s+sum(c_prec(m+2-2*s:m+1-s).*A.^(1:s)))*(y0s+sum(c_prec(m+2-3*s:m+1-2*s).*A.^(1:s)))+sum(c_prec(m+1-4*s:m+1-3*s).*A.^(0:s));
if leading_coeff_sign==1
        fprintf('y1s=(y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
        fprintf('    (y0s + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+1)*A) + ...\n')
        fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
    else
        fprintf('y1s=-((y0s + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
        fprintf('    (y0s + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+1)*A) + ...\n')
        fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I);\n\n')
    end
end
% Final trailing polynomial evaluation
y1s=y1s*A^p+sum(c_prec(1:p).*A.^(0:p-1));


fprintf('y1s = y1s*A^p + c(p)*A^(p-1) + c(p-1)*A^(p-2) + ... + c(1)*I;\n')
fprintf('where p = %i, and the trailing polynomial of degree p-1 is evaluated using\nthe Paterson-Stockmeyer method with k = floor(p/s) = %i and r = mod(p,s) = %i, resulting in:\n\n', p, floor(p/s), mod(p,s))

% Formatted Paterson-Stockmeyer breakdown
fprintf('    y1s = (...((y1s*A^r + c(p)*A^(r-1) + c(p-1)*A^(r-2) + ... + c(p-r+1)*I) * A^s ...\n')
fprintf('          + c(p-r)*A^(s-1) + ... + c(p-r-s+1)*I) * A^s ...\n')
fprintf('          ...\n')
fprintf('          + c(2*s)*A^(s-1) + ... + c(s+1)*I) * A^s ...\n')
fprintf('          + c(s)*A^(s-1) + ... + c(1)*I)\n\n')

% Generate non-linear system by equating the coefficients of the matrix powers
[cy1s,~]=coeffs(y1s,A);
system=cy1s-b(end:-1:1);

all_cvpa=struct2cell(vpasolve(system,c_prec));
all_cvpa=[all_cvpa{:}];

% Round solution to the chosen accuracy
switch precision
    case 'single'
        csolved_prec=single(all_cvpa);
    case 'double'
        csolved_prec=double(all_cvpa);
end

% Check stability
nsolutions=size(csolved_prec,1);
relative_error_norm=zeros(1,nsolutions);
for i=1:nsolutions
    % Just select real solutions whose leading coefficients are positive
    if isreal(csolved_prec(i,:)) && csolved_prec(i,m+1)>0
        indices_zeros=find(b(end:-1:1)==0);
        if isempty(indices_zeros)
            relative_error=subs(system./b(end:-1:1),c_prec,csolved_prec(i,:));
        else
            absolute_error=subs(system,c_prec,csolved_prec(i,:));
            relative_error=subs(system./b(end:-1:1),c_prec,csolved_prec(i,:));
            relative_error(indices_zeros)=absolute_error(indices_zeros);
        end
        switch precision
            case 'single'
                relative_error=single(abs(relative_error));
            case 'double'
                relative_error=double(abs(relative_error));
        end
        relative_error_norm(i)=norm(relative_error,inf);
    else % Reject complex solutions and real solutions whose leading coefficients are negative
        relative_error_norm(i)=maxreal;
    end
end

% Select the best solution
[er_min,icsol]=min(relative_error_norm);
if isempty(er_min)
    fprintf('No solutions were found\n');
    c_vpa=[];
    c_prec=[];
    er_min=[];
elseif er_min==maxreal
    fprintf('No real solutions were found\n');
    c_vpa=[];
    c_prec=[];
else
    all_cvpa=flip(all_cvpa,2);
    c_vpa=all_cvpa(icsol,:);
    c_prec=flip(csolved_prec(icsol,:));
    if er_min > 10*tol
        % Inaccurate case
        fprintf('Warning: The evaluation formulas are likely to be inaccurate:\n')
        fprintf('    The best solution has a relative error/tol ratio = %.3e\n', er_min/tol);
        fprintf('    Relative error: %.3e (Threshold 10*tol: %.3e)\n\n', er_min, 10*tol);
    else
        % Accurate case
        fprintf('The evaluation formulas are likely to be accurate:\n')
        fprintf('    The best solution has a relative error %.3e <= 10*tol (%.3e)\n', er_min, 10*tol);
        fprintf('    Relative error/tol ratio = %.3e\n\n', er_min/tol);
    end
end
end
