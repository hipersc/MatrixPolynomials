function [c_vpa,c_double,c_single,type_pol,er_min,erd_min,ers_min,savings,s,q]=MatrixPolEval2(pol,b,s,ndigits)
%
%--------------------------------------------------------------------------
% Description
%--------------------------------------------------------------------------
% MatrixPolEval2 computes the optimal coefficients involved in the 6s 
% formulation combined with the Paterson–Stockmeyer (PS) method. It 
% displays the general evaluation formulas and computes their coefficients,
% saving up to two matrix products (2M) with respect to PS method for 
% polynomials with a non-zero leading coefficient [1, 2, 3].
%
% [c_vpa,c_double,c_single,type_pol,er_min,erd_min,ers_min,savings,s,p]=...
%    MatrixPolEval2(pol,b,s,ndigits)
%
% This implementation considers evaluation formulas Y0, Y1, Y2 and Z2qs 
% (3)-(6) and (15) from [1, Sec. 2], respectively. 2M savings are possible 
% for degrees m = 18, 21, 24, 26, 27, 28 and m >= 30. 1M savings are 
% possible for m = 12, 13, 14, 19, 20, 22, 23, 25 and 29. Out cases, the PS
% method is recommended. The function filters out complex solutions and 
% real solutions with negative leading coefficients, as similar solutions 
% with positive leading coefficients are available. A warning is issued if 
% the numerical stability is low, based on the test in [3, Ex. 3.1].
%
%--------------------------------------------------------------------------
% Input arguments
%--------------------------------------------------------------------------
% - pol: string identifying the polynomial b (default: 'exp', i.e. 
%   exponential).
% - b: vector [b(1) b(2) ... b(m+1)] containing the coefficients of a 
%   polynomial b(A) of degree m, b(1)*I + b(2)*A + b(3)*A^2 + ... +
%   b(m+1)*A^m. The degree of the polynomial must be greater than or equal 
%   to 12 (default: b = 1./factorial(sym(0:12)); degree-12 Taylor 
%   approximation).
% - s: value s used in the formulation 6s [1, Sec.2], such that m-6s>= 0.
%   (default: minimal s giving the same cost to minimize memory overhead).
% - ndigits: number of significant digits used in the MATLAB variable 
%   precision arithmetic (vpa) function (default: 32).
%
%--------------------------------------------------------------------------
% Output arguments
%--------------------------------------------------------------------------
% - c_vpa, c_double, c_single: coefficients of the 6s + Paterson–Stockmeyer
%   formulation [1] computed in vpa, double, and single precision, 
%   respectively.
%   c_vpa = [c(1) ... c(q)| c(q+1) ... c(m+1)], where:
%    - c(1)...c(q) are the coefficients used to evaluate the Paterson–
%      Stockmeyer part of the formula.
%    - c(q+1)...c(m+1) are the coefficients corresponding to the 6s 
%      formulation.
% - type_pol: variant of the 4s formulation used (1, 2, or 3 [1, Sec. 2]) 
%   that yields the minimum error according to the stability test in 
%   [3, Ex. 3.1].
% - er_min, erd_min, and ers_min: they correspond, respectively, to the 
%   minimum value of the reconstruction maximum relative error [2, Ex. 3.1]
%   in vpa (er_min), double precision (erd_min), or single precision 
%   (ers_min). For a reliable polynomial evaluation, the target stability 
%   thresholds are:
%   - erd_min < 10 * 2^-53 (approx. 1.11e-15)
%   - ers_min < 10 * 2^-24 (approx. 5.96e-07)
%   Lower values indicate superior numerical accuracy.
% - savings: number of matrix products saved compared to the PS method.
% - s: value s finally used in the 6s formulation. If the argument s is not 
%   used as an input value, then it is calculated in such a way that, with 
%   the decomposition m=6s+q, the lowest possible cost is obtained, with 
%   the smallest s and the largest q.
% - q: corresponds to the number of coefficients in the residual polynomial 
%   block evaluated with the Paterson-Stockmeyer method.
%
% The output variables are also saved to a file (e.g., 'MatrixPolEval2_exp_
% m12_s2_q0_d32_f1.mat'), where 'exp' denotes the polynomial name, m, s, 
% and q are the corresponding parameters, d is the number of VPA digits 
% used, and f is the formulation number associated with type_pol.
%
%--------------------------------------------------------------------------
% Usage examples
%--------------------------------------------------------------------------
% Example 1: In this example, we show the procedure required to obtain the 
% elements corresponding to the degree-12 polynomial whose coefficients 
% are all equal to 1. Then, they will be used to compute the matrix 
% polynomial of a given matrix. In this case, the function uses the values 
% s = 2 and ndigits = 32 by default.
%
% b = ones(13,1);
% [c_vpa,c_double,c_single,type_pol,er_min,erd_min,ers_min,savings,s,q]=...
%    MatrixPolEval2('exp',b)
%
% After executing this instruction, the coefficients are obtained in vpa, 
% double, and single precision, together with the following output 
% arguments:
% type_pol = 2, er_min = 5.7016e-34, erd_min = 1.2969e-16, 
% ers_min = 8.5798e-08, savings = 1, s = 2, and q = 0.
%
% To compute the value of the above polynomial for the matrix 
% A = [0.1 -0.2; -0.1 0.3], we use the function z2qs as shown below:
%
% y2 = z2qs([0.1 -0.2; -0.1 0.3], c_double, s, type_pol)
%
% y2 =
%     1.1475   -0.3279
%    -0.1639    1.4754
%
% which coincides with the matrix obtained using the function polyvalm:
%
% B = polyvalm(ones(13,1), [0.1 -0.2; -0.1 0.3])
%
% B =
%     1.1475   -0.3279
%    -0.1639    1.4754
%
% Example 2: In this example, we show how to obtain the coefficients of the 
% Z2qs formulation corresponding to the degree-35 Taylor polynomial of the 
% exponential function with s = 5 and ndigits = 64. Then, they will used to
% compute the Taylor approximation of the matrix. Note that scaling and 
% squaring strategies may need to be employed to achieve higher accuracy.
%
% b = 1./factorial(sym(0:35));
% [c_vpa,c_double,c_single,type_pol,er_min,erd_min,ers_min,savings,s,q]=...
%    MatrixPolEval2('exp',b,5,64)
%
% In this case, we obtain:
% type_pol = 1, er_min = 5.6970e-65, erd_min = 1.3333e-16,
% ers_min = 8.2240e-08, savings = 2, s = 5, and q = 5.
%
% To compute the value of the above polynomial for the matrix 
% A = [1 2; 0 3], we use the function z2qs as shown below:
%
% y2 = z2qs([1 2; 0 3], c_double, s, type_pol)
%
% y2 =
%     2.7183   17.3673
%     0       20.0855
%
% Since 
% E = [exp(vpa(1))  exp(vpa(3)) - exp(vpa(1)); 
%      0            exp(vpa(3))] 
% is equal to exp(A), the relative error can be computed as:
%
% Er=double(norm(E - y2) / norm(E))
%
% Er =
%     3.025514947427634e-17
%
%--------------------------------------------------------------------------
% Copyright
%--------------------------------------------------------------------------
% Authors: Javier Ibáñez, Jorge Sastre and José Miguel Alonso
% Revised version: August 3, 2026.
%
% High Performance Scientific Computing Group (HiPerSC)
% Universitat Politècnica de València (Spain)
%
%--------------------------------------------------------------------------
% References
%--------------------------------------------------------------------------
% [1] J. Ibáñez, J. Sastre, J.M. Alonso, E. Defez, A MATLAB Tool for the 
%     Stable Generation of Matrix Polynomial Evaluation Schemes with 
%     Two-Product Savings, arXiv:2607.28286, 2026. 
%     https://arxiv.org/abs/2607.28286
% [2] J. Sastre, J. Ibáñez, Efficient Evaluation of Matrix Polynomials 
%     Beyond the Paterson-Stockmeyer Method, Mathematics 9(14), 1600, 2021,
%     https://doi.org/10.3390/math9141600
% [3] J. Sastre, Efficient evaluation of matrix polynomials, Linear 
%     Algebra Applications, 539, 2018, 229-250.
%     https://doi.org/10.1016/j.laa.2017.11.010
%
tic
st = dbstack;
name= st.name;

c_vpa=[];c_double=[];c_single=[];type_pol=[];er_min=[];savings=[];
% ==============================
% Default values handling
% ==============================

% If pol is not specified
if nargin < 1 || isempty(pol)
    pol = 'exp';
end

% If b is not specified
if nargin < 2 || isempty(b)
    b = 1./factorial(sym(0:12));
elseif ~isrow(b) % b must be a row vector
    b = b(:).';   
end

if b(end) == 0
    error('MatrixPolEval:ZeroLeadingCoefficient', ...
          'The leading coefficient of the polynomial must be non-zero.');
end

m=length(b);
m=m-1;

% If s is not specified
if nargin < 3 || isempty(s)
    if m<12
        error('MatrixPolEval:DegreeNotSupported',...
            'The polynomial degree must be greater than or equal to 12.')
    end
    if m==15 || m==16 || m==17
        warning('MatrixPolEval:DegreeNotSupported', ...
            'Matrix product savings are not available for the specified polynomial degree.\nPlease apply the Paterson-Stockmeyer method directly.')
    end
    if m<36
        s=floor(m/6);
    else
        s=floor(sqrt(m));
    end
q = m - 6*s;

% Selecting the smallest s that yields the same computational cost,
% in order to minimize memory overhead.
CostMatrixPolEval2 = s + 2 + ceil(q/s);
s = s - 1;
q = m - 6*s;
CostMatrixPolEval1k_s_lower = s + 2 + ceil(q/s);
while (CostMatrixPolEval2 >= CostMatrixPolEval1k_s_lower)&&(s>2)
    CostMatrixPolEval2 = CostMatrixPolEval1k_s_lower;
    s = s - 1;
    q = m - 6*s;
    CostMatrixPolEval1k_s_lower = s + 2 + ceil(q/s);
end
s=s+1;
end
q=m-6*s;

% Cost for the proposed method
CostMatrixPolEval2 = s + 2 + ceil(q/s);

% Optimal cost for the standard Paterson-Stockmeyer (PS) method
s_ps_ceil = ceil(sqrt(m));
s_ps_floor = floor(sqrt(m));
CostPS = min(ceil(m / s_ps_ceil) + s_ps_ceil, ceil(m / s_ps_floor) + s_ps_floor) - 2;

% Savings in the number of matrix products compared to standard PS
savings = CostPS - CostMatrixPolEval2;

% If ndigits is not specified
if nargin < 4 || isempty(ndigits)
    ndigits = 32;
end

% ==============================
% Set symbolic precision
% ==============================
digits(ndigits);

bq=b(q+2:end);
fprintf('\nDetermination of the coefficients of the 6s formulation for m=%d, s = %d, q = %d\n',m,s,q)

fprintf('\n1: Computation of the coefficients of the 6s formulation\n')
[csol,~]= coef_6s(bq);
t6s=toc;
fprintf('\tTotal time for the 6s formulation: %f seconds\n',t6s)

tic
ns6=length(csol);
er_min=realmax;

fprintf('\n2: Computation of the 4s coefficients for the three variants of the 4s formulation\n')
c_min=[];
for i=1:ns6
    pm=csol{i}(2*s+1:7*s);
    for f=1:3
        [csol4,er]=coef_4s(pm,f);
        if er>1
            continue
        end
        if ~isempty(csol4)
            ns4=length(csol4);
            for j=1:ns4
                c=[csol{i}(1:2*s),csol4{j},csol{i}(6*s+1:7*s)];
                [er,erd,ers]= test(bq,c,f);
                if er<er_min
                    er_min=er;
                    erd_min=erd;
                    ers_min=ers;
                    type_pol=f;
                    c_min=c;
                end
            end
        end
    end
end

c_vpa=[b(1:q+1),c_min];
c_double=double(c_vpa);
c_single=single(c_vpa);

filename = sprintf('%s_%s_m%d_s%d_q%d_d%d_f%d',name,pol,m,s,q,digits,type_pol);
save(filename,'c_vpa','c_double','c_single','type_pol','er_min','erd_min','ers_min','savings','q','s')

t4s=toc;

fprintf('\tTotal time for the 4s coefficient computation: %f seconds\n\n',t4s)
fprintf('\tTotal execution time (6s + 4s): %f seconds. Errors: vpa (%g), double (%g), single (%g)\n',...
        t4s+t6s,er_min,erd_min,ers_min)

fprintf('\nPOLYNOMIAL Pm EVALUATION FORMULAS, polynomial degree m = 6s + q = %i, s = %i, q = %i\n\n',m,s,q)
fprintf('To save up to 2 matrix products with respect to the Paterson-Stockmeyer method, follow these steps:\n')
fprintf('    1. Compute and store the matrix powers A^2, A^3,...A^s with s = %i. These will be reused in subsequent computations.\n',s)
fprintf('    2. Evaluate the following formulas using the precomputed matrix powers:\n\n')
fprintf('y0 = A^s*(c(m+1)*A^s + c(m)*A^(s-1) + ... + c(m+1-s+1)*A);\n\n')


if type_pol==1
    fprintf('y1 = (y0 + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
    fprintf('    (y0 + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+2)*A^2) + ...\n')
    fprintf('    c(m+1-3*s+1)*y0 + c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
elseif type_pol==2
    fprintf('y1 = (y0 + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A + c(m+1-2*s)*I) * ...\n')
    fprintf('    (y0 + c(m-2*s)*A^s + c(m-2*s-1)*A^(s-1) + ... + c(m+1-3*s+1)*A^2) + ...\n')
    fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
elseif type_pol==3
    fprintf('y1 = (y0 + c(m+1-s)*A^s + c(m-s)*A^(s-1) + ... + c(m+1-2*s+1)*A) * ...\n')
    fprintf('    (y0 + c(m+1-2*s)*A^s + c(m-2*s)*A^(s-1) + ... + c(m+1-3*s+1)*A) + ...\n')
    fprintf('    c(m+1-3*s)*A^s + c(m-3*s)*A^(s-1) + ... + c(m+1-4*s+1)*A + c(m+1-4*s)*I;\n\n')
end

fprintf('y2 = y1 * (y0 + c(m+1-4*s)*A^s + c(m-2*4)*A^(s-1) + ... + c(m+1-5*s+1)*A) + ...\n')
fprintf('    c(m+1-5*s)*A^s + c(m-5*s)*A^(s-1) + ... + c(m+1-6*s+1)*A + c(m+1-6*s)*I;\n\n')


% Final trailing polynomial evaluation
if q > 0
    k = floor(q/s);
    r = mod(q,s);
    fprintf('pm = y2*A^q + c(q)*A^(q-1) + c(q-1)*A^(q-2) + ... + c(1)*I;\n')
    fprintf('where q = %i, and the trailing polynomial is evaluated using\nthe Paterson-Stockmeyer method with k = floor(q/s) = %i and r = mod(q,s) = %i, resulting in:\n\n', q, k, r)
    if r~=0
        % Formatted Paterson-Stockmeyer breakdown (Simplified for the user)
        fprintf('    pm = (...((y2*A^r + c(q)*A^(r-1) + c(q-1)*A^(r-2) + ... + c(q-r+1)*I) ...\n')
        fprintf('          * A^s + c(q-r)*A^(s-1) + ... + c(q-r-s+1)*I) ...\n')
    else
        fprintf('    pm = (...((y2*A^s + c(q)*A^(s-1) + c(q-1)*A^(s-2) + ... + c(q-s+1)*I) ...\n')
        fprintf('          * A^s + c(q-s)*A^(s-1) + ... + c(q-2*s+1)*I) ...\n')
    end
    fprintf('          ...\n')
    fprintf('          * A^s + c(2*s)*A^(s-1) + ... + c(s+1)*I) ...\n')
    fprintf('          * A^s + c(s)*A^(s-1) + ... + c(1)*I\n\n')
else
% Direct assignment when q = 0
    fprintf('Final evaluation (m = 6s, q = 0):\n')
    fprintf('No trailing polynomial evaluation is required:\n')    
    fprintf('pm = y2\n\n')
end

u=eps('single')/2;
if ers_min > 10*u
    % Inaccurate case
    fprintf('\nWarning: The evaluation formulas are likely to be inaccurate for single precision:\n')
    fprintf('    The best solution has a relative error/u ratio = %.3e\n', ers_min/u);
    fprintf('    Relative error: %.3e (Threshold 10*u: %.3e)\n\n', ers_min, 10*u);
else
    % Accurate case
    fprintf('\nThe evaluation formulas are likely to be accurate for single precision:\n')
    fprintf('    The best solution has a relative error %.3e <= 10*u (%.3e)\n', ers_min, 10*u);
    fprintf('    Relative error/u ratio = %.3e\n\n', ers_min/u);
end

u=eps('double')/2;
if erd_min > 10*u
    % Inaccurate case
    fprintf('\nWarning: The evaluation formulas are likely to be inaccurate for double precision:\n')
    fprintf('    The best solution has a relative error/u ratio = %.3e\n', erd_min/u);
    fprintf('    Relative error: %.3e (Threshold 10*u: %.3e)\n\n', erd_min, 10*u);
else
    % Accurate case
    fprintf('\nThe evaluation formulas are likely to be accurate for double precision:\n')
    fprintf('    The best solution has a relative error %.3e <= 10*u (%.3e)\n', erd_min, 10*u);
    fprintf('    Relative error/u ratio = %.3e\n\n', erd_min/u);
end
end


% =====================================================================
% ======================== 6s FORMULATION ==============================
% =====================================================================

function [csol,er]= coef_6s(p)
% Computes the coefficients of the 6s formulation.

s=length(p)/6;
syms A 

for i=1:7*s
    eval(['c',int2str(i),'=sym(''c',int2str(i),''');']);
    eval(['c(i)=c',int2str(i),';']);
end

% Creation of the vectors y0, y1, and y2
y0=sum(c(6*s+1:7*s).*A.^(s+1:2*s));
y1=sum(c(2*s+1:6*s).*A.^(1:4*s));
[cy1_y02,~]=coeffs(y1-y0*y0,A);

y2=y1*(y0+sum(c(s+1:2*s).*A.^(1:s)))+sum(c(1:s).*A.^(1:s));
[cy2{1},~]=coeffs(y2,A);
cy2{1}=cy2{1}(end:-1:1)-p;
cy2m=cy2{1};

fprintf('\tCreation of the vectors y0, y1, and y2: %f seconds\n',toc)

% Step 1: Solve for c(6s+1:7s) in terms of c(5s+1:6s)
sol=solve(cy1_y02(1:s),c(6*s+1:7*s));
sol=struct2cell(sol);
sol=[sol{:}];

[n1,~]=size(sol);
aux=cy2{1};

for j=1:n1
    cy2{j}=aux;
    cs{j}(6*s+1:7*s)=sol(j,1:s);
    cy2{j}=subs(cy2{j},c(6*s+1:7*s),sol(j,1:s));
end

fprintf('\tSolving for c(6s+1:7s): %f seconds\n',toc)

% Elimination of the unknowns c(2s+1:6s)

j=1;
while j<=n1
    i=6*s;
    sol=solve(cy2{j}(i),c(i));
    [nsol,~]=size(sol);

    aux_cs=cs;
    aux_cy2=cy2; 
    l=0;

    for k=1:nsol
        if isreal(sol(k))
            cs{j}=aux_cs{j};
            cs{j}(i)=sol(k);
            cy2{j}=subs(aux_cy2{j},c(i),sol(k));
            l=l+1;
        end
    end

    if l>0
        for i=6*s-1:-1:2*s+1
            sol=solve(cy2{j}(i),c(i));
            cs{j}(i)=sol(1);
            cy2{j}=simplify(subs(cy2{j},c(i),sol(1)));
        end
        j=j+1;
    else
        n1=n1-1;
        for k=j:n1
            cs{k}=aux_cs{k+1};
            cy2{k}=aux_cy2{k+1};
        end
    end
end

fprintf('\tElimination of the unknowns c(2s+1:6s): %f seconds\n',toc)

% Numerical solution for c(s+1:2s)

k=1;
while k<=n1
    sol=vpasolve(cy2{k}(s+1:2*s),c(s+1:2*s),'Random',true);
    sol=struct2cell(sol);
    sol=[sol{:}];

    [nsol,~]=size(sol);
    h=[];

    for i=1:nsol
        if isreal(sol(i,:))
            h=[h,i];
        end
    end

    if ~isempty(h)
        n2=length(h);
        for j=1:n2
            cs{k,j}=cs{k};
            cy2{k,j}=cy2{k};
        end

        for j=1:n2
            for i=1:s
                cs{k,j}(s+i)=sol(h(j),i);
                cy2{k,j}=subs(cy2{k,j},c(s+i),sol(h(j),i));
            end
        end
        k=k+1;
    else
        n1=n1-1;
        for j=k:n1
            cs{j}=cs{j+1};
            cy2{j}=cy2{j+1};
        end
    end
end

fprintf('\tNumerical computation of c(s+1:2s): %f seconds\n',toc)

% Compute remaining coefficients

for k1=1:n1
    for k2=1:n2
        for i=2*s+1:6*s 
            cs{k1,k2}(2*s+1:6*s)=subs(cs{k1,k2}(2*s+1:6*s),...
                c(s+1:2*s),cs{k1,k2}(s+1:2*s));
        end
    end
end

fprintf('\tNumerical computation of c(2s+1:6s): %f seconds\n',toc)

for k1=1:n1
    for k2=1:n2
        for i=6*s+1:7*s
            for j=5*s+1:6*s
                cs{k1,k2}(i)=subs(cs{k1,k2}(i),c(j),cs{k1,k2}(j));
            end
        end
    end
end

fprintf('\tNumerical computation of c(6s+1:7s): %f seconds\n',toc)

for k1=1:n1
    for k2=1:n2
        for i=1:s 
            sf=solve(cy2{k1,k2}(i),c(i));
            cs{k1,k2}(i)=sf;
            cy2{k1,k2}=subs(cy2{k1,k2},c(i),sf);
        end
    end
end

fprintf('\tFinal computation of c(1:s) and storage of coefficients: %f seconds\n',toc)

fprintf('\tErrors in the solutions obtained during the initial 6s computation:\n')

n=0;
for k1=1:n1
    for k2=1:n2
        n=n+1;
        csol{n}=vpa(cs{k1,k2});
        er{n}=vpa(norm(subs(cy2m./p,c,csol{n}),inf));
        fprintf('\t\tSolution %d: error = %g\n',n,er{n})
    end
end

end


% =====================================================================
% ======================== 4s FORMULATION ==============================
% =====================================================================

function [cs,er]=coef_4s(p,f)
% Auxiliary function that computes the coefficients of the 4s formulation.

s=length(p)/5;
syms A 

for i=1:4*s
    eval(['c',int2str(i),'=sym(''c',int2str(i),''');']);
    eval(['c(i)=c',int2str(i),';']);
end

y0=sum(c(3*s+1:4*s).*A.^(s+1:2*s));

if f==1
    y1=(y0+sum(c(2*s+1:3*s).*A.^(1:s)))* ...
       (y0+sum(c(s+2:2*s).*A.^(2:s)))+ ...
       c(s+1)*y0+sum(c(1:s).*A.^(1:s));
elseif f==2
    y1=(y0+sum(c(2*s:3*s).*A.^(0:s)))* ...
       (y0+sum(c(s+1:2*s-1).*A.^(2:s)))+ ...
       sum(c(1:s).*A.^(1:s));
else
    y1=(y0+sum(c(2*s+1:3*s).*A.^(1:s)))* ...
       (y0+sum(c(s+1:2*s).*A.^(1:s)))+ ...
       sum(c(1:s).*A.^(1:s));
end

[cy1,~]=coeffs(y1,A);
cy1=cy1(end:-1:1)-p(1:4*s);
cy1=subs(cy1,c(3*s+1:4*s),p(4*s+1:5*s));
cy1m=cy1;

sol=vpasolve(cy1(1:3*s),c(1:3*s));
sol=struct2cell(sol);
sol=[sol{:}];

cs=[];er=[];

if ~isempty(sol)
    nsol=size(sol,1);
    k=0;
    for i=1:nsol
        if isreal(sol(i,:))
            k=k+1;
            cs{k}=sol(i,:);
            er=vpa(norm(subs(cy1m(1:3*s)./p(1:3*s),...
                c(1:3*s),cs{k}(1:3*s)),inf));
            if er>1
                fprintf('\tLarge error detected while solving expression %d (%f seconds). An alternative solution will be considered.\n',f,toc)
            end
        end
    end
end

end


% =====================================================================
% =========================== ERROR TEST ===============================
% =====================================================================

function [er,erd,ers]= test(b,cs,f)
% Auxiliary function that evaluates the error of the 6s–4s formulations.

s=length(cs)/6;
syms A 

for i=1:6*s
    eval(['c',int2str(i),'=sym(''c',int2str(i),''');']);
    eval(['c(i)=c',int2str(i),';']);
end

y0=sum(c(5*s+1:6*s).*A.^(s+1:2*s));

if f==1
    y1=(y0+sum(c(4*s+1:5*s).*A.^(1:s)))* ...
       (y0+sum(c(3*s+2:4*s).*A.^(2:s)))+ ...
       c(3*s+1)*y0+sum(c(2*s+1:3*s).*A.^(1:s));
elseif f==2
    y1=(y0+sum(c(4*s:5*s).*A.^(0:s)))* ...
       (y0+sum(c(3*s+1:4*s-1).*A.^(2:s)))+ ...
       sum(c(2*s+1:3*s).*A.^(1:s));
else
    y1=(y0+sum(c(4*s+1:5*s).*A.^(1:s)))* ...
       (y0+sum(c(3*s+1:4*s).*A.^(1:s)))+ ...
       sum(c(2*s+1:3*s).*A.^(1:s));
end

y2=y1*(y0+sum(c(s+1:2*s).*A.^(1:s)))+sum(c(1:s).*A.^(1:s));

[cy2,~]=coeffs(y2,A);
cy2=cy2(end:-1:1)-b;

er = double(norm(subs(cy2./b,c,cs),inf));
erd = double(norm(subs(cy2./double(b),c,double(cs)),inf));
ers = single(norm(subs(cy2./single(b),c,single(cs)),inf));

end