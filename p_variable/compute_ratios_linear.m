n_list = 1:6;

p = 1.5;
%p = 2;

yalmip('clear')

parfor i=1:size(n_list,2)
    n = n_list(i);
    rng(123456);
    run_tests(n,2*n,2*n,p);
end

function [x,fval] = run_tests(n,ell,m,p)

filename = ['data/',num2str(p),'_',num2str(n),'_',num2str(ell),'_',num2str(m),'_linear.mat'];

N = n * (ell + m);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 1000); 

[x, fval, exitflag, output, trials] = surrogateopt(@(X) ratio(X,p), -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X,p)

X = reshape(X, [n ell+m]);
A_trans = X(:,1:ell);
B_trans = X(:,ell+1:end);

if rank(B_trans) < n
    r = 1;
    return
end

p_star = p/(p-1);

res_st = linear_relaxation(A_trans, B_trans, p_star);
res_exact = exact(A_trans, B_trans, p_star);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end

function res = linear_relaxation(A_trans, B_trans, p)

    function res_norm = L_1_p_T_norm(X)
        cost = 0;
        for i = 1:size(X,1)
            factor = 0;
            for j = 1:size(X,2)
                factor = factor + abs(X(i,j));
            end
            cost = cost + factor^p;
        end
        res_norm = cost^(1/p);
    end

    options = optimoptions('fmincon','Display','none');
    B_trans_plus = pinv(B_trans);
    F = eye(size(B_trans,2)) - B_trans_plus * B_trans;
    
    W0 = zeros(size(F));

    cost_function = @(W) L_1_p_T_norm(B_trans_plus * A_trans + F * reshape(W, size(F)));
    
    [~,res] = fmincon(cost_function, W0,[],[],[],[],[],[],[],options);
end


function res = exact(A_trans, B_trans, p)
norm_nu = @(nu) ellipsotopeNorm(B_trans, A_trans*nu, p);

% Number of generators of Z1
m = size(A_trans, 2);

% Create list of all combinations of generators we have to check (i.e., the
% choices of the +- signs from above). Note that this is a better strategy
% than computing the vertices directly, since it takes requires less
% memory.
% The next two lines produce all m-combinations of +-1
combinations = dec2bin(0:2^m-1)-'0';
combinations = 2*(combinations - 0.5);

res = 0;
for iter = combinations'
    res = max(res, norm_nu(iter));
end
end

function res = ellipsotopeNorm(B_trans,x,p)
    p_norm = @(beta) norm(beta,p);
    options = optimoptions('fmincon','Display', 'none');
    [~,res] = fmincon(p_norm, pinv(B_trans)*x,[],[],B_trans,x,[],[],[],options);
end


