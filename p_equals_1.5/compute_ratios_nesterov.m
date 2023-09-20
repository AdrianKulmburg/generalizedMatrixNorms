p = 1.5;

n_list = 1:6;


for i=1:size(n_list,2)
    n = n_list(i);
    rng(123456);
    run_tests(n,2*n,2*n,p);
end

function [x,fval] = run_tests(n,ell,m,p)

filename = [num2str(n),'_',num2str(ell),'_',num2str(m),'_',num2str(p),'_nesterov.mat'];

N = n * (ell + m);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 1000); 

[x, fval, exitflag, output, trials] = surrogateopt(@(X) ratio(X), -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X)

X = reshape(X, [n ell+m]);
A_trans = X(:,1:ell);
B_trans = X(:,ell+1:end);

if rank(B_trans) < n
    r = 1;
    return
end
p_star = p/(p-1);
res_st = nesterov_relaxation(A_trans, B_trans);
res_exact = exact(A_trans, B_trans, p_star);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end

function res = nesterov_relaxation(A_trans, B_trans)
    A = A_trans';
    B = B_trans';
    
    X = sdpvar(size(A,1), size(A,1));
    
    
    cost = -trace(A * pinv(B) * pinv(B)' * A' * X);
    constraints = [X >= 0];
    for i=1:size(A,1)
        constraints = [constraints X(i,i) == 1];
    end
    
    options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    
    yalmipOptimizer = optimizer(constraints,cost,options,[],{X});
    
    [sol, exitflag] = yalmipOptimizer();
    
    res = sqrt(abs(trace(A * pinv(B) * pinv(B)' * A' * sol)));
end


function res = exact(A_trans, B_trans, p)

norm_nu = @(nu) ellipsotopeNorm(B_trans, A_trans*nu, p);

% Number of generators
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



