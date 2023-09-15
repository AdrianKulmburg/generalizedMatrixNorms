n_list = 1:7;


parfor i=1:size(n_list,2)
    n = n_list(i);
    run_tests(n,2*n,2*n);
end

function [x,fval] = run_tests(n,m1,m2)

filename = [num2str(n),'_',num2str(m1),'_',num2str(m2),'.mat'];

N = n * (m1 + m2);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 1000); 

[x, fval, exitflag, output, trials] = surrogateopt(@(X) ratio(X), -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X)

X = reshape(X, [n m1+m2]);
G = X(:,1:m1);
H = X(:,m1+1:end);

if rank(H) < n
    r = 1;
    return
end

res_st = nesterov_relaxation(G, H);
res_exact = exact(G, H);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end

function res = nesterov_relaxation(G, H)
    A = G';
    B = H';
    
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


function res = exact(G1, G2)
G = G1;
norm_nu = @(nu) ellipsotopeNorm(G2, G*nu, 2);

% Number of generators of Z1
m = size(G, 2);

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

function res = ellipsotopeNorm(G,x,p)
    p_norm = @(beta) norm(beta,p);
    options = optimoptions('fmincon','Display', 'none');
    [~,res] = fmincon(p_norm, pinv(G)*x,[],[],G,x,[],[],[],options);
end


