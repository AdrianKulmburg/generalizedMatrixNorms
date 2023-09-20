n_list = 1:6;
ell_list = 1:20;

for n=n_list
    parfor ell=ell_list
        if ell<n
            continue
        end
        rng(123456);
        run_tests(n,ell,ell)
    end
end


function [x,fval] = run_tests(n,ell,m)

filename = [num2str(n),'_',num2str(ell),'_',num2str(m),'.mat'];

N = n * (ell + m);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 10000);


[x, fval, exitflag, output, trials] = surrogateopt(@ratio, -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X)

X = reshape(X, [n ell+m]);
A_trans = X(:,1:ell);
B_trans = X(:,ell+1:end);

if rank(B_trans) < n
    r = 1;
    return
end


res_st = linearRelaxation(A_trans, B_trans);
res_exact = polyhedralMaximization(A_trans, B_trans);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end


function d = linearRelaxation(A_trans, B_trans)

nx = size(A_trans, 2);
ny = size(B_trans, 2);

Gamma = sdpvar(ny, nx, 'full');
constraints = [A_trans == B_trans*Gamma];
cost = norm(Gamma, Inf); 
options = sdpsettings('solver','linprog', 'verbose',0, 'allownonconvex',0);

% solve linear programming problem
yalmipOptimizer = optimizer(constraints, cost, options,[],{Gamma});

[optimizationResults, ~] = yalmipOptimizer();

d = norm(optimizationResults, Inf);
end  


function d = polyhedralMaximization(A_trans, B_trans)
% Solves the zonotope containment problem by computing the maximal value of
% the polyhedral norm.

n = size(B_trans,1);
Z = zonotope(zeros([n 1]), B_trans);
% Then, we compute the halfspace representation
Z = halfspace(Z);

% We then need the normalized halfspace-matrix
H_norm = Z.halfspace.H ./ Z.halfspace.K;

% Polyhedral norm at a point p
poly_norm = @(p) max([0 max(H_norm * p)]);

% Store signs for the main step (these are the signs that matter for the
% decision of the sign of the x_j).
M = sign(H_norm * A_trans);

n = size(M, 1);

d = 0;
% Iterate over each row of M
for i = 1:n
    mu = M(i,:); % Sign-combination
    maximum = poly_norm(A_trans*mu'); % Compute the resulting polyhedral norm
    d = max(d, maximum);
end
end
