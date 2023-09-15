n_list = 1:7;
m_list = 1:20;

for n=1:1
    parfor m=m_list
        if m<n
            continue
        end
        rng(123456);
        run_tests(n,m,m)
    end
end


function [x,fval] = run_tests(n,m1,m2)

filename = [num2str(n),'_',num2str(m1),'_',num2str(m2),'.mat'];

N = n * (m1 + m2);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 10000);


[x, fval, exitflag, output, trials] = surrogateopt(@ratio, -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X)

X = reshape(X, [n m1+m2]);
G = X(:,1:m1);
H = X(:,m1+1:end);

if rank(H) < n
    r = 1;
    return
end


res_st = linearRelaxation(G, H);
res_exact = zonotopeNormPolymax(G, H);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end


function d = linearRelaxation(G, H)

nx = size(G, 2);
ny = size(H, 2);

Gamma = sdpvar(ny, nx, 'full');
constraints = [G == H*Gamma];
cost = norm(Gamma, Inf); 
options = sdpsettings('solver','linprog', 'verbose',0, 'allownonconvex',0);

% solve linear programming problem
yalmipOptimizer = optimizer(constraints, cost, options,[],{Gamma});

[optimizationResults, ~] = yalmipOptimizer();

d = norm(optimizationResults, Inf);
end  


function d = polyhedralMaximization(G, H)
% Solves the zonotope containment problem by computing the maximal value of
% the polyhedral norm over Z1 w.r.t. Z2 (see also [2, Algorithm 2]).

n = size(H,1);
Z = zonotope(zeros([n 1]), H);
% Then, we compute the halfspace representation
Z = halfspace(Z);

% We then need the normalized halfspace-matrix
H_norm = Z.halfspace.H ./ Z.halfspace.K;

% Polyhedral norm at a point p
poly_norm = @(p) max([0 max(H_norm * p)]);

% Store signs for the main step (these are the signs that matter for the
% decision of the sign of the x_j in [2, Algorithm 2]).
M = sign(H_norm * G);

n = size(M, 1);

d = 0;
% Iterate over each row of M
for i = 1:n
    mu = M(i,:); % Sign-combination
    maximum = poly_norm(G*mu'); % Compute the resulting polyhedral norm
    d = max(d, maximum);
end
end
