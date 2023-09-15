p = 1.5;
%p = 2;
%p = 10;


n_list = 1:7;


parfor i=1:size(n_list,2)
    n = n_list(i);
    run_tests(n,2*n,2*n,p);
end

function [x,fval] = run_tests(n,m1,m2,p)

filename = [num2str(n),'_',num2str(m1),'_',num2str(m2),'_',num2str(p),'.mat'];

N = n * (m1 + m2);
Area = 10;


% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'iter',...
    'PlotFcn', [],...
    'MaxFunctionEvaluations', 1000); 

[x, fval, exitflag, output, trials] = surrogateopt(@(X) ratio(X,p), -Area.*ones([N 1])', Area.*ones([N 1])', [], options);

save(filename);

function r = ratio(X,p)

X = reshape(X, [n m1+m2]);
G = X(:,1:m1);
H = X(:,m1+1:end);

if rank(H) < n
    r = 1;
    return
end

res_st = linear_relaxation(G, H, p);
res_exact = exact(G, H, p);

if res_st >100*eps
    r = res_exact/res_st;
else
    r = res_exact/(100*eps);
end
end
end

function res = linear_relaxation(G, H, p)

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
    H_plus = pinv(H);
    F = eye(size(H,2)) - H_plus * H;
    
    W0 = zeros(size(F));

    cost_function = @(W) L_1_p_T_norm(H_plus * G + F * reshape(W, size(F)));
    
    [~,res] = fmincon(cost_function, W0,[],[],[],[],[],[],[],options);

end

function res = exact(G1, G2, p)
G = G1;
norm_nu = @(nu) ellipsotopeNorm(G2, G*nu, p);

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


