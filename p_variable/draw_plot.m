p = 1.5;
%p = 2;
%p = 10;


n_range = 1:7;
m_range = 2:2:14;

clf;
hold on

for m = m_range
    worst_case = 0;
    for n = n_range
        try
            filename = ['data/',num2str(n),'_',num2str(m),'_',num2str(m),'_',num2str(p),'.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points = plot(m,worst_case,'ok');
end

gamma_r = @(r) 2^(r/2)/sqrt(pi) * gamma((r+1)/2);

x = linspace(2, 14, 1000);
theory = plot(x, gamma_r(p)^(1/p)/gamma_r(1) * sqrt(x), 'k');

title("Worst-case approximation ratio for $q=$1, $p=$"+num2str(p), 'Interpreter', 'latex', 'FontSize', 13)
xlabel("$\ell$ (in this case, $\ell=m$, $2n = m$)", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Approximation Ratio", 'Interpreter', 'latex', 'FontSize', 13)

lgd = legend([data_points theory], {'Worst-Case Results', '$C_p\cdot\sqrt{\ell}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
lgd.FontSize = 15;

matlab2tikz(['p_equals_',num2str(p),'.tex'])