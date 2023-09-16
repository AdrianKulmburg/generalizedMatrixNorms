n_range = 1:7;
m_range = 2:20;

clf;
hold on

for m = m_range
    worst_case = 1;
    for n = n_range
        try
            filename = ['data/',num2str(n),'_',num2str(m),'_',num2str(m),'.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points = plot(m,worst_case,'ok');
end

x = linspace(2, 20, 1000);
theory = plot(x, sqrt(x), 'k');

title("Worst-case approximation ratio for $q=$1, $p=$1", 'Interpreter', 'latex', 'FontSize', 13)
xlabel("$\ell$", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Approximation Ratio", 'Interpreter', 'latex', 'FontSize', 13)

lgd = legend([data_points theory], {'$\rho_{\max}$', '$C_1\cdot\sqrt{\ell}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
lgd.FontSize = 15;

%matlab2tikz('p_equals_1.tex')