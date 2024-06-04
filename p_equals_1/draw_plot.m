n_range = 1:6;
ell_range = 2:20;

clf;
hold on

for ell = ell_range
    worst_case = 1;
    for n = n_range
        try
            filename = ['data/',num2str(n),'_',num2str(ell),'_',num2str(ell),'.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points = plot(ell,worst_case,'ok');
end

x = linspace(2, 20, 1000);
theory = plot(x, sqrt(x), 'k');

title("Worst-case approximation ratio for $q=1$, $p=1$", 'Interpreter', 'latex')
xlabel("$\ell$", 'Interpreter', 'latex')
ylabel("Approximation Ratio", 'Interpreter', 'latex')

lgd = legend([data_points theory], {'$\rho_{\max}$', '$C_1\cdot\sqrt{\ell}$'}, 'Interpreter', 'latex', 'Location', 'northwest');

%matlab2tikz('p_equals_1.tex')