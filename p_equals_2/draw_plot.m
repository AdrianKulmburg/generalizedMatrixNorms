
n_range = 1:6;
m_range = 2:2:12;

clf;
hold on

for m = m_range
    worst_case = 0;
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

gamma_r = @(r) 2^(r/2)/sqrt(pi) * gamma((r+1)/2);

x = linspace(2, 12, 1000);
theory = plot(x, 0*x+sqrt(pi/2), 'k');

title("Worst-case approximation ratio for $q=$1, $p=2$, Nesterov Algorithm", 'Interpreter', 'latex', 'FontSize', 13)
xlabel("$\ell$ (in this case, $\ell=m$, $2n = \ell$)", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Approximation Ratio", 'Interpreter', 'latex', 'FontSize', 13)
ylim([1 1.4])

lgd = legend([data_points theory], {'$\rho_{\max}$', '$\sqrt{\pi/2}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
lgd.FontSize = 15;

%matlab2tikz('p_equals_2_nesterov.tex')