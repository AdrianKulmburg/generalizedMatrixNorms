p = 1.5;


n_range = 1:6;
ell_range = 2:2:12;

clf;
hold on

for ell = ell_range
    worst_case = 0;
    for n = n_range
        try
            filename = ['data/',num2str(n),'_',num2str(ell),'_',num2str(ell),'_',num2str(p),'.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points = plot(ell,worst_case,'ok');
end

for ell = ell_range
    worst_case = 0;
    for n = n_range
        try
            filename = ['data/',num2str(n),'_',num2str(ell),'_',num2str(ell),'_',num2str(p),'_nesterov.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points_nesterov = plot(ell,worst_case,'xk');
end

gamma_r = @(r) 2^(r/2)/sqrt(pi) * gamma((r+1)/2);

x = linspace(2, 12, 1000);
theory = plot(x, gamma_r(p)^(1/p)/gamma_r(1) * sqrt(x), 'k');
theory_nesterov = plot(x, sqrt(pi/2) * x.^(1/p-1/2), 'k--');

title("Worst-case approximation ratio for $q=$1, $p=$"+num2str(p), 'Interpreter', 'latex', 'FontSize', 13)
xlabel("$\ell$ (in this case, $\ell=m$, $2n = \ell$)", 'Interpreter', 'latex', 'FontSize', 13)
ylabel("Approximation Ratio", 'Interpreter', 'latex', 'FontSize', 13)

ylim([1 5])

lgd = legend([data_points data_points_nesterov theory theory_nesterov], {'$\rho_{\max}$', '$\rho_{\max}$, Nesterov Algorithm', '$C_p\cdot\sqrt{\ell}$', '$\sqrt{\pi/2} \cdot m^{1/p-1/2}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
lgd.FontSize = 15;

%matlab2tikz(['p_equals_',num2str(p),'.tex'])