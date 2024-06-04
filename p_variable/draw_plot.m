p = 1.5;
%p = 2;


n_range = 1:6;
m_range = 2:2:12;

clf;
hold on

for m = m_range
    worst_case = 0;
    for n = n_range
        try
            filename = ['data/',num2str(p),'_',num2str(n),'_',num2str(m),'_',num2str(m),'_linear.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points_linear = plot(m,worst_case,'ok');
end

for m = m_range
    worst_case = 0;
    for n = n_range
        try
            filename = ['data/',num2str(p),'_',num2str(n),'_',num2str(m),'_',num2str(m),'_semidefinite.mat'];
            load(filename)
            worst_case = max(worst_case, 1/fval);
        catch
            continue
        end
    end
    data_points_semidefinite = plot(m,worst_case,'xk');
end

gamma_r = @(r) (2^(r/2)/sqrt(pi) * gamma((r+1)/2))^(1/r);

x = linspace(2, 12, 1000);
theory_linear = plot(x, gamma_r(p)/gamma_r(1) * sqrt(x), 'k');
if p>= 2
    theory_semidefinite = plot(x, gamma_r(p)/gamma_r(1) + zeros(size(x)), 'k--');
else
    theory_semidefinite = plot(x, gamma_r(2)*x.^(1/p-1/2)./gamma_r(1) + zeros(size(x)), 'k--');
end


title("Worst-case approximation ratio for $q=1$, $p="+num2str(p)+"$", 'Interpreter', 'latex')
xlabel("$m$ (in this case, $m=l$, $m = 2n$)", 'Interpreter', 'latex')
ylabel("Approximation Ratio", 'Interpreter', 'latex')

ylim([1 5])

if p>= 2
    lgd = legend([data_points_linear data_points_semidefinite theory_linear theory_semidefinite], {'Linear Relaxation', 'Semidefinite Relaxation', '$\sqrt{m}\cdot\gamma_p/\gamma_1$', '$\gamma_p/\gamma_1$'}, 'Interpreter', 'latex', 'Location', 'northwest');
else
    lgd = legend([data_points_linear data_points_semidefinite theory_linear theory_semidefinite], {'Linear Relaxation', 'Semidefinite Relaxation', '$\sqrt{m}\cdot\gamma_p/\gamma_1$', '$\sqrt{\pi/2} \cdot l^{1/p-1/2}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
end


%matlab2tikz(['p_equals_',num2str(p),'.tex'])