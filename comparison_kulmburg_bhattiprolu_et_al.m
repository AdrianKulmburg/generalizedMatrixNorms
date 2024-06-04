g_r = @(r) ((2.^(r/2)./sqrt(pi)) .* gamma((r+1)/2)).^(1./r);

holder_r = @(r) r./(r-1);


kulmburg = @(p,q) g_r(p) .* g_r(holder_r(q));
bhattiprolu_et_al = @(p,q) 1.00863./(g_r(holder_r(p)) .* g_r(q) * log(1+sqrt(2)));

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12);
ratio_pre = @(q) kulmburg(2, q) ./ bhattiprolu_et_al(2, q) - 1;
q_star = fsolve(ratio_pre, 1.35, options);


q_range = linspace(q_star, 2, 100);

p_q_values = [];

for q = q_range
    ratio = @(p) kulmburg(p, q) ./ bhattiprolu_et_al(p, q) - 1;
    p_star = fsolve(ratio, 3, options);


    p_q_values = [p_q_values [p_star; q]];

end

RTPH_Blue = [0 92 171]./255;
zone_upperLeft = polyshape([2 2 p_q_values(1,:)],[2 q_star p_q_values(2,:)]);
plot(zone_upperLeft, 'EdgeColor','none', 'FaceColor', RTPH_Blue, 'FaceAlpha', 1)
ylim([1 2])

axis square

title("", 'Interpreter', 'latex')
xlabel("$p$", 'Interpreter', 'latex')
ylabel("$q$", 'Interpreter', 'latex')





kulmburg2 = @(p) g_r(p) / g_r(1);

ratio_q_equals_1 = @(p) kulmburg2(p) ./ bhattiprolu_et_al(p,1) - 1;

p_star = fsolve(ratio_q_equals_1, 3.8, options) % Bottom bar will be added manually letter

matlab2tikz(['comparison_kulmburg_bhattiprolu_et_al.tex'])