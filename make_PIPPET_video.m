function make_PIPPET_video(filename, params, phibar_list, V_list, x_limits, y_limits, reset_at_end, filetype)

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

s = params.streams{1};

v = VideoWriter(filename,filetype);
open(v);
figure()
%axis tight manual 
set(gca,'nextplot','replacechildren'); 
hold on
hold off
phi_list = x_limits(1):diff(x_limits)/400:x_limits(2);

for i = 1:5:length(phibar_list)-1*reset_at_end

    plot(phi_list,gauss_distribution(phi_list, phibar_list(i), V_list(i)), 'LineWidth', 2);
    hold on
    plot(phi_list,10*s.expect_func(phi_list), ':k', 'LineWidth', 2);
    text(phibar_list(i), y_limits(1)/4, '$\mu_t$','Interpreter','latex', 'FontSize', 14, 'HorizontalAlignment', 'center')
    text(params.streams{1}.e_means(1), y_limits(1)/4, '$\phi_1$','Interpreter','latex', 'FontSize', 14, 'HorizontalAlignment', 'center')
    text(mean(x_limits), y_limits(1)/2, 'Phase $\phi$','Interpreter','latex', 'FontSize', 16, 'HorizontalAlignment', 'center')
    plot([phibar_list(i),phibar_list(i)], [0, y_limits(2)], 'k:')
    hold off
    xlim(x_limits)
    ylim(y_limits)
    writeVideo(v,getframe());
end

if reset_at_end
    plot(phi_list,gauss_distribution(phi_list, phibar_list(end-1), V_list(end-1)), 'LineWidth', 2);
    h = get(gca, 'Children');
    cx = get(h(1), 'Color');
    set(h(1), 'Color', 1-.3*(1-cx));
    hold on
    plot(phi_list,gauss_distribution(phi_list, phibar_list(end), V_list(end)), 'Color', cx, 'LineWidth', 2);
    plot(phi_list,10*s.expect_func(phi_list), ':k', 'LineWidth', 2);
    text(phibar_list(end-1), y_limits(1)/4, '$\mu_t$','Interpreter','latex', 'FontSize', 14, 'HorizontalAlignment', 'center')
    text(phibar_list(end), y_limits(1)/4, '$\mu_{t+}$','Interpreter','latex', 'FontSize', 14, 'HorizontalAlignment', 'center')
    %text(params.streams{1}.e_means(1), y_limits(1)/4, '$\phi_1$','Interpreter','latex', 'FontSize', 14)
    text(mean(x_limits), y_limits(1)/2, 'Phase $\phi$','Interpreter','latex', 'FontSize', 16, 'HorizontalAlignment', 'center')
    plot([phibar_list(end),phibar_list(end)], [0, y_limits(2)], 'k:')
    plot([phibar_list(end-1),phibar_list(end-1)], [0, y_limits(2)], 'k:')
    hold off
    xlim(x_limits)
    ylim(y_limits)
    legend('Pre-event', 'Post-event')
    writeVideo(v,getframe());
end
close(v)