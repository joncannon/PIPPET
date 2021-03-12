function make_PATIPPET_video(filename, params, phibar_list, V_list, x_limits, y_limits, filetype)

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));


event_times = params.streams{1}.event_times;
    e_means = params.streams{1}.e_means;
    expect_func = params.streams{1}.expect_func;

v = VideoWriter(filename,filetype);
open(v);
fig=figure()
%axis tight manual 
tiledlayout(5,1)

nexttile(5)
phi_list = x_limits(1):diff(x_limits)/400:x_limits(2);
    plot(phi_list, params.streams{1}.expect_func(phi_list), 'k')
    xlim([-.2, params.phimax])
    xlabel('Phase $\phi$','Interpreter','Latex')
    ylabel({'Expectation';'$\lambda(\phi)$'},'Interpreter','Latex')

    nexttile(1, [4,1])
    
set(gca,'nextplot','replacechildren'); 


i=1

    plot(phibar_list(1,i), phibar_list(2,i), 'ko')
    
hold on
            h=plot_ellipses(phibar_list(:,i), V_list(:,:,i), x_limits(2), y_limits, [1/4, 1/4]);
            h.LineColor = [0,.4,.8];
            %h.LineStyle = ':';

    xlim([-.2, params.phimax])

    ylim(y_limits)

    
    ylabel('Tempo $\theta$ (beats per sec)','Interpreter','Latex')
    xlabel('Phase $\phi$','Interpreter','Latex')
    for i = 1:length(e_means)
        plot([1,1]*e_means(i), [0,5], 'k:')
    end

    plot([-.2, params.phimax], [1,1]*params.true_speed, 'r')

    F = getframe(fig);
    hold off


for i = 1:33:size(phibar_list,2)
    


    plot(phibar_list(1,i), phibar_list(2,i), 'ko')
    
hold on
            h=plot_ellipses(phibar_list(:,i), V_list(:,:,i), x_limits(2), y_limits, [1/4, 1/4]);
            h.LineColor = [0,.4,.8];
            %h.LineStyle = ':';

    xlim([-.2, params.phimax])

    ylim(y_limits)

    
    ylabel('Tempo $\theta$ (beats per sec)','Interpreter','Latex')
    xlabel('Phase $\phi$','Interpreter','Latex')
    for i = 1:length(e_means)
        plot([1,1]*e_means(i), [0,5], 'k:')
    end

    plot([-.2, params.phimax], [1,1]*params.true_speed, 'r')

    F = getframe(fig);
    writeVideo(v,F);
    hold off
end

close(v)
end



function h = plot_ellipses(xbar, Sigma, phimax, ybounds, levels)
    gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));
    [X,Y]=meshgrid(-.2:0.05:phimax, ybounds(1):0.01:ybounds(2));
    z = zeros(size(X));
    for n = 1:size(X,1)
        for m = 1:size(X,2)
            z(n,m)=gauss2(X(n,m),Y(n,m), xbar, Sigma);
        end
    end
    
    top = max(max(z));
    bottom = min(min(z));
    [C,h] = contour(X,Y,z,top*levels + bottom+(1-levels));
end