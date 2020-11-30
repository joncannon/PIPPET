function [xbar_list, Sigma_list] = run_PATIPPET(params)

gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));

t_max = params.tmax;
dt = params.dt;
sigma_phi = params.sigma_phi;
sigma_theta = params.sigma_theta;

ybounds = [params.true_speed-.5, params.true_speed+.5];

t_list = 0:dt:ceil(t_max/dt)*dt;

xbar_list = zeros(2, length(t_list));
xbar_list(:,1) = params.xbar_0;

Sigma_list = zeros(2,2,length(t_list));
Sigma_list(:,:,1) = params.Sigma_0;

event_num = ones(1,params.n_streams);

for i=2:length(t_list)
    t = t_list(i);

    t_past = t_list(i-1);
    Sigma_past = Sigma_list(:,:,i-1);
    xbar_past = xbar_list(:,i-1);
    
    dxbar_sum = 0;
    dSigma_sum = 0;
    
    for j = 1:params.n_streams
        dxbar_sum = dxbar_sum + params.streams{j}.Lambda_hat(xbar_past, Sigma_past)*(params.streams{j}.x_hat(xbar_past, Sigma_past)-xbar_past);
    end
    
    dxbar = dt*([xbar_past(2); 0] - dxbar_sum);
    xbar = xbar_past+dxbar;
    
    for j = 1:params.n_streams
        dSigma_sum = dSigma_sum + params.streams{j}.Lambda_hat(xbar_past, Sigma_past)*(params.streams{j}.Sigma_hat(xbar, xbar_past, Sigma_past)-Sigma_past);
    end
    dSigma = dt*([sigma_phi^2 + 2*Sigma_past(1,2), Sigma_past(2,2); Sigma_past(2,2), sigma_theta^2] - dSigma_sum);
    Sigma = Sigma_past+dSigma;
    
    for j = 1:params.n_streams
        if event_num(j) <= length(params.streams{j}.event_times) && (t>=params.streams{j}.event_times(event_num(j)) & t_past<params.streams{j}.event_times(event_num(j)))
            xbar_tmp = params.streams{j}.x_hat(xbar, Sigma);
            Sigma_tmp = params.streams{j}.Sigma_hat(xbar_tmp, xbar, Sigma);
            Sigma = Sigma_tmp;
            xbar = xbar_tmp;
            event_num(j) = event_num(j)+1;
        end
    end

    xbar_list(:,i) = xbar;
    Sigma_list(:,:,i) = Sigma;
end

if params.display_phasetempo
    phi_list = -.2:dt:params.phimax;
    
    event_times = params.streams{1}.event_times;
    e_means = params.streams{1}.e_means;
    expect_func = params.streams{1}.expect_func;
    
    figure()
    

    subplot(5,1,5)
    plot(phi_list, expect_func(phi_list), 'k')
    xlim([-.2, params.phimax])
    xlabel('Phase $\phi$','Interpreter','Latex')
    ylabel({'Expectation';'$\tau(\phi)$'},'Interpreter','Latex')
    
    subplot(5,1, [1,2,3,4])
    

    h=plot_ellipses(xbar_list(:,1), Sigma_list(:,:,1), params.phimax, ybounds, [1/4, 1/4]);
    h.LineWidth = 2;
    h.LineColor = [1,0,0];

    hold on
    plot(xbar_list(1,:), xbar_list(2,:), 'ko-')
    
    for i=2:length(t_list)
        t = t_list(i);
        t_past = t_list(i-1);
        for j = 1:length(event_times)
            if t>=event_times(j) & t_past<event_times(j)

                h=plot_ellipses(xbar_list(:,i-1), Sigma_list(:,:,i-1), params.phimax, ybounds, [1/4, 1/4]);
                h.LineColor = [0,0,1];
                h.LineWidth = 2;
                h=plot_ellipses(xbar_list(:,i), Sigma_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
                h.LineWidth = 2;
                h.LineColor = [1,0,0];
            end


        end
        if floor(t/params.dt_ellipse)>floor(t_past/params.dt_ellipse)

            h=plot_ellipses(xbar_list(:,i), Sigma_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
            h.LineColor = [0,0,0];
        end
    end
    xlim([-.2, params.phimax])

    ylim(ybounds)

    
    ylabel('Tempo $\theta$ (beats per sec)','Interpreter','Latex')
    
    for i = 1:length(e_means)
        plot([1,1]*e_means(i), [0,5], 'b')
    end

    plot([-.2, params.phimax], [1,1]*params.true_speed, 'r')

    sgtitle(params.title)

end


function h = plot_ellipses(xbar, Sigma, phimax, ybounds, levels)
    [X,Y]=meshgrid(-.2:0.01:phimax, ybounds(1):0.01:ybounds(2));
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



if params.display_phase
    figure()
    subplot(1,5, [2,3,4,5])
    shadedErrorBar(t_list, xbar_list(1,:), 2*sqrt(Sigma_list(1,1,:)))
    ylim([0, t_max])
    hold on
    for j = 1:params.n_streams
        for i=1:length(params.streams{j}.event_times)
            width = .5;
            linespec = 'r';
            if params.streams{j}.highlight_event_indices(i)==0
                linespec = 'r-.';
            elseif params.streams{j}.highlight_event_indices(i)==2
                width = 1.5;
            end
            plot([1,1]*params.streams{j}.event_times(i), [0,t_max], linespec, 'LineWidth', width);
        end

        for i=1:length(params.streams{j}.e_means)
            width = .5;
            linespec = 'b';
            if params.streams{j}.highlight_expectations(i)==0
                linespec = 'b-.';
            elseif params.streams{j}.highlight_expectations(i)==2
                width = 1.5;
            end
            plot([0,t_max], [1,1]*params.streams{j}.e_means(i), linespec, 'LineWidth', width)
        end
    end
    xlabel('Time (sec)','Interpreter','Latex')
    

    subplot(1,5,1)
    for j = 1:params.n_streams
        plot(params.streams{j}.expect_func(t_list), t_list, 'k');
    end
    ylim([0, t_max])
    ylabel('Phase $\phi$','Interpreter','Latex')
    xlabel({'Expectation';'$\tau(\phi)$'},'Interpreter','Latex');
    set(gca,'Yticklabel',[])
    sgtitle(params.title)
    
end


if params.display_tempo
    figure()
    shadedErrorBar(t_list, xbar_list(2,:), 2*sqrt(Sigma_list(2,2,:)))
    hold on
    for j = 1:params.n_streams
        for i=1:length(params.streams{j}.event_times)
            width = .5;
            linespec = 'r';
            if params.streams{j}.highlight_event_indices(i)==0
                linespec = 'r-.';
            elseif params.streams{j}.highlight_event_indices(i)==2
                width = 1.5;
            end
            plot([1,1]*params.streams{j}.event_times(i), [0,t_max], linespec, 'LineWidth', width);
        end

    end
    plot([t_list(1), t_list(end)], [1,1]*params.true_speed, 'r')
    xlabel('Time (sec)','Interpreter','Latex')
    ylim([.6,1.4])

    ylabel('Tempo $\theta$','Interpreter','Latex')
    sgtitle(params.title)
    
end
end