function [mu_list, C_list, Pi_list] = run_PATIPPET(params)

gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));

t_max = params.tmax;
dt = params.dt;
sigma = params.sigma;
sigma_2 = params.sigma_2;

ybounds = [params.true_speed-.5, params.true_speed+.5];

t_list = 0:dt:ceil(t_max/dt)*dt;

mu_list = zeros(2, length(t_list));
mu_list(:,1) = params.mu_0;

Pi_list = zeros(2,2,length(t_list));
Pi_list(:,:,1) = params.C_0 + params.mu_0*(params.mu_0)';

C_list = zeros(2,2,length(t_list));
C_list(:,:,1) = params.C_0;

event_num = ones(1,params.n_streams);

for i=2:length(t_list)
    t = t_list(i);
    
    t_past = t_list(i-1);
    mu_past = mu_list(:,i-1);
    C_past = C_list(:,:,i-1);
    Pi_past = Pi_list(:,:,i-1);
    %Pi_past = C_past + mu_past * mu_past'
    %Pi_past == Pi_past_2

    dmu_sum = 0;
    dPi_sum = 0;
    
    for j = 1:params.n_streams
        dmu_sum = dmu_sum + params.streams{j}.Lambda_bar(mu_past, C_past)*(params.streams{j}.mu_bar(mu_past, C_past)-mu_past);
        dPi_sum = dPi_sum + params.streams{j}.Lambda_bar(mu_past, C_past)*(params.streams{j}.Pi_bar(mu_past, C_past)-Pi_past);
    end
%    dPi_sum
%    dmu_sum
    dmu = dt*([mu_past(2); 0] - dmu_sum);
    dPi = dt*([sigma^2 + 2*C_past(1,2) + 2*mu_past(1)*mu_past(2) + mu_past(2)^2*dt, C_past(2,2) + mu_past(2)^2; C_past(2,2) + mu_past(2)^2, sigma_2^2] );
%    dPi = dt*([sigma^2 + 2*Pi_past(1,2) + , Pi_past(2,2) ; Pi_past(2,2) , sigma_2^2] - dPi_sum)
%    dC_1 = dt*([sigma^2 + 2*C_past(1,2) , C_past(2,2) ; C_past(2,2) , sigma_2^2] )
    mu = mu_past+dmu;
    Pi = Pi_past+dPi;
%    'actual dC'
    
    % What's missing right now is a term indicating that Pi_11 is
    % increasing due to phase increase -- it's exactly mu(2)^2, but I can't
    % figure out why that's not included in Pi_12 term.

%    dC = Pi - (mu*mu') - C_past
    C = Pi - mu*mu'
%    A = mu_past*mu_past'
%    Z = mu*mu'
    
    for j = 1:params.n_streams
        if event_num(j) <= length(params.streams{j}.event_times) && (t>=params.streams{j}.event_times(event_num(j)) & t_past<params.streams{j}.event_times(event_num(j)))
            mu_tmp = params.streams{j}.mu_bar(mu_past, C_past);
            Pi = params.streams{j}.Pi_bar(mu_past, C_past);
            mu = mu_tmp;
            event_num(j) = event_num(j)+1;
            'ding'
        end
    end

    mu_list(:,i) = mu;
    Pi_list(:,:,i) = Pi;
    C_list(:,:,i) = Pi - mu*mu';
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
    xlabel('Phase \phi')
    ylabel({'\lambda(\phi_1)'})
    
    subplot(5,1, [1,2,3,4])
    

    h=plot_ellipses(mu_list(:,1), C_list(:,:,1), params.phimax, ybounds, [1/4, 1/4]);
    h.LineWidth = 2;
    h.LineColor = [1,0,0];

    hold on
    plot(mu_list(1,:), mu_list(2,:), 'ko-')
    
    for i=2:length(t_list)
        t = t_list(i);
        t_past = t_list(i-1);
        for j = 1:length(event_times)
            if t>=event_times(j) & t_past<event_times(j)

                h=plot_ellipses(mu_list(:,i-1), C_list(:,:,i-1), params.phimax, ybounds, [1/4, 1/4]);
                h.LineColor = [0,0,1];
                h.LineWidth = 2;
                h=plot_ellipses(mu_list(:,i), C_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
                h.LineWidth = 2;
                h.LineColor = [1,0,0];
            end


        end
        if floor(t/params.dt_ellipse)>floor(t_past/params.dt_ellipse)

            h=plot_ellipses(mu_list(:,i), C_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
            h.LineColor = [0,0,0];
        end
    end
    xlim([-.2, params.phimax])

    ylim(ybounds)

    
    ylabel('Tempo \theta (beats per sec)')
    
    for i = 1:length(e_means)
        plot([1,1]*e_means(i), [0,5], 'b')
    end

    plot([-.2, params.phimax], [1,1]*params.true_speed, 'r')

    sgtitle(params.title)

end


function h = plot_ellipses(mu, C, phimax, ybounds, levels)
    [X,Y]=meshgrid(-.2:0.01:phimax, ybounds(1):0.01:ybounds(2));
    z = zeros(size(X));
    for n = 1:size(X,1)
        for m = 1:size(X,2)
            z(n,m)=gauss2(X(n,m),Y(n,m), mu, C);
        end
    end
    
    top = max(max(z));
    bottom = min(min(z));
    [C,h] = contour(X,Y,z,top*levels + bottom+(1-levels));
end



if params.display_phase
    figure()
    subplot(1,5, [2,3,4,5])
    shadedErrorBar(t_list, mu_list(1,:), 2*sqrt(C_list(1,1,:)))
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
    xlabel('Time (sec)')
    

    subplot(1,5,1)
    for j = 1:params.n_streams
        plot(params.streams{j}.expect_func(t_list), t_list, 'k');
    end
    ylim([0, t_max])
    ylabel('Phase \phi')
    xlabel({'\lambda(\phi)'});
    set(gca,'Yticklabel',[])
    sgtitle(params.title)
    
end


if params.display_tempo
    figure()
    shadedErrorBar(t_list, mu_list(2,:), 2*sqrt(C_list(2,2,:)))
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
    xlabel('Time (sec)')
    ylim([.6,1.4])

    ylabel('Tempo \theta')
    sgtitle(params.title)
    
end
end