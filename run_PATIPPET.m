% © Jonathan Cannon, MIT, 2020
% Simulates PATIPPET model with specified parameters.

function [mu_list, V_list, surprisal_prepost, grad_surprisal] = run_PATIPPET(params)

gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));

t_max = params.tmax;
dt = params.dt;
sigma_phi = params.sigma_phi;
sigma_theta = params.sigma_theta;

ybounds = [params.true_speed-.6, params.true_speed+.6];

t_list = 0:dt:ceil(t_max/dt)*dt;

mu_list = zeros(2, length(t_list));
mu_list(:,1) = params.mu_0;

V_list = zeros(2,2,length(t_list));
V_list(:,:,1) = params.V_0;

surprisal_prepost = zeros([numel(t_list), params.n_streams, 2]);
grad_surprisal = zeros([numel(t_list),1]);

event_num = ones(1,params.n_streams);

for i=2:length(t_list)
    t = t_list(i);

    t_past = t_list(i-1);
    V_past = V_list(:,:,i-1);
    mu_past = mu_list(:,i-1);
    
    dmu_sum = 0;
    dV_sum = 0;
    
    grad_surprisal_sum = 0;
    
    for j = 1:params.n_streams
        dmu_sum = dmu_sum + params.streams{j}.Lambda(mu_past, V_past)*(params.streams{j}.mu_hat(mu_past, V_past)-mu_past);
    end
    
    dmu = dt*([mu_past(2); 0] - dmu_sum);
    mu = mu_past+dmu;
    
    for j = 1:params.n_streams
        dV_sum = dV_sum + params.streams{j}.Lambda(mu_past, V_past)*(params.streams{j}.V_hat(mu, mu_past, V_past)-V_past);
    end
    
    if params.TDDM
        extra_noise = mu(2);
    else
        extra_noise = 1;
    end
    
    dV = dt*([sigma_phi^2*extra_noise + 2*V_past(1,2), V_past(2,2); V_past(2,2), sigma_theta^2] - dV_sum);
    V = V_past+dV;
    
    for j = 1:params.n_streams
        if event_num(j) <= length(params.streams{j}.event_times) && (t>=params.streams{j}.event_times(event_num(j)) & t_past<params.streams{j}.event_times(event_num(j)))
            mu_tmp = params.streams{j}.mu_hat(mu, V);
            V_tmp = params.streams{j}.V_hat(mu_tmp, mu, V);
            V = V_tmp;
            mu = mu_tmp;
            event_num(j) = event_num(j)+1;
            surprisal_prepost(i,j,1) = -log(params.streams{j}.Lambda(mu_past, V_past)*dt);
            surprisal_prepost(i,j,2) = -log(params.streams{j}.Lambda(mu, V)*dt);
            grad_surprisal_sum = grad_surprisal_sum ...
                +(-log(params.streams{j}.Lambda(mu_past+.01, V_past)*dt) + log(params.streams{j}.Lambda(mu_past-.01, V_past)*dt))/.02;
        else
            surprisal_prepost(i,j,1) = -log(1-params.streams{j}.Lambda(mu_past, V_past)*dt);
            surprisal_prepost(i,j,2) = -log(1-params.streams{j}.Lambda(mu, V)*dt);
            grad_surprisal_sum = grad_surprisal_sum ...
                +(-log(1-params.streams{j}.Lambda(mu_past+.01, V_past)*dt) + log(1-params.streams{j}.Lambda(mu_past-.01, V_past)*dt))/.02;
        end
    end

    mu_list(:,i) = mu;
    V_list(:,:,i) = V;
    grad_surprisal(i) = grad_surprisal_sum;
end

if params.display_phasetempo
    phi_list = -.2:dt:params.phimax;
    
    event_times = params.streams{1}.event_times;
    e_means = params.streams{1}.e_means;
    expect_func = params.streams{1}.expect_func;
    
    figure()

    subplot(5,1,5)
    plot(phi_list, expect_func(phi_list), 'b')
    xlim([-.2, params.phimax])
    xlabel('Phase $\phi$','Interpreter','Latex')
    ylabel({'Expectation';'$\lambda(\phi)$'},'Interpreter','Latex')
    
    subplot(5,1, [1,2,3,4])
    

    h=plot_ellipses(mu_list(:,1), V_list(:,:,1), params.phimax, ybounds, [1/4, 1/4]);
    h.LineWidth = 2;
    h.LineColor = [1,0,0];

    hold on
    plot(mu_list(1,:), mu_list(2,:), 'ko-')
    
    for i=2:length(t_list)
        t = t_list(i);
        t_past = t_list(i-1);
        for j = 1:length(event_times)
            if t>=event_times(j) & t_past<event_times(j)

                h=plot_ellipses(mu_list(:,i-1), V_list(:,:,i-1), params.phimax, ybounds, [1/4, 1/4]);
                h.LineColor = [0,.7,0];
                h.LineWidth = 2;
                h=plot_ellipses(mu_list(:,i), V_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
                h.LineWidth = 2;
                h.LineColor = [1,0,0];
            end


        end
        if floor(t/params.dt_ellipse)>floor(t_past/params.dt_ellipse)

            h=plot_ellipses(mu_list(:,i), V_list(:,:,i), params.phimax, ybounds, [1/4, 1/4]);
            h.LineColor = [0,0,0];
            h.LineStyle = ':';
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





if params.display_phase
    figure()
    subplot(1,5, [2,3,4,5])
    shadedErrorBar(t_list, mu_list(1,:), 2*sqrt(V_list(1,1,:)))
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
        plot(params.streams{j}.expect_func(t_list), t_list, 'b');
    end
    ylim([0, t_max])
    ylabel('Phase $\phi$','Interpreter','Latex')
    xlabel({'Expectation';'$\lambda(\phi)$'},'Interpreter','Latex');
    set(gca,'Yticklabel',[])
    sgtitle(params.title)
    
end


if params.display_tempo
    figure()
    shadedErrorBar(t_list, mu_list(2,:), 2*sqrt(V_list(2,2,:)))
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


function h = plot_ellipses(mu, V, phimax, ybounds, levels)
    gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));
    [X,Y]=meshgrid(-.2:0.01:phimax, ybounds(1):0.01:ybounds(2));
    z = zeros(size(X));
    for n = 1:size(X,1)
        for m = 1:size(X,2)
            z(n,m)=gauss2(X(n,m),Y(n,m), mu, V);
        end
    end
    
    top = max(max(z));
    bottom = min(min(z));
    [C,h] = contour(X,Y,z,top*levels + bottom+(1-levels));
end