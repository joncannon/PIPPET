function [mu_list, C_list] = run_PIPPET(params)

t_max = params.tmax;
dt = params.dt;
sigma = params.sigma;

t_list = 0:dt:t_max;

mu_list = zeros(size(t_list));
mu_list(1) = params.mu_0;

C_list = zeros(size(t_list));
C_list(1) = params.C_0;

event_num = ones(1,params.n_streams);
for i=2:length(t_list)
    t = t_list(i);

    t_past = t_list(i-1);
    C_past = C_list(i-1);
    mu_past = mu_list(i-1);
    
    dmu_sum = 0;
    dC_sum = 0;
    
    for j = 1:params.n_streams
        dmu_sum = dmu_sum + params.streams{j}.Lambda_star(mu_past, C_past)*(params.streams{j}.mu_star(mu_past, C_past)-mu_past);
        dC_sum = dC_sum + params.streams{j}.Lambda_star(mu_past, C_past)*(params.streams{j}.C_star(mu_past, C_past)-C_past);
    end
    
    dmu = dt*(1 - dmu_sum);
    dC = dt*(sigma^2 - dC_sum);
    mu = mu_past+dmu;
    C = C_past+dC;
    
    for j = 1:params.n_streams
        if event_num(j) <= length(params.streams{j}.event_times) && (t>params.streams{j}.event_times(event_num(j)) & t_past<=params.streams{j}.event_times(event_num(j)))
            mu_tmp = params.streams{j}.mu_star(mu, C);
            C = params.streams{j}.C_star(mu, C);
            mu = mu_tmp;
            event_num(j) = event_num(j)+1;
        end
    end

    mu_list(i) = mu;
    C_list(i) = C;
end

if params.display
    figure()
    subplot(1,5, [2,3,4,5])
    shadedErrorBar(t_list, mu_list, 2*sqrt(C_list))
    ylim([0, t_max])
    hold on
    for j = 1:params.n_streams
        for i=1:length(params.streams{j}.event_times)
            plot([1,1]*params.streams{j}.event_times(i), [0,t_max], 'r')
        end

        for i=1:length(params.streams{j}.e_means)
            plot([0,t_max], [1,1]*params.streams{j}.e_means(i), 'b')
        end
    end
    xlabel('Time (s)')

    subplot(1,5,1)
    for j = 1:params.n_streams
        plot(log(params.streams{j}.expect_func(t_list)), t_list, 'k');
    end
    ylim([0, t_max])
    ylabel('Phase \phi')
    xlabel({'log likelihood';'log(\lambda(\phi))'});
    sgtitle(params.title)
end