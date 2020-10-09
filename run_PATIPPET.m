function [mu_list, C_list] = run_PATIPPET(params)

gauss2 = @(x,y, mean, V) exp(-.5 * (([x,y]' - mean)'*(V\([x,y]'-mean))))./ (sqrt((2*pi)^2*abs(det(V))));

t_max = params.tmax + .5;
dt = params.dt;
sigma = params.sigma;
sigma_2 = params.sigma_2;

t_list = 0:dt:t_max;

mu_list = zeros(2, length(t_list));
mu_list(:,1) = params.mu_0;

C_list = zeros(2,2,length(t_list));
C_list(:,:,1) = params.C_0;

event_num = ones(1,params.n_streams);
for i=2:length(t_list)
    t = t_list(i);

    t_past = t_list(i-1);
    C_past = C_list(:,:,i-1);
    mu_past = mu_list(:,i-1);
    
    dmu_sum = 0;
    dC_sum = 0;
    
    for j = 1:params.n_streams
        %dmu_sum = dmu_sum + params.streams{j}.Lambda_star(mu_past, C_past)*(params.streams{j}.mu_star(mu_past, C_past)-mu_past);
        %dC_sum = dC_sum + params.streams{j}.Lambda_star(mu_past, C_past)*(params.streams{j}.C_star(mu_past, C_past)-C_past);
    end
    
    dmu = dt*([mu_past(2); 0] - dmu_sum);
    dC = dt*([sigma^2 + 2*C_past(1,2), C_past(2,2); C_past(2,2), sigma_2^2] - dC_sum);
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

    mu_list(:,i) = mu;
    C_list(:,:,i) = C;
end

if params.display
    event_times = params.streams{1}.event_times;
    e_means = params.streams{1}.e_means;
    expect_func = params.streams{1}.expect_func;
    
    figure()
    subplot(5,1, [1,2,3,4])

    plot_ellipses(mu_list(:,1), C_list(:,:,1), t_max, 5)
    hold on

    for i=2:length(t_list)
        t = t_list(i);
        t_past = t_list(i-1);
        for j = 1:length(event_times)
            if t>=event_times(j) & t_past<event_times(j)

                plot_ellipses(mu_list(:,i), C_list(:,:,i), t_max, 1)
                plot_ellipses(mu_list(:,i+1), C_list(:,:,i+1), t_max, 5)
            end


        end
        if floor(t/.2)>floor(t_past/.2)
            plot_ellipses(mu_list(:,i), C_list(:,:,i), t_max, 1)
        end
    end
    xlim([-.2, t_max])

    ylim([0,2])

    
    ylabel('tempo (beats per sec)')
    
    for i = 1:length(e_means)
        plot([1,1]*e_means(i), [0,2])
    end

    plot([-.2, t_max], [1,1]*params.true_speed)



    subplot(5,1,5)
    plot(t_list, log(expect_func(t_list)))
    xlim([-.2, t_max])
    xlabel('phase')




end


function plot_ellipses(mu, C, tmax, N)
    [X,Y]=meshgrid(-.2:0.01:tmax, 0:0.01:2);
    z = zeros(size(X));
    for n = 1:size(X,1)
        for m = 1:size(X,2)
            z(n,m)=gauss2(X(n,m),Y(n,m), mu, C);
        end
    end
    %surf(X,Y, z)
    contour(X,Y,z,N);
end

end
