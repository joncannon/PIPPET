% © Jonathan Cannon, MIT, 2020
% Simulates PIPPET model with specified parameters.


function [phibar_list, V_list, tap_times] = run_PIPPET(params)

t_max = params.tmax;
dt = params.dt;
sigma_phi = params.sigma_phi;
eta_phi = params.eta_phi;

t_list = 0:dt:t_max;

phibar_list = zeros(size(t_list));
phibar_list(1) = params.phibar_0;

V_list = zeros(size(t_list));
V_list(1) = params.V_0;

event_num = ones(1,params.n_streams);
for i=2:length(t_list)
    t = t_list(i);

    t_past = t_list(i-1);
    V_past = V_list(i-1);
    phibar_past = phibar_list(i-1);
    
    dphibar_sum = 0;
    dV_sum = 0;
    
    for j = 1:params.n_streams
        dphibar_sum = dphibar_sum + params.streams{j}.T_hat(phibar_past, V_past)*(params.streams{j}.phi_hat(phibar_past, V_past)-phibar_past);
    end
    
    dphibar = dt*(1 - dphibar_sum) + sqrt(dt)*eta_phi*randn();
    phibar = phibar_past+dphibar;
    
    for j = 1:params.n_streams
        dV_sum = dV_sum + params.streams{j}.T_hat(phibar_past, V_past)*(params.streams{j}.V_hat(phibar, phibar_past, V_past)-V_past);
    end
    
    dC = dt*(sigma_phi^2 - dV_sum);
    C = V_past+dC;
    
    for j = 1:params.n_streams
        if event_num(j) <= length(params.streams{j}.perceived_event_times) && (t>params.streams{j}.perceived_event_times(event_num(j)) & t_past<=params.streams{j}.perceived_event_times(event_num(j)))
            phibar_tmp = params.streams{j}.phi_hat(phibar, C);
            C = params.streams{j}.V_hat(phibar_tmp, phibar, C);
            phibar = phibar_tmp;
            event_num(j) = event_num(j)+1;
        end
    end

    phibar_list(i) = phibar;
    V_list(i) = C;
end

tap_times = [];
tap_num = 0;
tap_thresh = params.tap_threshold;

for i = 2:length(phibar_list)
    if phibar_list(i) > tap_num+tap_thresh
        tap_times(end+1) = i*dt + (1-tap_thresh);
        tap_num = tap_num+1;
    end
end

if params.display
    figure()
    subplot(1,5, [2,3,4,5])
    shadedErrorBar(t_list, phibar_list, 2*sqrt(V_list))
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
            plot([1,1]*params.streams{j}.perceived_event_times(i), [0,t_max], linespec, 'LineWidth', width);
        end
        if params.tapping
            for i=1:length(tap_times)
                width = .5;
                linespec = 'k:';

                plot([1,1]*tap_times(i), [0,t_max], linespec, 'LineWidth', width);
            end
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