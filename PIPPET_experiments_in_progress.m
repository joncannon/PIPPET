% (C) Jonathan Cannon, MIT, 2020

%Tapping autocorr

n_events = 1000;
tap_thresh = .7;
event_times = 1:n_events;

expectation_variance_list = [.001, .01];
condition_list = {'NT', 'ASD'};

async_list = {[],[]};
AC_list = {[],[]};
alpha_list = [0,0];

for condition = 1:length(expectation_variance_list)
    params = PIPPET_params(...
        'means_unit', [1],...
        'variance_unit', [expectation_variance_list(condition)],...
        'tau_unit', [.02],...
        'tau_0', .0001,...
        'expected_cycles', n_events,...
        'expected_period', 1,...
        'event_times', event_times,...
        'display', false,...
        'eta_phi', .05,...
        'title', 'Tap along');
    [phibar_list, V_list] = run_PIPPET(params);

    tap_times = [];
    tap_num = 0;
    %subplot(1,5, [2,3,4,5])
    %hold on
    for i = 2:length(phibar_list)
        if phibar_list(i) > tap_num+tap_thresh
            tap_times(end+1) = i*params.dt + (1-tap_thresh);
            tap_num = tap_num+1;
            %plot(tap_times(end)*[1,1], [0,5])
        end
    end
    async =  tap_times(1:length(event_times)) - event_times;
    async_list{condition} = async
    AC_list{condition} = autocorr(async);
    
    figure()
    plot(async(1:end-1), async(2:end), 'o')
    title(strcat('Consecutive asynchronies: ', condition_list{condition}));
    xlim([-.25, .25])
    ylim([-.25, .25])
end



figure()
hold on
for condition = 1:length(expectation_variance_list)
    
    plot(AC_list{condition})
end
title('Asynchrony autocorrelation')
