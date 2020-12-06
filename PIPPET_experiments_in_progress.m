% (C) Jonathan Cannon, MIT, 2020


%Tapping: autocorr
%We introduce a tapping finger that aims to tap with the metronome but introduces motor noise. This tap is used as a second sensory stream for phase estimation. We calculate autocorrelation to determine the phase correction response as it corrects motor noise, and we compare it to the phase correction response to a phase shift.
n_events = 100;

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

e_means = 1:n_events;
e_vars = .0001*ones(1,n_events);
e_lambdas = .02*ones(1,n_events);
lambda_0 = .00001;

tap_vars = .0005*ones(1,n_events);
tap_lambdas = .02*ones(1,n_events);
tap_lambda_0 = .00001;
tap_noise = .02;

mu_i_list = @(mu, C) (mu/C + e_means./e_vars)./(1/C + 1./e_vars);
C_i_list = @(C) 1./(1/C + 1./e_vars);
Lambda_i_list = @(mu, C) e_lambdas .* gauss_distribution(mu, e_means, e_vars+C);
Lambda_star = @(mu, C) lambda_0 + sum(Lambda_i_list(mu, C));
mu_star = @(mu, C) (lambda_0*mu + sum(Lambda_i_list(mu, C) .* mu_i_list(mu,C)))/Lambda_star(mu,C);
C_star= @(mu,C) (lambda_0*C + sum(Lambda_i_list(mu, C) .* (C_i_list(C) + (mu_i_list(mu, C)-mu).^2)))/Lambda_star(mu,C);


tap_mu_i_list = @(mu, C) (mu/C + e_means./tap_vars)./(1/C + 1./tap_vars);
tap_C_i_list = @(C) 1./(1/C + 1./tap_vars);
tap_Lambda_i_list = @(mu, C) tap_lambdas .* gauss_distribution(mu, e_means, tap_vars+C);
tap_Lambda_star = @(mu, C) tap_lambda_0 + sum(tap_Lambda_i_list(mu, C));
tap_mu_star = @(mu, C) (tap_lambda_0*mu + sum(tap_Lambda_i_list(mu, C) .* tap_mu_i_list(mu,C)))/tap_Lambda_star(mu,C);
tap_C_star= @(mu,C) (tap_lambda_0*C + sum(tap_Lambda_i_list(mu, C) .* (tap_C_i_list(C) + (tap_mu_i_list(mu, C)-mu).^2)))/tap_Lambda_star(mu,C);


t_max = n_events+5;
dt = 0.01;
t_list = 0:dt:t_max;


template = zeros(length(t_list));
for i = 1:length(e_means)
    template = template + e_lambdas(i)*gauss_distribution(t_list, e_means(i), e_vars(i))';
end
template = template + lambda_0;


mu_list = zeros(size(t_list));
mu_list(1) = 0;

C_list = zeros(size(t_list));
C_list(1) = .0002;

sigma = 0.1;

event_times = e_means;
tap_times = [];

event_num = 1;
next_tap = -1;

for i=2:length(t_list)
    t = t_list(i);
    
    t_past = t_list(i-1);
    C_past = C_list(i-1);
    mu_past = mu_list(i-1);
    dmu = dt*(1 - Lambda_star(mu_past, C_past)*(mu_star(mu_past, C_past)-mu_past));
    dC = dt*(sigma^2 - Lambda_star(mu_past, C_past)*(C_star(mu_past, C_past)-C_past));
    mu = mu_past+dmu;
    C = C_past+dC;
    
    
    if event_num<=length(event_times) && (t>=event_times(event_num) & t_past<event_times(event_num))
        
        mu = mu_star(mu, C);
        C = C_star(mu, C);
        event_num = event_num+1;
    end
    
    if t>=next_tap & t_past<next_tap
        
        mu = tap_mu_star(mu, C);
        C = tap_C_star(mu, C);
        tap_times(end+1) = t;
    end
    
    if floor(mu_list(i-1)-.5) < floor(mu-.5)
        next_tap = t+.5+ tap_noise*randn();
    end
    
    mu_list(i) = mu;
    C_list(i) = C;
end

async = tap_times(1:length(event_times)) - event_times;
AC = autocorr(async)
alpha = 1-AC(2)


%Tapping: phase shift
%The phase correction response to a phase shift, with same parameters as above (minus motor noise)
n_events = 10;
phase_shift = 0.05;

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

e_means = 1:n_events;
e_vars = .0001*ones(1,n_events);
e_lambdas = .02*ones(1,n_events);
lambda_0 = .00001;

tap_vars = .0005*ones(1,n_events);
tap_lambdas = .02*ones(1,n_events);
tap_lambda_0 = 0;
tap_noise = 0;

mu_i_list = @(mu, C) (mu/C + e_means./e_vars)./(1/C + 1./e_vars);
C_i_list = @(C) 1./(1/C + 1./e_vars);
Lambda_i_list = @(mu, C) e_lambdas .* gauss_distribution(mu, e_means, e_vars+C);
Lambda_star = @(mu, C) lambda_0 + sum(Lambda_i_list(mu, C));
mu_star = @(mu, C) (lambda_0*mu + sum(Lambda_i_list(mu, C) .* mu_i_list(mu,C)))/Lambda_star(mu,C);
C_star= @(mu,C) (lambda_0*C + sum(Lambda_i_list(mu, C) .* (C_i_list(C) + (mu_i_list(mu, C)-mu).^2)))/Lambda_star(mu,C);


tap_mu_i_list = @(mu, C) (mu/C + e_means./tap_vars)./(1/C + 1./tap_vars);
tap_C_i_list = @(C) 1./(1/C + 1./tap_vars);
tap_Lambda_i_list = @(mu, C) tap_lambdas .* gauss_distribution(mu, e_means, tap_vars+C);
tap_Lambda_star = @(mu, C) tap_lambda_0 + sum(tap_Lambda_i_list(mu, C));
tap_mu_star = @(mu, C) (tap_lambda_0*mu + sum(tap_Lambda_i_list(mu, C) .* tap_mu_i_list(mu,C)))/tap_Lambda_star(mu,C);
tap_C_star= @(mu,C) (tap_lambda_0*C + sum(tap_Lambda_i_list(mu, C) .* (tap_C_i_list(C) + (tap_mu_i_list(mu, C)-mu).^2)))/tap_Lambda_star(mu,C);


t_max = n_events + 2;
dt = 0.001;
t_list = 0:dt:t_max;


template = zeros(length(t_list));
for i = 1:length(e_means)
    template = template + e_lambdas(i)*gauss_distribution(t_list, e_means(i), e_vars(i))';
end
template = template + lambda_0;


mu_list = zeros(size(t_list));
mu_list(1) = 0;

C_list = zeros(size(t_list));
C_list(1) = .0002;

sigma = 0.1;

event_times = e_means;
event_times(5:end) = event_times(5:end) - phase_shift;
tap_times = [];

event_num = 1;
next_tap = -1;

for i=2:length(t_list)
    t = t_list(i);
    
    t_past = t_list(i-1);
    C_past = C_list(i-1);
    mu_past = mu_list(i-1);
    dmu = dt*(1 - Lambda_star(mu_past, C_past)*(mu_star(mu_past, C_past)-mu_past));
    dC = dt*(sigma^2 - Lambda_star(mu_past, C_past)*(C_star(mu_past, C_past)-C_past));
    mu = mu_past+dmu;
    C = C_past+dC;
    
    
    if event_num<=length(event_times) && (t>=event_times(event_num) & t_past<event_times(event_num))
        
        mu = mu_star(mu, C);
        C = C_star(mu, C);
        event_num = event_num+1;
    end
    
    if t>=next_tap & t_past<next_tap
        
        mu = tap_mu_star(mu, C);
        C = tap_C_star(mu, C);
        tap_times(end+1) = t;
    end
    
    if floor(mu_list(i-1)-.5) < floor(mu-.5)
        next_tap = t+.5+ tap_noise*randn();
    end
    
    mu_list(i) = mu;
    C_list(i) = C;
end

async = tap_times(1:length(event_times)) - event_times
alpha = 1-(async(6)/async(5));



