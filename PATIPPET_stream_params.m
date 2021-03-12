% © Jonathan Cannon, MIT, 2020
% Creates parameters (including auxilliary functions) for a single event stream in the PATIPPET model.
% Inputs:
%   means_unit:         one repetition of a pattern of expected event times phi_i
%   variance_unit:      one repetition of a pattern of expected event timing variances v_i
%   lambda_unit:           one repetition of a pattern of event expectation strengths lambda_i
%   lambda_0
%   expected_cycles:    number of cycles that the patterns repeat to create
%                           temporal expectation template
%   expected_period:    period with which the patterns repeat
%   event_times:        time points at which events actually occur in this
%                           stream
%   highlight_expectations:     Display weights for lines marking expected timepoints
%   highlight_event_indices:    Display weights for lines marking event times
%   corrected:          whether PATIPPET model is corrected to normalize
%                           event rate over phase (rather than over time) 
% Output:
%   p:   parameter set for stream

function p = PATIPPET_stream_params(means_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, corrected)

p = struct();

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

p.e_means = [];
p.event_times = event_times;
p.lambda_0 = lambda_0;
p.highlight_event_indices = highlight_event_indices;
p.highlight_expectations = [];

for i=1:expected_cycles
    p.e_means = [p.e_means, means_unit + (i-1)*expected_period];
    p.highlight_expectations = [p.highlight_expectations, highlight_expectations];
end

p.e_vars = repmat(variance_unit, [1,expected_cycles]);
p.e_lambdas = repmat(lambda_unit, [1,expected_cycles]);



mu_i = @(mu, V, i) (inv(V) + [1./p.e_vars(i), 0;0,0])\(V\mu + [p.e_means(i)./p.e_vars(i);0]);
K_i = @(V, i) inv(inv(V) + [1./p.e_vars(i), 0;0,0]);
Lambda_i = @(mu, V, i) p.e_lambdas(i) .* gauss_distribution(mu(1), p.e_means(i), p.e_vars(i)+V(1,1));


function out = Lambda_uncorrected(mu, V)
    out = p.lambda_0;
    for i = 1:length(p.e_means)
         out = out + Lambda_i(mu, V, i);
    end
end

function out = V_hat_uncorrected(mu_new, mu_old, V)

    V_sum = p.lambda_0 * (V + (mu_old-mu_new)*(mu_old-mu_new)');
    for i = 1:length(p.e_means)
        mu_i_tmp = mu_i(mu_old, V, i);
        K_i_tmp = K_i(V, i);
        V_sum = V_sum + Lambda_i(mu_old, V, i) * (K_i_tmp + (mu_i_tmp-mu_new)*(mu_i_tmp-mu_new)');
    end
    out = V_sum/Lambda(mu_old, V);
end
   
function out = mu_hat_uncorrected(mu, V)
    mu_sum = p.lambda_0*mu;
    for i = 1:length(p.e_means)
         mu_sum = mu_sum + Lambda_i(mu, V, i) * mu_i(mu, V, i);
    end
    out = mu_sum/Lambda(mu, V);
end


function out = Lambda(mu, V)
    out = p.lambda_0 * mu(2);
    for i = 1:length(p.e_means)
        mu_i_tmp = mu_i(mu, V, i);
        out = out + Lambda_i(mu, V, i)*mu_i_tmp(2);
    end
end

function out = V_hat(mu_new, mu_old, V)

    V_sum = p.lambda_0 * (mu_old(2)*(V + (mu_old-mu_new)*(mu_old-mu_new)') + (mu_old - mu_new)*V(2,:) +V(:,2)*(mu_old - mu_new)');
    for i = 1:length(p.e_means)
        mu_i_tmp = mu_i(mu_old, V, i);
        K_i_tmp = K_i(V, i);
        V_sum = V_sum + Lambda_i(mu_old, V, i) * (mu_i_tmp(2)*(K_i_tmp + (mu_i_tmp-mu_new)*(mu_i_tmp-mu_new)') + (mu_i_tmp - mu_new)*K_i_tmp(2,:) + K_i_tmp(:,2)*(mu_i_tmp - mu_new)');
    end
    out = V_sum/Lambda(mu_old, V);
end
   
function out = mu_hat(mu, V)
    mu_sum = p.lambda_0*(V(:,2)+mu*mu(2));
    for i = 1:length(p.e_means)
        K_i_tmp = K_i(V, i);
        mu_i_tmp = mu_i(mu, V, i);
        mu_sum = mu_sum + Lambda_i(mu, V, i) * (K_i_tmp(:,2)+mu_i_tmp*mu_i_tmp(2));
    end
    out = mu_sum/Lambda(mu, V);
end

if corrected
    p.Lambda = @Lambda;
    p.mu_hat = @mu_hat;
    p.V_hat = @V_hat;
else
    p.Lambda = @Lambda_uncorrected;
    p.mu_hat = @mu_hat_uncorrected;
    p.V_hat = @V_hat_uncorrected;
end

end