% © Jonathan Cannon, MIT, 2020
% Creates parameters (including auxilliary functions) for a single event stream in the PIPPET model.
% Inputs:
%   means_unit:         one repetition of a pattern of expected event times mu_i
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
% Output:
%   p:   parameter set for stream

function p = PIPPET_stream_params(means_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, eta_e)

p = struct();

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

p.e_means = [];
p.event_times = event_times;
p.perceived_event_times = event_times + randn(size(event_times))*eta_e;
p.lambda_0 = lambda_0;
p.highlight_event_indices = highlight_event_indices;
p.highlight_expectations = [];

for i=1:expected_cycles
    p.e_means = [p.e_means, means_unit + (i-1)*expected_period];
    p.highlight_expectations = [p.highlight_expectations, highlight_expectations];
end

p.e_vars = repmat(variance_unit, [1,expected_cycles]);
p.e_lambdas = repmat(lambda_unit, [1,expected_cycles]);

mu_i_list = @(mu, V) (mu/V + p.e_means./p.e_vars)./(1/V + 1./p.e_vars);
K_i_list = @(V) 1./(1/V + 1./p.e_vars);
Lambda_i_list = @(mu, V) p.e_lambdas .* gauss_distribution(mu, p.e_means, p.e_vars+V);

p.Lambda = @(mu, V) lambda_0 + sum(Lambda_i_list(mu, V));
p.mu_hat = @(mu, V) (lambda_0*mu + sum(Lambda_i_list(mu, V) .* mu_i_list(mu,V)))/p.Lambda(mu,V);
p.V_hat= @(mu_new, mu_old, V) (lambda_0*(V+(mu_old-mu_new).^2) + sum(Lambda_i_list(mu_old, V) .* (K_i_list(V) + (mu_i_list(mu_old, V)-mu_new).^2)))/p.Lambda(mu_old,V);