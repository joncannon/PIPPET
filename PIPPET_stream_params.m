% © Jonathan Cannon, MIT, 2020
% Creates parameters (including auxilliary functions) for a single event stream in the PIPPET model.
% Inputs:
%   means_unit:         one repetition of a pattern of expected event times phi_i
%   variance_unit:      one repetition of a pattern of expected event timing variances v_i
%   tau_unit:           one repetition of a pattern of event expectation strengths tau_i
%   tau_0
%   expected_cycles:    number of cycles that the patterns repeat to create
%                           temporal expectation template
%   expected_period:    period with which the patterns repeat
%   event_times:        time points at which events actually occur in this
%                           stream
%   highlight_expectations:     Display weights for lines marking expected timepoints
%   highlight_event_indices:    Display weights for lines marking event times
% Output:
%   p:   parameter set for stream

function p = PIPPET_stream_params(means_unit, variance_unit, tau_unit, tau_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, eta_e)

p = struct();

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

p.e_means = [];
p.event_times = event_times;
p.perceived_event_times = event_times + randn(size(event_times))*eta_e;
p.tau_0 = tau_0;
p.highlight_event_indices = highlight_event_indices;
p.highlight_expectations = [];

for i=1:expected_cycles
    p.e_means = [p.e_means, means_unit + (i-1)*expected_period];
    p.highlight_expectations = [p.highlight_expectations, highlight_expectations];
end

p.e_vars = repmat(variance_unit, [1,expected_cycles]);
p.e_taus = repmat(tau_unit, [1,expected_cycles]);

phibar_i_list = @(phibar, V) (phibar/V + p.e_means./p.e_vars)./(1/V + 1./p.e_vars);
K_i_list = @(V) 1./(1/V + 1./p.e_vars);
T_i_list = @(phibar, V) p.e_taus .* gauss_distribution(phibar, p.e_means, p.e_vars+V);

p.T_hat = @(phibar, V) tau_0 + sum(T_i_list(phibar, V));
p.phi_hat = @(phibar, V) (tau_0*phibar + sum(T_i_list(phibar, V) .* phibar_i_list(phibar,V)))/p.T_hat(phibar,V);
p.V_hat= @(phibar_new, phibar_old, V) (tau_0*(V+(phibar_old-phibar_new).^2) + sum(T_i_list(phibar_old, V) .* (K_i_list(V) + (phibar_i_list(phibar_old, V)-phibar_new).^2)))/p.T_hat(phibar_old,V);