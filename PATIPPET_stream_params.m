% © Jonathan Cannon, MIT, 2020
% Creates parameters (including auxilliary functions) for a single event stream in the PATIPPET model.
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
%   corrected:          whether PATIPPET model is corrected to normalize
%                           event rate over phase (rather than over time) 
% Output:
%   p:   parameter set for stream

function p = PATIPPET_stream_params(means_unit, variance_unit, tau_unit, tau_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices, corrected)

p = struct();

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

p.e_means = [];
p.event_times = event_times;
p.tau_0 = tau_0;
p.highlight_event_indices = highlight_event_indices;
p.highlight_expectations = [];

for i=1:expected_cycles
    p.e_means = [p.e_means, means_unit + (i-1)*expected_period];
    p.highlight_expectations = [p.highlight_expectations, highlight_expectations];
end

p.e_vars = repmat(variance_unit, [1,expected_cycles]);
p.e_taus = repmat(tau_unit, [1,expected_cycles]);



xbar_i = @(xbar, Sigma, i) (inv(Sigma) + [1./p.e_vars(i), 0;0,0])\(Sigma\xbar + [p.e_means(i)./p.e_vars(i);0]);
K_i = @(Sigma, i) inv(inv(Sigma) + [1./p.e_vars(i), 0;0,0]);
T_i = @(xbar, Sigma, i) p.e_taus(i) .* gauss_distribution(xbar(1), p.e_means(i), p.e_vars(i)+Sigma(1,1));


function out = T_hat_uncorrected(xbar, Sigma)
    out = p.tau_0;
    for i = 1:length(p.e_means)
         out = out + T_i(xbar, Sigma, i);
    end
end

function out = Sigma_hat_uncorrected(xbar_new, xbar_old, Sigma)

    Sigma_sum = p.tau_0 * (Sigma + (xbar_old-xbar_new)*(xbar_old-xbar_new)');
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar_old, Sigma, i);
        K_i_tmp = K_i(Sigma, i);
        Sigma_sum = Sigma_sum + T_i(xbar_old, Sigma, i) * (K_i_tmp + (xbar_i_tmp-xbar_new)*(xbar_i_tmp-xbar_new)');
    end
    out = Sigma_sum/T_hat(xbar_old, Sigma);
end
   
function out = x_hat_uncorrected(xbar, Sigma)
    x_sum = p.tau_0*xbar;
    for i = 1:length(p.e_means)
         x_sum = x_sum + T_i(xbar, Sigma, i) * xbar_i(xbar, Sigma, i);
    end
    out = x_sum/T_hat(xbar, Sigma);
end


function out = T_hat(xbar, Sigma)
    out = p.tau_0 * xbar(2);
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar, Sigma, i);
        out = out + T_i(xbar, Sigma, i)*xbar_i_tmp(2);
    end
end

function out = Sigma_hat(xbar_new, xbar_old, Sigma)

    Sigma_sum = p.tau_0 * (xbar_old(2)*(Sigma + (xbar_old-xbar_new)*(xbar_old-xbar_new)') + (xbar_old - xbar_new)*Sigma(2,:) +Sigma(:,2)*(xbar_old - xbar_new)');
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar_old, Sigma, i);
        K_i_tmp = K_i(Sigma, i);
        Sigma_sum = Sigma_sum + T_i(xbar_old, Sigma, i) * (xbar_i_tmp(2)*(K_i_tmp + (xbar_i_tmp-xbar_new)*(xbar_i_tmp-xbar_new)') + (xbar_i_tmp - xbar_new)*K_i_tmp(2,:) + K_i_tmp(:,2)*(xbar_i_tmp - xbar_new)');
    end
    out = Sigma_sum/T_hat(xbar_old, Sigma);
end
   
function out = x_hat(xbar, Sigma)
    x_sum = p.tau_0*(Sigma(:,2)+xbar*xbar(2));
    for i = 1:length(p.e_means)
        K_i_tmp = K_i(Sigma, i);
        xbar_i_tmp = xbar_i(xbar, Sigma, i);
        x_sum = x_sum + T_i(xbar, Sigma, i) * (K_i_tmp(:,2)+xbar_i_tmp*xbar_i_tmp(2));
    end
    out = x_sum/T_hat(xbar, Sigma);
end

if corrected
    p.T_hat = @T_hat;
    p.x_hat = @x_hat;
    p.Sigma_hat = @Sigma_hat;
else
    p.T_hat = @T_hat_uncorrected;
    p.x_hat = @x_hat_uncorrected;
    p.Sigma_hat = @Sigma_hat_uncorrected;
end

end