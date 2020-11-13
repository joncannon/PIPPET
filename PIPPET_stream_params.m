function p = PIPPET_stream_params(events_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices)

p = struct();

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

p.e_means = [];
p.event_times = event_times;
p.lambda_0 = lambda_0;
p.highlight_event_indices = highlight_event_indices;
p.highlight_expectations = [];

for i=1:expected_cycles
    p.e_means = [p.e_means, events_unit + (i-1)*expected_period];
    p.highlight_expectations = [p.highlight_expectations, highlight_expectations];
end

p.e_vars = repmat(variance_unit, [1,expected_cycles]);
p.e_lambdas = repmat(lambda_unit, [1,expected_cycles]);

mu_i_list = @(mu, C) (mu/C + p.e_means./p.e_vars)./(1/C + 1./p.e_vars);
K_i_list = @(C) 1./(1/C + 1./p.e_vars);
Lambda_i_list = @(mu, C) p.e_lambdas .* gauss_distribution(mu, p.e_means, p.e_vars+C);
p.Lambda_bar = @(mu, C) lambda_0 + sum(Lambda_i_list(mu, C));
p.mu_bar = @(mu, C) (lambda_0*mu + sum(Lambda_i_list(mu, C) .* mu_i_list(mu,C)))/p.Lambda_bar(mu,C);
p.C_bar= @(mu,C) (lambda_0*C + sum(Lambda_i_list(mu, C) .* (K_i_list(C) + (mu_i_list(mu, C)-mu).^2)))/p.Lambda_bar(mu,C);