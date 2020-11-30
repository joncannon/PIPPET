function p = PIPPET_stream_params(means_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices)

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

phibar_i_list = @(phibar, V) (phibar/V + p.e_means./p.e_vars)./(1/V + 1./p.e_vars);
K_i_list = @(V) 1./(1/V + 1./p.e_vars);
Lambda_i_list = @(phibar, V) p.e_lambdas .* gauss_distribution(phibar, p.e_means, p.e_vars+V);

p.Lambda_hat = @(phibar, V) lambda_0 + sum(Lambda_i_list(phibar, V));
p.phi_hat = @(phibar, V) (lambda_0*phibar + sum(Lambda_i_list(phibar, V) .* phibar_i_list(phibar,V)))/p.Lambda_hat(phibar,V);
p.V_hat= @(phibar_new, phibar_old, V) (lambda_0*(V+(phibar_old-phibar_new).^2) + sum(Lambda_i_list(phibar_old, V) .* (K_i_list(V) + (phibar_i_list(phibar_old, V)-phibar_new).^2)))/p.Lambda_hat(phibar_old,V);