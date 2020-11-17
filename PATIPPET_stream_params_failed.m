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



mu_i = @(mu, C, i) (inv(C) + [1./p.e_vars(i), 0;0,0])\(C\mu + [p.e_means(i)./p.e_vars(i);0]);
K_i = @(C, i) inv(inv(C) + [1./p.e_vars(i), 0;0,0]);
Lambda_i = @(mu, C, i) p.e_lambdas(i) .* gauss_distribution(mu(1), p.e_means(i), p.e_vars(i)+C(1,1));

third_moment_operator = @(mu, C) [C(1,1), 2*C(1,2); 2*C(1,2), 3*C(2,2)]*mu(2) + [2*C(1,2), C(2,2); C(2,2), 0]*mu(1) + mu(2)*mu*mu';




function LS = Lambda_bar(mu, C)
    LS = p.lambda_0 * mu(2);
    for i = 1:length(p.e_means)
        mu_i_tmp = mu_i(mu, C, i);
        LS = LS + Lambda_i(mu, C, i)*mu_i_tmp(2);
    end
end

function PS = Pi_bar(mu, C)
    S = p.lambda_0*third_moment_operator(mu, C);
    for i = 1:length(p.e_means)
        mu_i_tmp = mu_i(mu, C, i);
        K_i_tmp = K_i(C, i);
        S = S + Lambda_i(mu, C, i) * third_moment_operator(mu_i_tmp, K_i_tmp);
    end
    PS = S/Lambda_bar(mu, C);
end
   
function MS = mu_bar(mu, C)
    S = p.lambda_0*(C(:,2)+mu*mu(2));
    for i = 1:length(p.e_means)
        K_i_tmp = K_i(C, i);
        mu_i_tmp = mu_i(mu, C, i);
        S = S + Lambda_i(mu, C, i) * (K_i_tmp(:,2)+mu_i_tmp*mu_i_tmp(2));
    end
    MS = S/Lambda_bar(mu, C);
end

p.Lambda_bar = @Lambda_bar;
p.mu_bar = @mu_bar;
p.Pi_bar = @Pi_bar;
end
