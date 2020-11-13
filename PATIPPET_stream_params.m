function p = PATIPPET_stream_params(events_unit, variance_unit, lambda_unit, lambda_0, expected_cycles, expected_period, event_times, highlight_expectations, highlight_event_indices)

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



mu_i = @(mu, C, i) (inv(C) + [1./p.e_vars(i), 0;0,0])\(C\mu + [p.e_means(i)./p.e_vars(i);0]);
K_i = @(C, i) inv(inv(C) + [1./p.e_vars(i), 0;0,0]);
Lambda_i = @(mu, C, i) p.e_lambdas(i) .* gauss_distribution(mu(1), p.e_means(i), p.e_vars(i)+C(1,1));

function LS = Lambda_bar(mu, C)
    LS = p.lambda_0;
    for i = 1:length(p.e_means)
         LS = LS + Lambda_i(mu, C, i);
    end
end

function CS = C_bar(mu, C)
    S = p.lambda_0*C;
    for i = 1:length(p.e_means)
        S = S + Lambda_i(mu, C, i) * (K_i(C, i) + (mu_i(mu, C, i)-mu)*(mu_i(mu, C, i)-mu)');
    end
    CS = S/Lambda_bar(mu, C);
end
   
function MS = mu_bar(mu, C)
    S = p.lambda_0*mu;
    for i = 1:length(p.e_means)
         S = S + Lambda_i(mu, C, i) * mu_i(mu, C, i);
    end
    MS = S/Lambda_bar(mu, C);
end

p.Lambda_bar = @Lambda_bar;
p.mu_bar = @mu_bar;
p.C_bar = @C_bar;

end