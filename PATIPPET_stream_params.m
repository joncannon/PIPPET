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



xbar_i = @(xbar, Sigma, i) (inv(Sigma) + [1./p.e_vars(i), 0;0,0])\(Sigma\xbar + [p.e_means(i)./p.e_vars(i);0]);
K_i = @(Sigma, i) inv(inv(Sigma) + [1./p.e_vars(i), 0;0,0]);
Lambda_i = @(xbar, Sigma, i) p.e_lambdas(i) .* gauss_distribution(xbar(1), p.e_means(i), p.e_vars(i)+Sigma(1,1));






function Lambda_hat = Lambda_hat_uncorrected(xbar, Sigma)
    Lambda_hat = p.lambda_0;
    for i = 1:length(p.e_means)
         Lambda_hat = Lambda_hat + Lambda_i(xbar, Sigma, i);
    end
end

function Sigma_hat = Sigma_hat_uncorrected(xbar_new, xbar_old, Sigma)

    Sigma_sum = p.lambda_0 * (Sigma + (xbar_old-xbar_new)*(xbar_old-xbar_new)');
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar_old, Sigma, i);
        K_i_tmp = K_i(Sigma, i);
        Sigma_sum = Sigma_sum + Lambda_i(xbar_old, Sigma, i) * (K_i_tmp + (xbar_i_tmp-xbar_new)*(xbar_i_tmp-xbar_new)');
    end
    Sigma_hat = Sigma_sum/Lambda_hat(xbar_old, Sigma);
end
   
function x_hat = x_hat_uncorrected(xbar, Sigma)
    x_sum = p.lambda_0*xbar;
    for i = 1:length(p.e_means)
         x_sum = x_sum + Lambda_i(xbar, Sigma, i) * xbar_i(xbar, Sigma, i);
    end
    x_hat = x_sum/Lambda_hat(xbar, Sigma);
end







function Lambdahat = Lambda_hat(xbar, Sigma)
    Lambdahat = p.lambda_0 * xbar(2);
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar, Sigma, i);
        Lambdahat = Lambdahat + Lambda_i(xbar, Sigma, i)*xbar_i_tmp(2);
    end
end

function Sigmahat = Sigma_hat(xbar_new, xbar_old, Sigma)

    Sigma_sum = p.lambda_0 * (xbar_old(2)*(Sigma + (xbar_old-xbar_new)*(xbar_old-xbar_new)') + (xbar_old - xbar_new)*Sigma(2,:) +Sigma(:,2)*(xbar_old - xbar_new)');
    for i = 1:length(p.e_means)
        xbar_i_tmp = xbar_i(xbar_old, Sigma, i);
        K_i_tmp = K_i(Sigma, i);
        Sigma_sum = Sigma_sum + Lambda_i(xbar_old, Sigma, i) * (xbar_i_tmp(2)*(K_i_tmp + (xbar_i_tmp-xbar_new)*(xbar_i_tmp-xbar_new)') + (xbar_i_tmp - xbar_new)*K_i_tmp(2,:) + K_i_tmp(:,2)*(xbar_i_tmp - xbar_new)');
    end
    Sigmahat = Sigma_sum/Lambda_hat(xbar_old, Sigma);
end
   
    function xhat = x_hat(xbar, Sigma)
    x_sum = p.lambda_0*(Sigma(:,2)+xbar*xbar(2));
    for i = 1:length(p.e_means)
        K_i_tmp = K_i(Sigma, i);
        xbar_i_tmp = xbar_i(xbar, Sigma, i);
        x_sum = x_sum + Lambda_i(xbar, Sigma, i) * (K_i_tmp(:,2)+xbar_i_tmp*xbar_i_tmp(2));
    end
    xhat = x_sum/Lambda_hat(xbar, Sigma);
end

if corrected
    p.Lambda_hat = @Lambda_hat;
    p.x_hat = @x_hat;
    p.Sigma_hat = @Sigma_hat;
else
    p.Lambda_hat = @Lambda_hat_uncorrected;
    p.x_hat = @x_hat_uncorrected;
    p.Sigma_hat = @Sigma_hat_uncorrected;
end

end