function out = expectation_func(x_list, e_means, e_vars, e_taus, tau_0)

gauss_distribution = @(x, mean, v) exp(-.5 * ((x - mean).^ 2) ./ v)./ (sqrt(2*pi*v));

out = tau_0 + zeros(size(x_list));
    
for i = 1:length(e_means)
    out = out + e_taus(i)*gauss_distribution(x_list, e_means(i), e_vars(i));
end