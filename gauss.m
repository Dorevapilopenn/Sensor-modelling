function f = gauss(x, x_max, sigma)
% Gaussian function with standard deviation sigma, centered at x_max

    f = exp(-0.5 * ((x - x_max) / sigma).^2);
end