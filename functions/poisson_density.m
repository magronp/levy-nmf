function [p] = poisson_density(x,lambda)

p = lambda.^x ./ gamma(x+1) .* exp(-lambda);

end