function [integral] = spatialIntegral2Gaussians(gaussian1, gaussian2, stepSize)
%%
product = gaussian1 * gaussian2;
% product = conv2(gaussian1, gaussian2, 'same');

integral = sum(sum(product, 2), 1) * stepSize ^ 2;
end