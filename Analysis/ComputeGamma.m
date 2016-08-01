function gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu, sigma)
%%
% SpaceMin, SpaceMax, NPoints, cortical surface
% nx, number of basis functions
% sigma, covariance matrix of gaussian basis function
%% Generate data
% create a cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
%% Compute gamma
% Equation 
gamma = zeros(nx, nx);

for m = 1 : nx
    for n = 1 : nx
        gamma(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma);
    end
end
end