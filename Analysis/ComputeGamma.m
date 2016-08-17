function gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu, sigma)
%% Compute Gamma
% To implement Equation (21). Freestone et al., 2011, NeuroImage
% Miao Cao

% ~~~~~~~~~~~~~~~
% input parameter list:
% parameters used to create a  2-D cortical surface
% SpaceMin - Edge of surface on the negative side
% SpaceMax - Edge of surface on the positive side
% NPoints - number of points along each dimension

% nx - number of Gaussian basis functions
% mu - centres of Gaussians
% sigma - sigma of Gaussians

% ~~~~~~~~~~~~~~~
% output parameter list:
% gamma - gamma matrix

%% Generate data
% ~~~~~~~~~~~~~~~
% create a 2-D surface

x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);

%% Compute gamma
% ~~~~~~~~~~~~~~~
gamma = zeros(nx, nx); % Gamma is a matrix with dimensions, nx * nx

for m = 1 : nx % compute pair-wise inner production of two Gaussians
    for n = 1 : nx
        gamma(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma);
    end
end

end