%% Compute Gamma
% compute gamma, Equation (21), Freestone et al., 2011 NeuroImage
%%
clc, clear, close all
%% Generate data
SpaceMin = -4; SpaceMax = 4; NPoints = 401;
% create a cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
%% define Gaussian basis functions
% define mu and sigma for gaussian basis functions

nx = 16; % number of basis functions

mu = [-3 3; -3 1; -3 -1; -3 -3;
    -1 3; -1 1; -1 -1; -1 -3;
    1 3; 1 1; 1 -1; 1 -3;
    3 3; 3 1; 3 -1; 3 -3];

sigma = [0.2 0; 0 0.2];

% define gaussian basis functions
guassians = zeros(NPoints, NPoints, nx);
for n = 1 : nx
    guassians(:,:, n) = Define2DGaussian_AnisotropicKernel(mu(n, 1), mu(n, 2), sigma, NPoints, SpaceMin, SpaceMax);
end

% phi = Define2DGaussian(mu_phi(1), mu_phi(2), sigma_phi(1,1), 0, NPoints, SpaceMin, SpaceMax);
%% Gamma
gamma = zeros(nx, nx);

for m = 1 : nx
    for n = 1 : nx
        gamma(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma);
    end
end
