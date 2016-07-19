%% Compute Gamma
% compute gamma, Equation (21), Freestone et al., 2011 NeuroImage
%%
clc, clear, close all
%% Generate data
SpaceMin = -10; SpaceMax = 10; NPoints = 1001;
% create a cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
%% define Gaussian basis functions
% define mu and sigma for gaussian basis functions

nx = 16; % number of basis functions

numRow = sqrt(nx); % number of gaussians for each colomn
numCol = nx / numRow; % number of columns

widthSpace = SpaceMax - SpaceMin;
widthCentre = widthSpace / (numCol*2);

mu = zeros(nx, 2); % centres of each gaussian
for m = 1 : numRow
    for n = 1 : numCol
        mu(n + numCol*(m-1), :) = [(SpaceMin - widthCentre + m*widthCentre*2) (SpaceMin - widthCentre + n*widthCentre*2)];
    end
end
sigma = [3 0; 0 3]; % covariance matrix

% define gaussian basis functions
guassians = zeros(NPoints, NPoints, nx);
for n = 1 : nx
    guassians(:,:, n) = Define2DGaussian_AnisotropicKernel(mu(n, 1), mu(n, 2), sigma, NPoints, SpaceMin, SpaceMax);
end
%% plot
figure, clf, shg; imagesc(squeeze(sum(guassians, 3))), colorbar;
%% Gamma
gamma = zeros(nx, nx);

for m = 1 : nx
    for n = 1 : nx
        gamma(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma);
    end
end
