%% Compute Gamma
% compute gamma, Equation (21), Freestone et al., 2011 NeuroImage
%%


clc
clear
close all

%% figure save path
% ~~~~~~~~~~~~~~~


figurePath = '../Figures/';

%% Generate data
% ~~~~~~~~~~~~~~~


SpaceMin = -10; SpaceMax = 10; NPoints = 501;
% create a cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);

%% define Gaussian basis functions
% ~~~~~~~~~~~~~~~
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
sigma = [2 0; 0 2]; % covariance matrix

% define gaussian basis functions
gaussians = zeros(NPoints, NPoints, nx);
for n = 1 : nx
    gaussians(:,:, n) = Define2DGaussian_AnisotropicKernel(mu(n, 1), mu(n, 2), sigma, NPoints, SpaceMin, SpaceMax);
end

%% plot all basis functions in the neural field
% ~~~~~~~~~~~~~~~


figure, clf, shg; imagesc(squeeze(sum(gaussians, 3))), colorbar; title('Guassian basis functions in the field');

%% Compute gamma - analytic
% ~~~~~~~~~~~~~~~


gamma_analytic = zeros(nx, nx);

% cycle through each pair of basis funciton to generate a nx times nx
% matrix - Gamma
for m = 1 : nx
    for n = 1 : nx
        gamma_analytic(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma); % analytically integral of two Gaussians in calculated as inner product of two Gaussians. See Appendix D, Freestone et al. NeuroImage 2011
    end
end

%% Compute gamma - numeric
% ~~~~~~~~~~~~~~~


gamma_numeric = zeros(nx, nx);

for m = 1 : nx
    for n = 1 : nx
        gamma_numeric(m, n) = spatialIntegral2Gaussians(X, Y, NPoints, mu(m, :)', sigma, mu(n, :)', sigma); % numerically the integral over space
    end
end

%% compare analytic and numeric results - residual
% ~~~~~~~~~~~~~~~


residual = gamma_analytic -gamma_numeric;

% visualise the residual
fig = figure; shg, clf;
subplot(3,1,1);
imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], gamma_analytic); colorbar; title('Analytic'); % plot analytic
subplot(3,1,2);
imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], gamma_numeric), colorbar; title('Numeric'); % plot numeric
subplot(3,1,3);
imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], residual), colorbar; title('Diff(analytic, numeric)');
suptitle('Compute Gamma - Compare analytic and numeric');

filename =[figurePath 'ComputeGamma_Check_AnalyticNumeric_SpatialRes_' num2str(NPoints) '.pdf'];
print(fig, '-dpdf', filename);
