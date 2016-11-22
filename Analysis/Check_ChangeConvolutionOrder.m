%% Analytic check of swap order of convolution and inner-product
% ~~~~~~~~~~~~~~~



clc


clear

close all

%% Generate data
% ~~~~~~~~~~~~~~~




SpaceMin = -10; SpaceMax = 10; NPoints = 301;

% random data, fData
x = linspace(SpaceMin, SpaceMax, NPoints);

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x);

fData = randn(NPoints);

% Define phi Gaussian basis function
%
mu_phi = [0 0];                   % center of basis function

sigma_phi = [1 0.2; 0.2 1];           % variance-covariance matrix of phi


phi = Define2DGaussian_AnisotropicKernel(mu_phi(1), mu_phi(2), sigma_phi, NPoints, SpaceMin, SpaceMax);
% phi = Define2DGaussian_3(mu_phi(1), mu_phi(2), sigma_phi(1,1), 0, NPoints, SpaceMin, SpaceMax);


% Define psi Gaussian kernel
mu_psi = [0 0];              % centre of basis function

sigma_psi = [2 1; 1 2];        % variance-covariance matrix of phi

psi = Define2DGaussian_AnisotropicKernel(mu_psi(1), mu_psi(2), sigma_psi, NPoints, SpaceMin, SpaceMax);
% psi = Define2DGaussian_3(mu_psi(1), mu_psi(2), sigma_psi(1,1), 0, NPoints, SpaceMin, SpaceMax);


% psi = mvnpdf([X(:) Y(:)], mu_psi, sigma_psi);
% psi = reshape(psi, len_x, len_x);


%% compute 2 ways
% ~~~~~~~~~~~~~~~



% E1 Phi@(Psi*fData), @ is defined as inner product of two functions (matrices)

% convolution
convE1 = conv2(psi, fData, 'same');

% integral
E1 = sum(sum(phi.*convE1*stepSize^2, 2), 1);


% E2 (Phi*Psi)@fData

% convolution
convE2 = conv2(phi, psi, 'same');

% integral
E2 = sum(sum(convE2.*fData*stepSize^2, 2), 1);

residual = E1 - E2
