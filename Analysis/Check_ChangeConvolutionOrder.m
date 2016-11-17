%% Analytic test of inner product
% ~~~~~~~~~~~~~~~



clc
clear
close all
%% Generate data
% ~~~~~~~~~~~~~~~




SpaceMin = -4; SpaceMax = 4; NPoints = 401;
% random data, fData
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
fData = randn(NPoints);

% Define phi Gaussian basis function
%
mu_phi = [0.5 0.5];                   % center of basis function
sigma_phi = [0.02 0; 0 0.02];           % variance-covariance matrix of phi

phi = Define2DGaussian(mu_phi(1), mu_phi(2), sigma_phi(1,1), 0, NPoints, SpaceMin, SpaceMax);
% phi = Define2DGaussian(mu_phi(1), mu_phi(2), sigma_phi(1,1), 0, NPoints, SpaceMin, SpaceMax);

% Define psi Gaussian kernel
mu_psi = [-.5 -.5];              % centre of basis function
sigma_psi = [0.04 0; 0 0.04];        % variance-covariance matrix of phi
psi = Define2DGaussian(mu_psi(1), mu_psi(2), sigma_psi(1,1), 0, NPoints, SpaceMin, SpaceMax);
% psi = Define2DGaussian(mu_psi(1), mu_psi(2), sigma_psi(1,1), 0, NPoints, SpaceMin, SpaceMax);


% psi = mvnpdf([X(:) Y(:)], mu_psi, sigma_psi);
% psi = reshape(psi, len_x, len_x);


%% compute 2 ways
% ~~~~~~~~~~~~~~~



% E1 Phi@(Psi*fData), @ is defined as inner product of two functions (matrices)
convE1 = conv2(psi, fData, 'same');
E1 = sum(sum(phi.*convE1, 2), 1);

% E2 (Phi*Psi)@fData
convE2 = conv2(phi, psi, 'same');
E2 = sum(sum(convE2.*fData, 2), 1);

conv2_convPhiPsi = convE2 * stepSize^2;
%% Plot
% ~~~~~~~~~~~~~~~



figure, surf(X, Y, convE1); colorbar; title('Psi convolves with random data');
