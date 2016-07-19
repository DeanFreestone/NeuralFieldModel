%% Analytic test of inner product
clc, clear, close all

%% Generate data
SpaceMin = -1; SpaceMax = 1; NPoints = 101;

x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);

[X,Y] = meshgrid(x, x);

% phi Gaussian basis function
%
var_phi = 0.02;
mu_phi = [0.5 0.5];
sigma_phi = [var_phi 0; 0 var_phi];
phi =Define2DGaussian_2(mu_phi(1), mu_phi(2), var_phi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_phi);


% psi Gaussian basis function
var_psi = 0.04;
mu_psi = [-0.5 -0.5];
sigma_psi = [var_psi 0; 0 var_psi];
psi = Define2DGaussian_2(mu_psi(1), mu_psi(2), var_psi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_psi);


%% mvnpdf
% phi Gaussian basis function
% mvnpdf
phi_mvnpdf = mvnpdf([X(:) Y(:)], mu_phi, sigma_phi);
phi_mvnpdf = reshape(phi_mvnpdf, NPoints, NPoints);

% psi Gaussian basis function
psi_mvnpdf = mvnpdf([X(:) Y(:)], mu_psi, sigma_psi);
psi_mvnpdf = reshape(psi_mvnpdf, NPoints, NPoints);

%% DefineGaussian3 - basic equation without coefficient
% phi_3
phi_3 =Define2DGaussian_3(mu_phi(1), mu_phi(2), var_phi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_phi);

% psi_3
psi_3 = Define2DGaussian_3(mu_psi(1), mu_psi(2), var_psi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_psi);

%% DefineGaussian1 - original
% phi_3
phi_1 =Define2DGaussian(mu_phi(1), mu_phi(2), var_phi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_phi);

% psi_3
psi_1 = Define2DGaussian(mu_psi(1), mu_psi(2), var_psi, 0, NPoints, SpaceMin, SpaceMax)/sqrt(2*pi*var_psi);


%% Plot
figure;
surf(X, Y, phi_3); hold on;
surf(X, Y, psi_3); title('phi and psi - 2DCreateGaussian-basic equation without coefficient'); colorbar;

figure;
surf(X, Y, phi_1); hold on;
surf(X, Y, psi_1); title('phi and psi - 2DCreateGaussian-original'); colorbar;


figure;
surf(X, Y, phi_mvnpdf); hold on;
surf(X, Y, psi_mvnpdf); title('phi and psi - mvnpdf'); colorbar;
%% compute 2 ways

phi_psi_conv2 = conv2(phi, psi, 'same');

%% Analytical check of convolution of two gaussians
%
phi_psi_1 = pi*var_psi*var_phi*Define2DGaussian(mu_phi(1)+mu_psi(1), mu_phi(2)+mu_psi(2), var_phi + var_psi, 0, NPoints, SpaceMin, SpaceMax) / (var_phi+var_psi);

phi_psi_2 = pi*var_psi*var_phi*Define2DGaussian_2(mu_phi(1)+mu_psi(1), mu_phi(2)+mu_psi(2), var_phi + var_psi, 0, NPoints, SpaceMin, SpaceMax) / (var_phi+var_psi);

phi_psi_3 = pi*var_psi*var_phi*Define2DGaussian_3(mu_phi(1)+mu_psi(1), mu_phi(2)+mu_psi(2), var_phi + var_psi, 0, NPoints, SpaceMin, SpaceMax) / (var_phi+var_psi);

figure, imagesc(phi_psi_conv2); colorbar; title('phi and psi - conv2');
figure, imagesc(phi_psi_1); colorbar; title('phi and psi - Analytic - Define2DGaussian-orginal');
figure, imagesc(phi_psi_2); colorbar; title('phi and psi - Analytic - Define2DGaussian-basic equation');

phi_psi_conv2(51,51)/phi_psi_2(51,51)
