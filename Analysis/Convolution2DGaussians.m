%% Analytic test of inner product
clc, clear, close all
%% Generate data
SpaceMin = -4; SpaceMax = 4; NPoints = 401;
% random data, fData
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
fData = randn(NPoints);

% define phi Gaussian basis function
%
mu_phi = [0.5 0.5];                   % center of basis function
sigma_phi = [0.02 0; 0 0.02];           % variance-covariance matrix of phi
phi =Define2DGaussian_AnisotropicKernel(mu_phi(1), mu_phi(2), sigma_phi, NPoints, SpaceMin, SpaceMax);
% phi = Define2DGaussian(mu_phi(1), mu_phi(2), sigma_phi(1,1), 0, NPoints, SpaceMin, SpaceMax);


% define psi Gaussian basis function
mu_psi = [-.5 -.5];              % centre of basis function
sigma_psi = [0.04 0; 0 0.04];        % variance-covariance matrix of phi
psi = Define2DGaussian_AnisotropicKernel(mu_psi(1), mu_psi(2), sigma_psi, NPoints, SpaceMin, SpaceMax);
% psi = Define2DGaussian(mu_psi(1), mu_psi(2), sigma_psi(1,1), 0, NPoints, SpaceMin, SpaceMax);
%% Analytic check of convolution of two gaussians
% r is the location vector, specifically a row vector (consistent with the equation in Dean's paper)
mu = (mu_phi + mu_psi)';
var_phi = sigma_phi(1,1); var_psi = sigma_psi(1,1);
% if two gaussian centres are symmetric, addition of two centres will be
% (0, 0). Therefore, mu=0.

exponential = zeros(NPoints, NPoints);
CovMat = (sigma_psi + sigma_phi);
for m = 1 : NPoints
    for n = 1 : NPoints
        r = [X(m, n); Y(m, n)];
        
        %         exponent(i, j) = exp(-1*((r - mu)'*(r-mu))/(var_phi + var_psi));
        exponential(m, n) = exp(-((r - mu)'/CovMat*(r-mu)));
    end
end

coefficient = (pi*var_phi*var_psi)/(var_phi + var_psi);
convE2_equivalent = coefficient*exponential;                       % analytic solution

figure, imagesc(conv2_convPhiPsi); colorbar; title('phi-psi - conv2');
figure, imagesc(convE2_equivalent); colorbar; title('phi-psi - analytic');


%% plot the error
%
figure, imagesc(conv2_convPhiPsi - convE2_equivalent); colorbar; title('Diff(Con2 - analytic)');

% analytic
figure, surf(X, Y, convE2_equivalent); colorbar; title('analytic conv(phi,psi)');

% numerical Conv2
figure, surf(X, Y, conv2_convPhiPsi); colorbar; title('numerical (conv2) conv(phi,psi)');
