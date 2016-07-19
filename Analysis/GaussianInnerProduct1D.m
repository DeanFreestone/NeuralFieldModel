%% Analytic test of inner product
clc, clear, close all
%% Generate data
SpaceMin = -1; SpaceMax = 1; NPoints = 100;
% random data, fData
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
len_x = length(x);
fData = randn(1, len_x);

% phi Gaussian basis function
%
mu_phi = 0.5; var_phi = 0.04;
phi = normpdf(x, mu_phi, var_phi);
% phi = mvnpdf([X(:) Y(:)], mu_phi, sigma_phi);
% phi = reshape(phi, len_x, len_x);

% psi Gaussian kernel
mu_psi = -0.5; var_psi = 0.04;
psi = normpdf(x, mu_psi, var_psi);
% psi = mvnpdf([X(:) Y(:)], mu_psi, sigma_psi);
% psi = reshape(psi, len_x, len_x);

% plot phi and psi
% figure; plot(x, phi); hold on; plot(x, psi);

%% compute 2 ways or change the order of computing

% E1 Phi@(Psi*fData), @ is defined as inner product of two functions (matrices)
convE1 = conv2(psi, fData, 'same');
E1 = sum(phi.*convE1*stepSize, 2);

% E2 (Phi*Psi)@fData
convE2 = conv2(phi, psi, 'same');
E2 = sum(convE2.*fData*stepSize, 2);

% Plot
figure, plot(x, fData); figure, plot(x, phi); figure, plot(x, psi);
figure; plot(x, convE1); figure; plot(x, convE2*stepSize);
%% analytical check of convolution of two gaussians (1-D)
mu = (mu_phi + mu_psi); %

% From Dean's paper
coefficient = sqrt((pi*var_phi*var_psi)/(var_phi + var_psi));
exponential = exp(-1*((x - mu).*(x - mu))/(var_phi + var_psi));

% From a paper, convolution of 2 uni-variate gaussians
% coefficient = sqrt(1/(2*pi*(var_phi + var_psi)));
% exponential = exp(-1*((x - mu).*(x - mu))/(2*(var_phi + var_psi)));
convE2_equivalent = coefficient*exponential;


% plot
figure, plot(x, convE2*stepSize); title('phi-psi - conv2');
figure, plot(x, convE2_equivalent); title('phi-psi - analytic');
