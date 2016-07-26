%% Compute Psi
clc, clear, close all
%% parameters and variables are pre-defined here
SpaceMin = -10; SpaceMax = 10; NPoints = 401;
Ts = 0.0001;
nx = 16; % number of states basis functions
theta = [1, 1.2, 1.5]; % number of connectivity kernel basis functions

%
numRow = sqrt(nx); % number of gaussians for each colomn
numCol = nx / numRow; % number of columns

widthSpace = SpaceMax - SpaceMin;
widthCentre = widthSpace / (numCol*2);

mu_phi = zeros(nx, 2); % centres of each gaussian
for m = 1 : numRow
    for n = 1 : numCol
        mu_phi(n + numCol*(m-1), :) = [(SpaceMin - widthCentre + m*widthCentre*2) (SpaceMin - widthCentre + n*widthCentre*2)];
    end
end

sigma_phi = [4 0; 0 4]; % covariance matrix of phi

sigma_psi = [2, 3, 4];

mu_psi = [1 1;
    1 1;
    1 1];
%% Compute gamma
Gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu_phi, sigma_phi); % compute gamma based on the function
%% Compute Psi - analytic
% now form the matrix
% these are the coefficients for the analytic convolution of psi and phi
% But, we haven't figure out covariance matrix here.
psi_phi_coefficient(1) = pi*sigma_psi(1)^2*sigma_phi(1, 1)^2 / (sigma_psi(1)^2+sigma_phi(1, 1)^2);
psi_phi_coefficient(2) = pi*sigma_psi(2)^2*sigma_phi(1, 1)^2 ./ (sigma_psi(2)^2+sigma_phi(1, 1)^2);
psi_phi_coefficient(3) = pi*sigma_psi(3)^2*sigma_phi(1, 1)^2 ./ (sigma_psi(3)^2+sigma_phi(1, 1)^2);

% compute the convolution between phi and psi

for n=1 : nx
    for m=1:length(theta)
        % these guys here are used with the LS algorithm for estimating
        % theta and xi
        mu = mu_phi(n, :) + mu_psi(m, :) + 2*mu_psi(m, :);
        psi_phi = psi_phi_coefficient(m)*Define2DGaussian_AnisotropicKernel(mu(1), mu(2), [sigma_psi(m)^2 0; 0 sigma_psi(m)^2]+sigma_phi^2, NPoints, SpaceMin, SpaceMax);
        
        psi_phi_basis(m, n, :, :) = psi_phi(:, :);
        %         theta_psi_phi_basis(nn,n,:) = theta(nn)*psi_phi_basis(nn,n,:);
    end
end

% Ts_invGamma_theta_phi_psi = Ts*(Gamma\squeeze(theta_psi_phi_basis(1,:,:) ...
%     + theta_psi_phi_basis(2,:,:) ...
%     + theta_psi_phi_basis(3,:,:)));

Ts_invGamma_phi_psi(1,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(1, :, :, :)));
Ts_invGamma_phi_psi(2,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(2, :, :, :)));
Ts_invGamma_phi_psi(3,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(3, :, :, :)));
%% Compute Psi - numeric



