%% Compute Psi
% This script is to analytically and numerically check the function
% ComputePsi
% Miao Cao


clc
clear
close all

%% parameters and variables are pre-defined here
% ~~~~~~~~~~~~~~~



SpaceMin = -10; SpaceMax = 10; NPoints = 301;
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);

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

sigma_phi = [2 0; 0 2]; % covariance matrix of phi

sigma_psi = [0.6, 0.8, 2];

mu_psi = [1 1;
    1 1;
    1 1];

%% Compute gamma
% ~~~~~~~~~~~~~~~


Gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu_phi, sigma_phi); % compute gamma based on the function

%% Compute Psi - analytic
% ~~~~~~~~~~~~~~~
% now form the matrix
% these are the coefficients for the analytic convolution of psi and phi
% But, we haven't figure out covariance matrix here.



psi_phi_coefficient(1) = pi*sigma_psi(1)*sigma_phi(1, 1) / (sigma_psi(1)+sigma_phi(1, 1));
psi_phi_coefficient(2) = pi*sigma_psi(2)*sigma_phi(1, 1) / (sigma_psi(2)+sigma_phi(1, 1));
psi_phi_coefficient(3) = pi*sigma_psi(3)*sigma_phi(1, 1) / (sigma_psi(3)+sigma_phi(1, 1));

% compute the convolution between phi and psi
psi_phi_basis = zeros(length(theta), nx, NPoints, NPoints);
for m=1:length(theta)
    for n=1 : nx
        % these guys here are used with the LS algorithm for estimating
        % theta and xi
        mu = mu_phi(n, :) + mu_psi(m, :) + 2*mu_psi(m, :);
        psi_phi = psi_phi_coefficient(m)*Define2DGaussian_AnisotropicKernel(mu(1), mu(2), [sigma_psi(m) 0; 0 sigma_psi(m)]+sigma_phi, NPoints, SpaceMin, SpaceMax);
        
        psi_phi_basis(m, n, :, :) = psi_phi(:, :);
        %         theta_psi_phi_basis(nn,n,:) = theta(nn)*psi_phi_basis(nn,n,:);
    end
end

% Ts_invGamma_theta_phi_psi = Ts*(Gamma\squeeze(theta_psi_phi_basis(1,:,:) ...
%     + theta_psi_phi_basis(2,:,:) ...
%     + theta_psi_phi_basis(3,:,:)));

Ts_invGamma_phi_psi = zeros(length(theta), nx, NPoints, NPoints); % initialise the matrix of fields

inv_Gamma = inv(Gamma);

for m = 1 : length(theta) % cycle through each row of the matrix of fields
    
    fieldVector = squeeze(psi_phi_basis(m, :, :, :)); % a
    inv_Gamma_fieldVector = zeros(size(fieldVector));
    
    for p = 1 : nx
        for q = 1 : nx
            inv_Gamma_fieldVector(p, :, :) = inv_Gamma_fieldVector(p, :, :) + inv_Gamma(p, q) .* fieldVector(q, :, :);
        end
    end
    
    Ts_invGamma_phi_psi(m,:,:,:) = Ts * inv_Gamma_fieldVector;
end

% Ts_invGamma_phi_psi(1,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(1, :, :, :)));
% Ts_invGamma_phi_psi(2,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(2, :, :, :)));
% Ts_invGamma_phi_psi(3,:,:,:) = Ts*(Gamma\squeeze(psi_phi_basis(3, :, :, :)));

%% Compute Psi - numeric
% ~~~~~~~~~~~~~~~


% Compute convolution
psi_phi_basis_numeric = zeros(length(theta), nx, NPoints, NPoints);
for m=1:length(theta)
    for n=1 : nx
        % these guys here are used with the LS algorithm for estimating
        % theta and xi
        phi_gaussian = Define2DGaussian_AnisotropicKernel(mu_phi(n, 1), mu_phi(n, 2), sigma_phi, NPoints, SpaceMin, SpaceMax);
        psi_gaussian = Define2DGaussian_AnisotropicKernel(3*mu_psi(m, 1), 3*mu_psi(m, 2), [sigma_psi(m) 0; 0 sigma_psi(m)], NPoints, SpaceMin, SpaceMax);
        
        psi_phi_numeric(:, :) = conv2(phi_gaussian, psi_gaussian, 'same') .* stepSize^2;
        
        psi_phi_basis_numeric(m, n, :, :) = psi_phi_numeric(:, :);
    end
end

Ts_invGamma_phi_psi_numeric = zeros(length(theta), nx, NPoints, NPoints); % initialise the matrix of fields

inv_Gamma = inv(Gamma);

for n = 1 : length(theta) % cycle through each row of the matrix of fields
    
    fieldVector = squeeze(psi_phi_basis_numeric(n, :, :, :)); %
    inv_Gamma_fieldVector = zeros(size(fieldVector));
    
    for p = 1 : nx
        for q = 1 : nx
            inv_Gamma_fieldVector(p, :, :) = inv_Gamma_fieldVector(p, :, :) + inv_Gamma(p, q) .* fieldVector(q, :, :);
        end
    end
    
    Ts_invGamma_phi_psi_numeric(n,:,:,:) = Ts * inv_Gamma_fieldVector;
end
%% compare difference
% ~~~~~~~~~~~~~~~



for m=1:length(theta)
    for n=1 : nx
        figure, shg, clf;
        subplot(3,1,1);
        imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], squeeze(psi_phi_basis(m, n, :, :))); colorbar; title('analytic');
        subplot(3,1,2);
        imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], squeeze(psi_phi_basis_numeric(m, n, :, :))), colorbar; title('numeric');
        subplot(3,1,3);
        diffMat = squeeze(psi_phi_basis(m, n, :, :)) - squeeze(psi_phi_basis_numeric(m, n, :, :));
        imagesc([SpaceMin SpaceMax], [SpaceMax SpaceMin], diffMat), colorbar; title('Diff(analytic, numeric)');
        suptitle(['n:' num2str(n) ' m:' num2str(m)]);
    end
end
%% compare
xT = -10 : 0.1 : 10;
offset = 2;
mu_f = 0; mu_g = 1;
f = exp(-1*(xT - mu_f).^2);
g = 0.5*exp(-1*0.3*(xT - mu_g).^2);
g2 = 0.5*exp(-1*0.3*(xT - mu_g - offset).^2);
c1 = conv(f, g, 'same');
c2 = conv(f, g2, 'same');
figure;
subplot(2,1,1);
plot(xT, f); hold on; plot(xT, g); hold on; plot(xT, c1); hold on; plot(xT, c2);
subplot(2,1,2);
plot(xT+offset, c1); xlim([-10 10])
