function [v_tplus1, v_t] = ReducedModel_ComputeFieldVtPlus1(x_t, nx, sigma_phi, theta, nTheta, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints)
%% Reduced model
% Compute field at time (T+1)
% To implement Equation (25), Freestone et al., 2011, NeuroImage
% Miao Cao


clc
clear
close all

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
% SpaceMin = -10; SpaceMax = 10; NPoints = 501;

x = linspace(SpaceMin, SpaceMax, NPoints);

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x);

%% Parameters
% ~~~~~~~~~~~~~~~


Ts = 0.0001; % time step

% ~~~~~~~~~~~~~~~
% field basis functions
% nx = 16; number of Gaussian basis functions

% sigma_phi = 0.8; % variance of gaussian basis functions

mu_phi = []; % centres of Gaussian basis functions. But for now let's leave it empty for function ComputePsi to generate centres that are uniformly distributed on the surface.


% ~~~~~~~~~~~~~~~
% connectivity kernel
% theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

% nTheta = 3; % number of connectivity kernel basis functions

mu_psi = [1 1;
    1 1;
    1 1]; % centres of Gaussian basis functions of connectivity kernel
% vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

% ~~~~~~~~~~~~~~~
% parameters for firing rate function
slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

%% Create Phi basis functions
% ~~~~~~~~~~~~~~~


phi_basisFunctions = CreatePhiBasisFunctions(SpaceMin, SpaceMax, NPoints, nx, mu_phi, sigma_phi); % use this function to create Gaussian basis functions

%% Compute Psi
% ~~~~~~~~~~~~~~~


psi = ComputePsi(SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi); % compute psi with function ComputePsi

%% Firing rate function
% ~~~~~~~~~~~~~~~


% Firing rate function. Equation (4) Freestone et al., 2011, NeuroImage

% x_t = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.
phi_fields = zeros(size(phi_basisFunctions));

for m = 1 : nx % to compute phi * x(t)
    phi_fields(:,:, m) = phi_basisFunctions(:,:, m) * x_t(m);
end
v_t = sum(phi_fields, 3); % v, mean membrane potential field, at time t.

% firing rate function
firingRate_phi_xt = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - v_t))); % firing rate sigmoidal function, field

%% integral over 2-D space
% ~~~~~~~~~~~~~~~


% compute x(T+1)
ingtegralProduct = zeros(nx, nTheta);
for pNX = 1 : nx
    for qNTheta = 1 : nTheta
        product_psi_firingRate = squeeze(psi(qNTheta, pNX, :, :)) * firingRate_phi_xt;
        ingtegralProduct(pNX, qNTheta) = sum(sum(product_psi_firingRate * stepSize^2, 2), 1);
    end
end
x_tplus1 = ingtegralProduct * theta; % finally times theta (vector) and get x(t+1)

%% Compute field at time (T+1)
% ~~~~~~~~~~~~~~~


for m = 1 : nx % to compute phi * x(t)
    fields_tplus1(:,:, m) = phi_basisFunctions(:,:, m) * x_tplus1(m);
end
v_tplus1 = squeeze(sum(fields_tplus1, 3));

end
