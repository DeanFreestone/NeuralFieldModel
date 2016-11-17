%% Compute state vector x(t+1)
% To implement Equation (25), Freestone et al., 2011, NeuroImage
% Miao Cao


clc
clear
close all

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 301;

x = linspace(SpaceMin, SpaceMax, NPoints);

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x);

%% Parameters used to compute psi
% ~~~~~~~~~~~~~~~

tau = 0.01; % synaptic time constant

Ts = 0.0001; % time step

ks = 1- Ts*(1/tau); % time constant parameter

nx = 121; % number of Gaussian basis functions

theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

nTheta = 3; % number of connectivity kernel basis functions

% ~~~~~~~~~~~~~~~
% parameters for firing rate function

slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

% ~~~~~~~~~~~~~~~
mu_phi = []; % centres of Gaussian basis functions. But for now let's leave it empty for function ComputePsi to generate centres that are uniformly distributed on the surface.

sigma_phi = [2 0; 0 2]; % variance-covariance matrix of Gaussian basis function of field decomposition

% centres of Gaussian basis functions of connectivity kernel

mu_psi = [0 0; 0 0; 0 0]; % centres of basis functions of connectivity kernel

vector_Sigma_Psi = [0.6 0; 0.8 0; 2 0]; % width of Gaussian basis functions of connectivity kernel

%% Create Phi basis functions
% ~~~~~~~~~~~~~~~


phi_basisFunctions = CreatePhiBasisFunctions(SpaceMin, SpaceMax, NPoints, nx, mu_phi, sigma_phi); % use this function to create Gaussian basis functions

%% Compute Psi
% ~~~~~~~~~~~~~~~


psi = ComputePsi(SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi); % compute psi with function ComputePsi

%% Firing rate function
% ~~~~~~~~~~~~~~~


% Firing rate function. Equation (4) Freestone et al., 2011, NeuroImage

Xt = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.
phi_fields = zeros(size(phi_basisFunctions));

for m = 1 : nx % to compute phi * x(t)
    
    phi_fields(:,:, m) = phi_basisFunctions(:,:, m) * Xt(m);
    
end
v_t = sum(phi_fields, 3); % v, mean membrane potential field, at time t.

% firing rate function
firingRate_phi_xt = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - v_t))); % firing rate sigmoidal function, field

%% integrate over 2-D space
% ~~~~~~~~~~~~~~~



ingtegralProduct = zeros(nx, nTheta);

for pNX = 1 : nx
    for qNTheta = 1 : nTheta
        
        product_psi_firingRate = squeeze(psi(qNTheta, pNX, :, :)) * firingRate_phi_xt; % Psi times firing(Xt)
        
        ingtegralProduct(pNX, qNTheta) = sum(sum(product_psi_firingRate * stepSize^2, 2), 1);
        
    end
end

XtPlus1 = ks * Xt + ingtegralProduct * theta; % finally times theta (vector) and get x(t+1)
