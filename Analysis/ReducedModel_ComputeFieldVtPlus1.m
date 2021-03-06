function [VtPlus1, Vt, XtPlus1] = ReducedModel_ComputeFieldVtPlus1(Xt, tau, Ts, nx, sigma_phi, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints)
%% Reduced model
% Compute neural field (post-synaptic membrane potential) at time (T+1)
% Miao Cao


% Parameter list:
% Xt - State vector at time T
% nx - Number of field basis functions (field decomposition)
% sigma_phi - width of field basis functions (field decomposition)
% theta - scale of basis functions of connectivity kernel
% nTheta - number of basis functions of connectivity kernel
% mu_psi - widths of basis functions of connectivity kernel
% vector_Sigma_Psi - a vector of widiths of Gaussian basis functions of
% connvectivity kernel (spatial decomposition)
% SpaceMin - the negative edge of cortical surface/neural field
% SpaceMax - the posive edge of cortical surface/neural field
% NPoints - number of points in each row or column of cortical
% surface/neural field (spatial resolution)

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
% a 2-D surface

x = linspace(SpaceMin, SpaceMax, NPoints); % coordinates of discretisation of the surface

stepSize = x(2)-x(1); % distance between two adjacent points in the 2-D surface

[X, Y] = meshgrid(x, x); % x and y coordinates of the 2-D surface


%% Parameters
% ~~~~~~~~~~~~~~~


% tau = 0.01; % synaptic time constant

% Ts = 0.0001; % time step

ks = 1- Ts*(1/tau); % time constant parameter

% ~~~~~~~~~~~~~~~
% field basis functions
% nx = 16; number of Gaussian basis functions

% sigma_phi = 0.8; % variance of gaussian basis functions

mu_phi = []; % centres of Gaussian basis functions. But for now let's leave it empty for function ComputePsi to generate centres that are uniformly distributed on the surface.


% ~~~~~~~~~~~~~~~
% connectivity kernel
% theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

% nTheta = 3; % number of connectivity kernel basis functions

% vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

% ~~~~~~~~~~~~~~~
% parameters for firing rate function
slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

%% Create Phi basis functions
% ~~~~~~~~~~~~~~~


[phi_basisFunctions, mu_phi] = CreatePhiBasisFunctions(SpaceMin, SpaceMax, NPoints, nx, mu_phi, sigma_phi); % use this function to create Gaussian basis functions

%% Compute Psi
% ~~~~~~~~~~~~~~~


psi = ComputePsi(SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi); % compute PS with function ComputePsi

%% Firing rate function
% ~~~~~~~~~~~~~~~


% Firing rate function. Equation (4) Freestone et al., 2011, NeuroImage

% x_t = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.
phi_fields = zeros(size(phi_basisFunctions));

for m = 1 : nx % to compute phi * x(t)
    
    phi_fields(:,:, m) = phi_basisFunctions(:,:, m) * Xt(m); % field basis function times state vector
    
end

Vt = sum(phi_fields, 3); % v, mean membrane potential field, at time t.

% firing rate function
firingRate_phi_Vt = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - Vt))); % firing rate sigmoidal function, field

%% integral over 2-D space
% ~~~~~~~~~~~~~~~


% compute x(T+1)
ingtegralProduct = zeros(nx, nTheta); % a nx times nTheta matrix of integrations

for pNX = 1 : nx % cycle through field basis functions
    
    for qNTheta = 1 : nTheta % cycle through basis functions of connectivity kernel
        
        product_psi_firingRate = squeeze(psi(qNTheta, pNX, :, :)) * firingRate_phi_Vt; % product of Psi and field after firing rate function
        
        ingtegralProduct(pNX, qNTheta) = sum(sum(product_psi_firingRate * stepSize^2, 2), 1); % integrate over space
        
    end
    
end

XtPlus1 = ks * Xt + ingtegralProduct * theta; % finally times theta (vector) and get x(t+1)

%% Compute field at time (T+1)
% ~~~~~~~~~~~~~~~


for m = 1 : nx % to compute phi * x(t)
    
    fields_tplus1(:,:, m) = phi_basisFunctions(:,:, m) * XtPlus1(m); % field basis functions times state vector
    
end

VtPlus1 = squeeze(sum(fields_tplus1, 3)); % field at time T+1

end
