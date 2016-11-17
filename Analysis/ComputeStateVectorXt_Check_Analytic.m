%% Compute state vector x(t+1)
% To implement Equation (25), Freestone et al., 2011, NeuroImage
% Miao Cao


clc
clear
close all

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 201;
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);

%% Parameters used to compute psi
% ~~~~~~~~~~~~~~~


Ts = 0.0001; % time step
nx = 16; % number of Gaussian basis functions
theta = [1, 1.2, 1.5]'; % number of connectivity kernel basis functions
nTheta = 3;

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


psi = ComputePsi(X, Y, SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi); % compute psi with function ComputePsi

%% Firing rate function
% ~~~~~~~~~~~~~~~


% Firing rate function. Equation (4) Freestone et al., 2011, NeuroImage

x_t = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.
phi_fields = zeros(size(phi_basisFunctions));

for m = 1 : nx % to compute phi * x(t)
    phi_fields(:,:, m) = phi_basisFunctions(:,:, m) * x_t(m);
end
v = sum(phi_fields, 3); % v, mean membrane potential field, at time t.

% firing rate function
firingRate_phi_xt = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - v))); % firing rate sigmoidal function, field

%% integrate over 2-D space
% ~~~~~~~~~~~~~~~

x_tplus1 = []; % here is x(t+1)
ingtegralProduct = zeros(nx, nTheta);
for pNX = 1 : nx
    for qNTheta = 1 : nTheta
        product_psi_firingRate = squeeze(psi(qNTheta, pNX, :, :)) * firingRate_phi_xt;
        ingtegralProduct(pNX, qNTheta) = sum(sum(product_psi_firingRate * stepSize^2, 2), 1);
    end
end
x_tplus1 = ingtegralProduct * theta; % finally times theta (vector) and get x(t+1)
%% Numerical check. Implement equation (12)
% implementationf of Equation (12)


tau = 0.01; % synaptic time constant

v_tplus1 = []; % field at T+1

v_t = rand(NPoints, NPoints); % initialise a random field v at time point T

ks = 1- Ts*(1/tau); % time constant parameter

errorPart = zeros(NPoints, NPoints); % set error part to zero for now

%% integral part


% initialise integral part
integralPart = zeros(NPoints, NPoints);
% firing rate function
firingRate_v_t = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - v_t)));

% integral. convolution or integral
for m = 1 : NPoints
    for n = 1 : NPoints
        r = [X(m, n), Y(m, n)]; % location r vector
        
        % define connectivity kernel at location r
        % connectivity kernel, a sum of three gaussian kernels
        theta = [10, -8, 0.5]';
        sigma = [0.6 0.8 2];
        for p = 1 : 3
            gaussians(:,:, p) = Define2DGaussian_AnisotropicKernel(r(1), r(2), [sigma(p) 0; 0 sigma(p)], NPoints, SpaceMin, SpaceMax) * theta(p);
        end
        w = squeeze(sum(gaussians, 3));
        
        integralPart = integralPart + w.*v_t;
    end
end

%% v(t+1)


v_tplus1 = ks * v_t + Ts * integralPart + errorPart; % calculate v(t+1)

%% check

figure, imagesc(v_t), colorbar, title('v(t)');
figure, imagesc(integralPart), colorbar, title('integral part');
figure, imagesc(v_tplus1), colorbar, title('v(t+1)');

