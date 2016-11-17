function VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, tau, Ts, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints)
%% Complete/Full model
% Compute neural field (post-synaptic membrane potential) at time (T+1)
% Based on the equation (12), Freestone et al. 2011 NeuroImage
% But we ignore error part in equation (12) for now.
% Miao Cao


% Parameter list:
% Vt - Neural field at time T
% theta - scale of basis functions of connectivity kernel
% nTheta - number of basis functions of connectivity kernel
% mu_psi - the centre of Gaussian basis function of connectivity kernel
% vector_Sigma_Psi - a vector of widiths of Gaussian basis functions of
% connvectivity kernel (spatial decomposition)
% SpaceMin - the negative edge of cortical surface/neural field
% SpaceMax - the posive edge of cortical surface/neural field
% NPoints - number of points in each row or column of cortical
% surface/neural field (spatial resolution)

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints); % a 2-D surface

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x); % x and y coordinates of the 2-D surface

%% model parameters
% ~~~~~~~~~~~~~~~



% tau = 0.01; % synaptic time constant

VtPlus1 = []; % field at T+1

% V at time t

%  v_t = ones(NPoints, NPoints);

% rng(0,'twister');
%
% lower = -6;
% upper = 6;
% v_t = (upper-lower).*rand(NPoints, NPoints) + lower;

% v_t = rand(NPoints, NPoints); % initialise a random field v at time point T



% Ts = 0.0001; % time step, temporal resolution

% nx = 16; % number of Gaussian basis functions

% theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

% vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

% nTheta = 3; % number of connectivity kernel basis functions

% ~~~~~~~~~~~~~~~
% parameters for firing rate function
slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

ks = 1- Ts*(1/tau); % time constant parameter

errorPart = zeros(NPoints, NPoints); % set error part to zero for now
%% integral part


% initialise integral part
integralPart = zeros(NPoints, NPoints);


% firing rate function
firingRate_v_t = 1 ./ ( 1 + exp(slope_sigmoidal*(v0 - Vt)));


% Compute connectivity kernel, decomposed into three basis functions
gaussians = zeros(NPoints, NPoints, nTheta);

for m = 1 : nTheta
    
    % covariance matrix of this basis function of connectivity kernel
    covMat_Psi = [vector_Sigma_Psi(m, 1) vector_Sigma_Psi(m, 2); vector_Sigma_Psi(m, 2) vector_Sigma_Psi(m, 1)];
    
    gaussians(:,:, m) = Define2DGaussian_AnisotropicKernel(mu_psi(m, 1), mu_psi(m, 2), covMat_Psi, NPoints, SpaceMin, SpaceMax) * theta(m); % define each Gaussian basis function
    
end
w = squeeze(sum(gaussians, 3)); % connectivity kernel

% convolution of connectivity kernel and neural field (after firing rate function)

%% extend the field to 3*3 fields/tiles, to solve errors induced by boundary condition
% ~~~~~~~~~~~~~~~


NPoints_ExtendedField = 3 * NPoints; % number of discretisation points in the extended field

extendedField = zeros(NPoints_ExtendedField, NPoints_ExtendedField); % extended neural field



% sysmetry boundary condition

% topLeft = rot90(firingRate_v_t, -2);
% top = flipud(firingRate_v_t);
% topRight = rot90(firingRate_v_t, -2);
% left = fliplr(firingRate_v_t);
% centre = firingRate_v_t;
% right = fliplr(firingRate_v_t);
% bottomLeft = rot90(firingRate_v_t, -2);
% bottom = flipud(firingRate_v_t);
% bottomRight = rot90(firingRate_v_t, -2);

% boundary condition - 9 identical copies into 9 tiles

topLeft = firingRate_v_t;
top = firingRate_v_t;
topRight = firingRate_v_t;
left = firingRate_v_t;
centre = firingRate_v_t;
right = firingRate_v_t;
bottomLeft = firingRate_v_t;
bottom = firingRate_v_t;
bottomRight = firingRate_v_t;


% top row
extendedField(1 : NPoints, 1 : NPoints) = topLeft; extendedField(1 : NPoints, NPoints+1 : 2 * NPoints) = top; extendedField(1 : NPoints, 2 * NPoints +1 : 3 * NPoints) = topRight;
% middle row
extendedField(NPoints+1 : 2 * NPoints, 1 : NPoints) = left; extendedField(NPoints+1 : 2 * NPoints, NPoints+1 : 2 * NPoints) = centre; extendedField(NPoints+1 : 2 * NPoints, 2 * NPoints +1 : 3 * NPoints) = right;
% bottom row
extendedField(2 * NPoints +1 : 3 * NPoints, 1 : NPoints) = bottomLeft; extendedField(2 * NPoints +1 : 3 * NPoints, NPoints+1 : 2 * NPoints) = bottom; extendedField(2 * NPoints +1 : 3 * NPoints, 2 * NPoints +1 : 3 * NPoints) = bottomRight;

%% convolution
% ~~~~~~~~~~~~~~~


integralPart = conv2(extendedField, w); % convolution of extended field and connectivity kernel

integralPart = integralPart(NPoints+1 : 2 * NPoints, NPoints+1 : 2 * NPoints); % get the field
%% V(t+1), neural field at time T+1

VtPlus1 = ks * Vt + Ts * integralPart + errorPart; % calculate v(t+1), neural field at time T+1

end
