function VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints)
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



tau = 0.01; % synaptic time constant

VtPlus1 = []; % field at T+1

% V at time t

%  v_t = ones(NPoints, NPoints);

% rng(0,'twister');
%
% lower = -6;
% upper = 6;
% v_t = (upper-lower).*rand(NPoints, NPoints) + lower;

% v_t = rand(NPoints, NPoints); % initialise a random field v at time point T



Ts = 0.0001; % time step, temporal resolution

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
for p = 1 : nTheta
    gaussians(:,:, p) = Define2DGaussian_AnisotropicKernel(mu_psi(p, 1), mu_psi(p, 2), [vector_Sigma_Psi(p) 0; 0 vector_Sigma_Psi(p)], NPoints, SpaceMin, SpaceMax) * theta(p); % define each Gaussian basis function
end
w = squeeze(sum(gaussians, 3)); % connectivity kernel

% convolution
integralPart = conv2(w, firingRate_v_t, 'same');

% % integral. convolution or integral
% for m = 1 : NPoints
%     for n = 1 : NPoints
%         r = [X(m, n), Y(m, n)]; % location r vector
%
%         % connectivity kernel, a sum of Gaussian basis functions
%         for p = 1 : nTheta
%             gaussians(:,:, p) = Define2DGaussian_AnisotropicKernel(r(1), r(2), [vector_Sigma_Psi(p) 0; 0 vector_Sigma_Psi(p)], NPoints, SpaceMin, SpaceMax) * theta(p);
%         end
%         w = squeeze(sum(gaussians, 3)); % connectivity kernel
%
%         % define connectivity kernel at location r
%         integralPart = integralPart + w.*firingRate_v_t;
%     end
% end

%% V(t+1), neural field at time T+1

VtPlus1 = ks * Vt + Ts * integralPart + errorPart; % calculate v(t+1), neural field at time T+1

end
