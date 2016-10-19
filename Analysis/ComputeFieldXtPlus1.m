%% Numerical check. Implement equation (12)
% implementationf of Equation (12)

clc
clear
close

%% Spatial parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 201;

x = linspace(SpaceMin, SpaceMax, NPoints);

stepSize = x(2)-x(1);

[X, Y] = meshgrid(x, x);

%%

tau = 0.01; % synaptic time constant

v_tplus1 = []; % field at T+1


% V at time t

%  v_t = ones(NPoints, NPoints);
 
rng(0,'twister');

lower = -6;
upper = 6;
v_t = (upper-lower).*rand(NPoints, NPoints) + lower;

% v_t = rand(NPoints, NPoints); % initialise a random field v at time point T




ks = 1- Ts*(1/tau); % time constant parameter

errorPart = zeros(NPoints, NPoints); % set error part to zero for now


Ts = 0.0001; % time step

nx = 16; % number of Gaussian basis functions

theta = [1, 1.2, 1.5]'; % number of connectivity kernel basis functions

nTheta = 3;

% ~~~~~~~~~~~~~~~
% parameters for firing rate function
slope_sigmoidal = 0.56; % slope of sigmoidal activation function

v0 = 1.8; % Firing threshold

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
