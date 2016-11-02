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
gaussians = zeros(NPoints, NPoints, nTheta);

for m = 1 : nTheta
    
    % covariance matrix of this basis function of connectivity kernel
    covMat_Psi = [vector_Sigma_Psi(m, 1) vector_Sigma_Psi(m, 2); vector_Sigma_Psi(m, 2) vector_Sigma_Psi(m, 1)];
    
    gaussians(:,:, m) = Define2DGaussian_AnisotropicKernel(mu_psi(m, 1), mu_psi(m, 2), covMat_Psi, NPoints, SpaceMin, SpaceMax) * theta(m); % define each Gaussian basis function
    
end
w = squeeze(sum(gaussians, 3)); % connectivity kernel

% convolution of connectivity kernel and neural field (after firing rate function)

%%
integralPart = conv2(w, firingRate_v_t, 'same');



% spatial parameters
% ~~~~~~~~~~~
Delta = 0.5;                          % space step for the spatial discretisation
Delta_squared = Delta^2;
SpaceMaxPeriodicField = 30;                    % maximum space in mm
SpaceMinPeriodicField = -30;         % minimum space in mm
NPointsInPeriodicField = (SpaceMaxPeriodicField-SpaceMinPeriodicField)/Delta+1;
NPointsInField = (NPointsInPeriodicField-1)/3 + 1;
r = linspace(SpaceMinPeriodicField/3,SpaceMaxPeriodicField/3,NPointsInField);      % define space

% temporal parameters
% ~~~~~~~~~~~~~
% Ts = 1e-3;          % sampling period (s)
% T = 500;            % maximum time (ms)

% disturbance paramters
% ~~~~~~~~~~~~~
% sigma_gamma = 1.3;          % parameter for covariance of disturbance
% gamma_weight = 0.1;            % variance of disturbance

% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
FirstThirdEnd = (NPointsInPeriodicField-1)/3;
SecondThirdEnd = 2*FirstThirdEnd+1;

mm=1;
Sigma_gamma = zeros(NPointsInField^2,NPointsInField^2);   % create disturbance covariance matrix
% figure
template = zeros(NPointsInField,NPointsInField);
for n=1:NPointsInField
    for nn=1:NPointsInField
        temp = gamma_weight*Define2DGaussian(r(n),r(nn), sigma_gamma^2, 0, NPointsInPeriodicField, SpaceMinPeriodicField,SpaceMaxPeriodicField);
        
        topleft =[zeros(1,FirstThirdEnd+1) ; [zeros(FirstThirdEnd,1) temp(1:FirstThirdEnd,1:FirstThirdEnd)]];         
        top = [zeros(1,FirstThirdEnd+1) ;temp(1:FirstThirdEnd,FirstThirdEnd+1:SecondThirdEnd)];
        topright = [zeros(1,FirstThirdEnd+1) ; [temp(1:FirstThirdEnd,SecondThirdEnd+1:end) zeros(FirstThirdEnd,1)]]; 
        
        left = [zeros(FirstThirdEnd+1,1) temp(FirstThirdEnd+1:SecondThirdEnd,1:FirstThirdEnd)];
        middle = temp(FirstThirdEnd+1:SecondThirdEnd,FirstThirdEnd+1:SecondThirdEnd);
        right = [temp(FirstThirdEnd+1:SecondThirdEnd,SecondThirdEnd+1:end) zeros(FirstThirdEnd+1,1)];
        
        bottomleft = [zeros(FirstThirdEnd+1,1) [temp(SecondThirdEnd+1:end,1:FirstThirdEnd) ; zeros(1,FirstThirdEnd)]];
        bottom = [temp(SecondThirdEnd+1:end,FirstThirdEnd+1:SecondThirdEnd) ; zeros(1,FirstThirdEnd+1)];
        bottomright = [[temp(SecondThirdEnd+1:end,SecondThirdEnd+1:end) zeros(FirstThirdEnd,1)] ; zeros(1,FirstThirdEnd+1)];
        
        temp2 = middle + topleft + top + topright + left + right + bottom + bottomleft + bottomright;
        
        Sigma_gamma(:,mm) = temp2(:);
        mm=mm+1;
        
%         clim = [-200 -1];         % for log
%         clim = [0 0.08];
% 
%         subplot(3,3,1),imagesc(topleft,clim)
%         subplot(3,3,2),imagesc(top,clim)
%         subplot(3,3,3),imagesc(topright,clim)
%         
%         subplot(3,3,4),imagesc(left,clim)
%         subplot(3,3,5),imagesc(middle,clim)
%         subplot(3,3,6),imagesc(right,clim)
%         
%         subplot(3,3,7),imagesc(bottomleft,clim)
%         subplot(3,3,8),imagesc(bottom,clim)
%         subplot(3,3,9),imagesc(bottomright,clim)
% % 
%  %       imagesc(log10(temp2))
%         drawnow
%         clim = [0 0.08];
% % 
%         subplot(3,3,1),imagesc(topleft,clim)
%         subplot(3,3,2),imagesc(top,clim)
%         subplot(3,3,3),imagesc(topright,clim)
%         
%         subplot(3,3,4),imagesc(left,clim)
%         subplot(3,3,5),imagesc(middle,clim)
%         subplot(3,3,6),imagesc(right,clim)
%         
%         subplot(3,3,7),imagesc(bottomleft,clim)
%         subplot(3,3,8),imagesc(bottom,clim)
%         subplot(3,3,9),imagesc(bottomright,clim)
% % 
% %         imagesc(log10(temp2))
%         drawnow
     end
 end


%% V(t+1), neural field at time T+1

VtPlus1 = ks * Vt + Ts * integralPart + errorPart; % calculate v(t+1), neural field at time T+1

end
