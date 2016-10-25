%% Compare complete model and reduced model
%
% Miao Cao


clc
clear
close all
%% add function path
% ~~~~~~~~~~~~~~~


addpath(genpath('./Functions/'));

%% Parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 101;

% field basis function parameters
nx = 121; % number of Gaussian basis function of field decomposition

sigma_phi = 1.5; % width of Gaussian basis function of field decomposition

% connectivity kernel parameters
theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

nTheta = 3; % number of connectivity kernel basis functions

mu_psi = [0 0; 0 0; 0 0]; % centres of basis functions of connectivity kernel

vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

%% Compare models with Mexican-hat connectvity kernels
% ~~~~~~~~~~~~~~~


Xt = randn(nx, 1, 'single')*10; % x(t), state vector at time t. Set as rand numbers for now.
figure('units','normalized','outerposition',[0 0 1 1]);

for t = 1 : 100
    % reduced model
    [ReducedModel_VtPlus1, Vt, XtPlus1] = ReducedModel_ComputeFieldVtPlus1(Xt, nx, sigma_phi, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);
    shg, imagesc(ReducedModel_VtPlus1), colorbar; title(['t:' num2str(t)]);
    Xt = XtPlus1;
    pause(0.1)
end

% compare
