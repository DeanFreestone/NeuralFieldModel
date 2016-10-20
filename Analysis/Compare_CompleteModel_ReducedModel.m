%% Compare complete model and reduced model
%
% Miao Cao


clc
clear
close all

%% Parameters
% ~~~~~~~~~~~~~~~


% parameters to create a 2-D cortical surface
SpaceMin = -10; SpaceMax = 10; NPoints = 501;

nx = 16; % number of Gaussian basis function of field decomposition

sigma_phi = 0.8; % width of Gaussian basis function of field decomposition

% connectivity kernel
theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

nTheta = 3; % number of connectivity kernel basis functions

vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

%% Models
% ~~~~~~~~~~~~~~~

Xt = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.

% reduced model
[ReducedModel_VtPlus1, Vt] = ReducedModel_ComputeFieldVtPlus1(Xt, nx, sigma_phi, theta, nTheta, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);

% complete model
CompleteModel_VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, theta, nTheta, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);

% compare


%%
