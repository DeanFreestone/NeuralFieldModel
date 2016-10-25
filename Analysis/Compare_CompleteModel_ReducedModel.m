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
SpaceMin = -10; SpaceMax = 10; NPoints = 201;

% field basis function parameters
nx = 225; % number of Gaussian basis function of field decomposition

sigma_phi = 1.0; % width of Gaussian basis function of field decomposition

% connectivity kernel parameters
theta = [10, -8, 0.5]'; % scale Gaussian basis functions of connectivity kernel

nTheta = 3; % number of connectivity kernel basis functions

mu_psi = [0 0; 0 0; 0 0]; % centres of basis functions of connectivity kernel

vector_Sigma_Psi = [0.6 0.8 2]; % width of Gaussian basis functions of connectivity kernel

%% Compare models with Mexican-hat connectvity kernels
% ~~~~~~~~~~~~~~~


Xt = randn(nx, 1, 'single')*10; % x(t), state vector at time t. Set as rand numbers for now.

% reduced model
[ReducedModel_VtPlus1, Vt] = ReducedModel_ComputeFieldVtPlus1(Xt, nx, sigma_phi, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);


% complete model
CompleteModel_VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, theta, nTheta,mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);

% residual of two models
residual = CompleteModel_VtPlus1 - ReducedModel_VtPlus1;

% calculate squared error
sqrError = sum(sum(residual.^2, 2), 1)

% compare
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
imagesc(Vt), colorbar; title('V(t)');
subplot(2,2,2);
imagesc(ReducedModel_VtPlus1), colorbar; title('Reduced V(t+1)');
subplot(2,2,3);
imagesc(CompleteModel_VtPlus1), colorbar; title('Full V(t+1)');
subplot(2,2,4);
imagesc(CompleteModel_VtPlus1 - ReducedModel_VtPlus1), colorbar; title('Residual');
suptitle('Mexican Hat kernel');

%% Compare models with Gabor-kernel connectvity kernels
% ~~~~~~~~~~~~~~~


% connectivity kernel - Gabor
theta = [5 -5]'; % scale Gaussian basis functions of connectivity kernel

nTheta = 2; % number of connectivity kernel basis functions

mu_psi = [-0.5 0; 0 0.5]; % centres of basis functions of connectivity kernel

vector_Sigma_Psi = [0.8 0.8]; % width of Gaussian basis functions of connectivity kernel



Xt = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.

% reduced model
[ReducedModel_VtPlus1, Vt] = ReducedModel_ComputeFieldVtPlus1(Xt, nx, sigma_phi, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);


% complete model
CompleteModel_VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, theta, nTheta, mu_psi, vector_Sigma_Psi, SpaceMin, SpaceMax, NPoints);


% compare
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
imagesc(Vt), colorbar; title('V(t)');
subplot(2,2,2);
imagesc(ReducedModel_VtPlus1), colorbar; title('Reduced V(t+1)');
subplot(2,2,3);
imagesc(CompleteModel_VtPlus1), colorbar; title('Full V(t+1)');
subplot(2,2,4);
imagesc(CompleteModel_VtPlus1 - ReducedModel_VtPlus1), colorbar; title('Residual');
suptitle('Gabor');
