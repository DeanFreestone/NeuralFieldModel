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

%% Models
% ~~~~~~~~~~~~~~~

x_t = randn(nx, 1, 'single'); % x(t), state vector at time t. Set as rand numbers for now.

% reduced model
[ReducedModel_VtPlus1, Vt] = ReducedModel_ComputeFieldVtPlus1(x_t, SpaceMin, SpaceMax, NPoints);

% complete model
CompleteModel_VtPlus1 = CompleteModel_ComputeFieldVtPlus1(Vt, SpaceMin, SpaceMax, NPoints);

% compare
