% Dean Freestone. 
% test correlation analysis for estimation of kernel support

clc
clear
close all

UsePBC = true;
UseIso = true;

UseLinear = false;
UseLinearized = false;
UseNonlinear = true;



tic

% for plotting
% ~~~~~~~
FS_Label = 8;          % fontsize for the axis label
FS_Tick = 8;                % fontsize for the ticks
MS = 10;                     % marker size
LW = 1;
plotwidth = 4.4;        % cm
plotheight = 3.5;

% parameters

%% spatial parameters
% ~~~~~~~~~~~
Delta = 0.5;                          % space step for the spatial discretisation
Delta_squared = Delta^2;
SpaceMax = 10;                    % maximum space in mm
SpaceMin = -SpaceMax;         % minimum space in mm
NPoints = (SpaceMax-SpaceMin)/Delta+1;
NPoints_total = NPoints^2;
r = linspace(SpaceMin,SpaceMax,NPoints);      % define space

EstimationSpaceMax = 10;
EstimationSpaceMin = -10;

%% temporal parameters
% ~~~~~~~~~~~~~~~~
Ts = 1e-3;              % sampling period (s)
T = 20000;              % maximum time (ms)          % 2 seconds = 100 seconds

%% kernel parameters
% ~~~~~~~~~~~~~~
if UseIso
    theta(1) = 100;%80.0;           % local kernel amplitude
    theta(2) = -80;             % surround kernel amplitude
    theta(3) = 5;               % lateral kernel amplitude
    theta(4) = 0;%15;               % anisotropi amplitude

    cmax = 25.5;
    cmin = -10.0;
else
    theta(1) = 80.0;           % local kernel amplitude
    theta(2) = -80;             % surround kernel amplitude
    theta(3) = 5;               % lateral kernel amplitude
    theta(4) = 15;               % anisotropi amplitude
    
    cmax = 12.5;
    cmin = -12.5;
end
sigma_psi(1) = 1.8;     % local kernel width
sigma_psi(2) = 2.4;     % surround kernel width
sigma_psi(3) = 6;       % lateral kernel width
sigma_psi(4) = 22;       % anisotropic width

psi_0 = Define2DGaussian(0,0, sigma_psi(1)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_1 = Define2DGaussian(0,0, sigma_psi(2)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_2 = Define2DGaussian(0,0, sigma_psi(3)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_3 = Define2DGaussian(-3,0, sigma_psi(4)^2, 0,NPoints,SpaceMin,SpaceMax);

w = theta(1)*psi_0 + theta(2)*psi_1 + theta(3)*psi_2 + theta(4)*psi_3;       % the kernel

psi_0_large = Define2DGaussian(0,0, sigma_psi(1)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
psi_1_large = Define2DGaussian(0,0, sigma_psi(2)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
psi_2_large = Define2DGaussian(0,0, sigma_psi(3)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
psi_3_large = Define2DGaussian(-3,0, sigma_psi(4)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);

w_large = theta(1)*psi_0_large + theta(2)*psi_1_large + theta(3)*psi_2_large + theta(4)*psi_3_large;       % the large kernel

filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/Kernel.pdf';
figure('units','centimeters','position',[0 0 plotwidth plotheight],'filename',filename,...
    'papersize',[plotheight, plotwidth],'paperorientation','landscape','renderer','painters')  

imagesc(r,r,w,[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar

%% sensor parameters
% ~~~~~~~~~~~~~~~
sensor_index = 1:3:41;
NSensors_xy = length(sensor_index);
NSensors = NSensors_xy^2;
sigma_m = 0.9;              % sensor width
m = Define2DGaussian(0, 0, sigma_m^2, 0, NPoints, SpaceMin, SpaceMax);
m_large = Define2DGaussian(0, 0, sigma_m^2, 0, 2*NPoints-1,2*SpaceMin,2*SpaceMax);
delta_y = r(sensor_index(2)) - r(sensor_index(1));

%% observation noise characteristics
% ~~~~~~~~~~~~~~~~~~~~
sigma_varepsilon = 0.1;                                                                     
Sigma_varepsilon = sigma_varepsilon^2*eye(NSensors);        % observation covariance matrix

% create the observation noise
varepsilon = mvnrnd(zeros(1,NSensors),Sigma_varepsilon,T);
    
%% sigmoid parameters
% ~~~~~~~~~~~~~~~~
% f_max = 1;                  % maximum firing rate
varsigma = 0.56;         % sigmoid slope
v_0 = 1.8;                    % firing threshold

%% synaptic kernel parameter
% ~~~~~~~~~~~~~~~~~~~~
tau = 0.01;                   % synaptic time constant
zeta = 1/tau;                 % inverse synaptic time constant
xi = 1-Ts*zeta;           % coefficient for the discrete time model

%% disturbance paramters
% ~~~~~~~~~~~~~
sigma_gamma = 1.3;              % parameter for covariance of disturbance
gamma_weight = 0.1;            % variance of disturbance

if UsePBC
    SphericalBoundary               % run the script that generates the covariance matrix
else
    mm=1;
    Sigma_gamma = zeros(NPoints^2,NPoints^2);   % create disturbance covariance matrix
    for n=1:NPoints
        for nn=1:NPoints
            temp = gamma_weight*Define2DGaussian(r(n),r(nn), sigma_gamma^2, 0,NPoints,SpaceMin,SpaceMax);
            Sigma_gamma(:,mm) = temp(:);
            mm=mm+1;
        end
    end
end

gamma = gamma_weight*Define2DGaussian(0,0, sigma_gamma^2, 0,NPoints,SpaceMin,SpaceMax);
e = mvnrnd(zeros(1,NPoints^2),Sigma_gamma,T);

%% set initial condition
% ~~~~~~~~~~~~~~~~
v = zeros(T,NPoints,NPoints);
y = zeros(T,NSensors_xy,NSensors_xy);

max_field_init = gamma_weight;
min_field_init = -gamma_weight;
InitialCondition = min_field_init + (max_field_init - min_field_init)*rand(NPoints,NPoints);
v(1,:,:) = InitialCondition;

if UsePBC
    v_temp = padarray(InitialCondition,size(InitialCondition),'circular');
    y_full = conv2(v_temp,m_large,'valid')*Delta_squared;      % the full field filtered by the larger sensor kernel
    y_full = y_full(2:end-1,2:end-1);
else
    v_temp = InitialCondition;
    y_full = conv2(m,v_temp,'same')*Delta_squared;      % the full field filtered by the smaller sensor kernel
end

varepsilon_t = reshape(varepsilon(1,:,:), NSensors_xy, NSensors_xy);
y(1,:,:) = y_full(sensor_index,sensor_index) + varepsilon_t;

%% Generate data
for t=1:T-1
    
    y_t = squeeze(y(t,:,:));    
    v_t = squeeze(v(t,:,:));

    % calculate the firing rate
    if UseLinear
        f = varsigma*v_t;
    elseif UseLinearized
        f = ( 2 + varsigma*(v_t - v_0) ) / 4;
    elseif UseNonlinear
        f = 1./(1+exp( varsigma*(v_0-v_t) ));           % calc firing rate using sigmoid
    end
    
    if UsePBC
        f = padarray(f,size(f),'circular');
        g = conv2(f,w_large,'valid')*Delta_squared;   
        g = g(2:end-1,2:end-1);
    else
        g = conv2(w,f,'same')*Delta_squared;   
    end

    % conv firing rate with spatial kernel
    e_t = reshape(e(t,:,:),NPoints,NPoints);
    v_t_plus1 = Ts*g + xi*v_t + e_t;  % update field 
    v(t+1,:,:) = v_t_plus1;

    if UsePBC
        v_t_plus1 = padarray(v_t_plus1,size(v_t_plus1),'circular');
        y_full = conv2(v_t_plus1,m_large,'valid') * Delta_squared;                  % field filtered by observation kernel at t+1
        y_full = y_full(2:end-1,2:end-1);
    else
        y_full = conv2(m,v_t_plus1,'same') * Delta_squared; 
    end

    % filter field with sensor kernel and get observations
    varepsilon_tplus1 = reshape(varepsilon(t+1,:,:),NSensors_xy,NSensors_xy);
    y_tplus1 = y_full(sensor_index,sensor_index) + varepsilon_tplus1;                 % discretize sensor spacing
    y(t+1,:,:) = y_tplus1;
end

toc
