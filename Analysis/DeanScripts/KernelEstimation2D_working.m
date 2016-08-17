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

CheckDerivation = false;

tic

% for plotting
% ~~~~~~~
FS_Label = 10;          % fontsize for the axis label
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

%% spatial kernel parameters
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
sigma_psi(4) = 2;       % anisotropic width

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
% Sigma_varepsilon = sigma_varepsilon^2*eye(NSensors);        % observation covariance matrix
Sigma_varepsilon = 0.1*eye(NSensors);        % observation covariance matrix

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
    mm = 1;
    Sigma_gamma = zeros(NPoints^2,NPoints^2);   % create disturbance covariance matrix
    for n=1:NPoints
        for nn=1:NPoints
            temp = gamma_weight*Define2DGaussian(r(n),r(nn), sigma_gamma^2, 0,NPoints,SpaceMin,SpaceMax);
            Sigma_gamma(:,mm) = temp(:);
            mm = mm+1;
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

%% Initialize the correlation matrices to zeros for speed
R_y_tplus1y = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);
R_yy = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);
R_yy_tplus1 = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);

%% 
if UseLinearized || UseNonlinear
    if UsePBC
        const_field = (2-varsigma*v_0)*ones(3*NPoints,3*NPoints)/4;
        w_conv_const = Ts*conv2(const_field,w_large,'valid')*Delta_squared; 
        w_conv_const = w_conv_const(2:end-1,2:end-1);       % trim result to the field space
    else
        const_field = (2-varsigma*v_0)*ones(NPoints,NPoints)/4;
        w_conv_const = Ts*conv2(w, const_field,'same')*Delta_squared; 
    end

    if UsePBC
        w_conv_const = padarray(w_conv_const,size(w_conv_const),'circular');                                 % set up the periodic boundary
        term1_temp1 = conv2(w_conv_const,m_large,'valid')*Delta_squared;
        term1_temp1 = term1_temp1(2:end-1,2:end-1);             % trim result
    else
        term1_temp1 = conv2(m, w_conv_const,'same')*Delta_squared;
    end

    c_1 = term1_temp1(sensor_index,sensor_index);
%     c_1_temp =  Ts * sum(sum(w_large)) * Delta_squared * sum(sum(m_large)) * Delta_squared * (2-varsigma*v_0)/4;
    c_2 = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);
end

% figure
%% Generate data
for t=1:T-1
    
    y_t = squeeze(y(t,:,:));    
    v_t = squeeze(v(t,:,:));

    R_yy(t,:,:) = xcorr2(y_t);

    % calculate the firing rate
    if UseLinear
        f = varsigma*v_t;
    elseif UseLinearized
        c_2(t,:,:) = xcorr2(c_1,y_t);
%         c_2(t) = c_1*sum(sum(y_t));
        f = ( 2 + varsigma*(v_t - v_0) ) / 4;
    elseif UseNonlinear
        c_2(t,:,:) = xcorr2(c_1,y_t);
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
%     imagesc(v_t_plus1)
%     colorbar
%     drawnow
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

    R_y_tplus1y(t,:,:) = xcorr2(y_tplus1,y_t);
    R_yy_tplus1(t,:,:) = xcorr2(y_t,y_tplus1);
    
    %%
    if UseLinearized && CheckDerivation
    % this bit is used to check if the derivation is correct
%     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        % calc the RHS of eq 11
        if UsePBC
            g = padarray(g,size(g),'circular');                                 % set up the periodic boundary
            m_conv_g_full_temp = Ts*conv2(g,m_large,'valid') * Delta_squared;         
            m_conv_g_full_temp = m_conv_g_full_temp(2:end-1,2:end-1);
        else
            m_conv_g_full_temp = Ts*conv2(m,g,'same') * Delta_squared;         
        end
        m_conv_g_full = m_conv_g_full_temp(sensor_index,sensor_index);
        RHS = xcorr2(m_conv_g_full,y_t);

        %calc the LHS of eq 11
        varepsilon_t = reshape(varepsilon(t,:,:),NSensors_xy,NSensors_xy);
        R_varepsilon_y = xcorr2(varepsilon_tplus1,y_t);
        if UsePBC
            e_t = padarray(e_t,size(e_t),'circular');                                 % set up the periodic boundary
            m_conv_disturbance_full = conv2(e_t,m_large,'valid') * Delta_squared;
            m_conv_disturbance_full = m_conv_disturbance_full(2:end-1,2:end-1);
        else
            m_conv_disturbance_full = conv2(m,e_t,'same') * Delta_squared;
        end
        m_conv_disturbance = m_conv_disturbance_full(sensor_index,sensor_index);
        R_m_conv_disturb_y = xcorr2(m_conv_disturbance,y_t);
        if UsePBC               % do not comment this line. It is used later.
            v_t = padarray(v_t,size(v_t),'circular');                                 % set up the periodic boundary
            y_no_noise_full = conv2(v_t,m_large,'valid')*Delta_squared;
            y_no_noise_full = y_no_noise_full(2:end-1,2:end-1);
        else
            y_no_noise_full = conv2(m,v_t,'same')*Delta_squared;
        end
        y_no_noise = y_no_noise_full(sensor_index,sensor_index);
        LHS = squeeze(R_y_tplus1y(t,:,:)) - xi*xcorr2(y_no_noise,y_t) - R_m_conv_disturb_y - R_varepsilon_y ;

    % ~~~~~~~~~~~~~~~~~~~~~~
    % full version of term 1
    % ~~~~~~~~~~~~~~~~~~~~~~
        if UsePBC
            const_field = (2-varsigma*v_0)*ones(3*NPoints,3*NPoints)/4;
            w_conv_const = Ts*conv2(const_field,w_large,'valid' )*Delta_squared;    
            w_conv_const = w_conv_const(2:end-1,2:end-1);
        else
            const_field = (2-varsigma*v_0)*ones(NPoints,NPoints)/4;
            w_conv_const = Ts*conv2(w, const_field,'same')*Delta_squared;    
        end
        
        if UsePBC
            w_conv_const = padarray(w_conv_const,size(w_conv_const),'circular');                                 % set up the periodic boundary
            term1_temp1 = conv2(w_conv_const,m_large,'valid' )*Delta_squared;
            term1_temp1 = term1_temp1(2:end-1,2:end-1);
        else
            term1_temp1 = conv2(m, w_conv_const,'same')*Delta_squared;
        end
        term1_temp2 = term1_temp1(sensor_index,sensor_index);
        term1 = xcorr2(term1_temp2,y_t);

    % simplified version of term 1
    %     term1_temp1 = sum(sum(m))*Delta_squared * Ts*sum(sum(w))*Delta_squared*(2-varsigma*v_0)/4;
    %     term1_temp2 = term1_temp1*ones(size(y_t));
    %     term1 = xcorr2(term1_temp2,y_t);
    % ~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~~~~~~~~~~~
    % full version of term2
    % ~~~~~~~~~~~~~~~~~~~~~~
        % v is already set up as with PBC as req.
        if UsePBC
            w_conv_v = (varsigma/4)*Ts*conv2(v_t,w_large,'valid')*Delta_squared;        
            w_conv_v = w_conv_v(2:end-1,2:end-1);
            w_conv_v = padarray(w_conv_v,size(w_conv_v),'circular');                                 % set up the periodic boundary
            term2_temp1 = conv2(w_conv_v,m_large,'valid')*Delta_squared;
            term2_temp1 = term2_temp1(2:end-1,2:end-1);
        else
            w_conv_v = (varsigma/4)*Ts*conv2(w,v_t,'same')*Delta_squared;        
            term2_temp1 = conv2(m,w_conv_v,'same')*Delta_squared;
        end
        term2_temp2 = term2_temp1(sensor_index,sensor_index);
        term2 = xcorr2(term2_temp2,y_t);

    %  simplified version of term2
    % ~~~~~~~~~~~~~~~~~~~~~~
        R_y_no_noise_y = xcorr2(y_no_noise,y_t);
    % %     trim to the central part
    %     midpoint = ceil(size(R_y_no_noise_y,1)/2);
    %     startpoint = midpoint - floor(size(y_t,1)/2);
    %     endpoint = midpoint + floor(size(y_t,1)/2);
    %     R_y_no_noise_y = R_y_no_noise_y(startpoint:endpoint,startpoint:endpoint);

        if UsePBC
    %         R_y_no_noise_y = padarray(R_y_no_noise_y,size(R_y_no_noise_y),'circular');                                 % set up the periodic boundary
        end
        w_sampled = w(sensor_index,sensor_index);
        reversed_w = flipud(fliplr(w_sampled));
        term2_simp = (Ts*varsigma/4) * conv2(R_y_no_noise_y,reversed_w,'same')*delta_y^2;
    %     term2_simp = (Ts*varsigma/4) * conv2(R_y_no_noise_y,flipud(fliplr(w)))*delta_y^2;

    %     midpoint = ceil(size(term2_simp,1)/2);
    %     startpoint = midpoint - floor(size(R_y_no_noise_y/2)/2);
    %     endpoint = midpoint + floor(size(R_y_no_noise_y/2)/2);
    %     term2_simp = term2_simp(startpoint:endpoint,startpoint:endpoint);
    % ~~~~~~~~~~~~~~~~~~~~~~

        subplot(231),imagesc(RHS),colorbar,axis square,title('RHS eq. 11')
        subplot(232),imagesc(term1+term2),colorbar,axis square,title('term1 + term2')
        subplot(233),imagesc(RHS - (term1+term2)),colorbar,axis square,title('error')
        subplot(234),imagesc(LHS),colorbar,axis square,title('LHS eq. 11')
        subplot(235),imagesc(LHS-RHS),colorbar,axis square,title('LHS - RHS eq. 11')
%         subplot(234),imagesc(term2),colorbar,axis square,title('term2')
%         subplot(235),imagesc(term2_simp),colorbar,axis square,title('term2\_simp')
        drawnow
        pause(0.1) 
    end
%     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end

%% Initialise for speed
R_xcorr = squeeze(mean(R_y_tplus1y(500:end-1,:,:),1));
R_acorr = squeeze(mean(R_yy(500:end-1,:,:),1));

R_xcorr_reverse = squeeze(mean(R_yy_tplus1(500:end-1,:,:),1));

temp = zeros(size(R_acorr));
temp(ceil(size(temp,1)/2),ceil(size(temp,2)/2)) = sigma_varepsilon^2*NSensors;

F_R_acorr = fft2(R_acorr - temp);       % substracting the observation noise covairance

if UseLinear
    
    F_R_xcorr = fft2(R_xcorr);
    
elseif UseLinearized || UseNonlinear
    
    mean_c_2 = squeeze(mean(c_2(500:end-1,:,:)));
    F_R_xcorr = fft2(R_xcorr - mean_c_2);
%         F_R_xcorr = fft2(R_xcorr );

end

W_est = F_R_xcorr ./ F_R_acorr - xi; 
w_est = fftshift(ifft2((W_est))) / (delta_y^2 * varsigma*Ts);

if UseLinearized || UseNonlinear
    w_est = w_est*4;
end

%% plot results
if UsePBC
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullNonlinearPBC.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullLinearizedPBC.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullLinearPBC.pdf';
    end
else
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullNonlinear.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullLinearized.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/FFTKernelEstimateFullLinear.pdf';
    end
end
figure('units','centimeters','position',[0 0 plotwidth plotheight],'filename',filename,...
    'papersize',[plotheight, plotwidth],'paperorientation','landscape','renderer','painters')
nu = linspace(0,2*delta_y*2,size(W_est,1));
imagesc(nu,nu,abs(W_est))
xlabel('Spatial Freq','fontsize',FS_Label)
ylabel('Spatial Freq','fontsize',FS_Label)
set(gca,'xtick',[0 1 2],'ytick',[0 1 2],'fontsize',FS_Tick)
axis square
axis xy
colorbar
colormap hot
drawnow

if UsePBC
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullNonlinearPBC.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullLinearizedPBC.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullLinearPBC.pdf';
    end
else
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullNonlinear.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullLinearized.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelEstimateFullLinear.pdf';
    end
end
figure('units','centimeters','position',[0 0 plotwidth plotheight],'filename',filename,...
    'papersize',[plotheight, plotwidth],'paperorientation','landscape','renderer','painters')   
r_est = linspace(-20,20,size(w_est,1));
imagesc(r_est,r_est,(w_est),[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar
% colormap hot

toc

if UsePBC
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionNonlinearPBC.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionLinearizedPBC.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionLinearPBC.pdf';
    end
else
    if UseNonlinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionNonlinear.pdf';
    elseif UseLinearized
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionLinearized.pdf';
    elseif UseLinear
        filename = '/Users/dean/Projects/BrainIDE/ltx/EMBCCorrelationAnalysisPaper/ltx/figures/KernelCrossSectionLinear.pdf';
    end
end
figure('units','centimeters','position',[0 0 plotwidth plotheight],'filename',filename,...
    'papersize',[plotheight, plotwidth],'paperorientation','landscape','renderer','painters')

plot(r_est,w_est(14,:),'r'), hold
plot(r,w(21,:),'k')

xlim([-10,10])
ylim([cmin,cmax])
leg = legend('$\hat{w}(\mathbf{r})$','$w(\mathbf{r})$');
set(leg,'interpreter','latex','box','off')
xlabel('Space','fontsize',FS_Label)
ylabel('Connection Strength','fontsize',FS_Label)
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10.0],'fontsize',FS_Tick)

%% estimate disturbance support

MMGamma_est = fft2(R_acorr/NSensors) - sigma_varepsilon^2 ...
    - xi*fft2(R_xcorr_reverse/NSensors)...
    - fft2(w_est).*fft2(mean_c_2/NSensors+varsigma*R_xcorr_reverse/NSensors)*Ts/4;

mmgamma_est = ifft2(MMGamma_est);

gamma_conv_m_conv_m = conv2(gamma,conv2(m_large,m_large,'same')*Delta_squared,'same')*Delta_squared;

figure
subplot(121)
imagesc(gamma_conv_m_conv_m)
axis square
colorbar
subplot(122)
imagesc(r_est,r_est,abs(mmgamma_est))
xlim([-10,10])
ylim([-10,10])
axis square
colorbar

% fft2(mean_c_2) - F_R_xcorr .* fft2(R_xcorr_reverse) ./ F_R_acorr;
% mmgamma_est = fftshift(ifft2(MMGamma_est)) / delta_y^2;
% 
% mmgamma = conv2(m,conv2(m,gamma,'same'),'same')*Delta_squared^2;
% figure,subplot(121),imagesc(r,r,mmgamma),axis square,title('true'),colorbar
% subplot(122),imagesc(r_est,r_est,real(mmgamma_est),[0,0.3]),axis square,title('estimate'),colorbar
% xlim([-10 10])
% ylim([-10 10])

%%
tempxxx = zeros(T,41,41);
for n=1:T
    tempxxx(n,:,:) = cov(squeeze(v(n,:,:)));
    
end
figure,imagesc(squeeze(mean(tempxxx,1)))
figure,imagesc(Sigma_gamma)