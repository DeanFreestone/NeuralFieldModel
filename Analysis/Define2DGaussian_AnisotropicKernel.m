function pdf = Define2DGaussian_AnisotropicKernel(xmu, ymu, sigma, NPoints, SpaceMin, SpaceMax)
%% Function to define a 2-D anisotropic Gaussian pdf
% Parameter list:
% xmu - x location of the centre of Gaussian
% ymu - y location of the centre of Gaussian
% sigma - the variance-covariance matrix of Gaussian
% NPoints - spatial resolution, or number of discretisation (each row or column) of the surface
% SpaceMin - negative edge of the surface
% SpaceMax - positive edge of the surface



A = inv(sigma);                                 % the inverse covariance matrix

x = linspace(SpaceMin, SpaceMax, NPoints);    % xmu-PlotWidth*maxsd:0.1:xmu+PlotWidth*maxsd; % location of points at which x is calculated
[X, Y] = meshgrid(x, x);                                         % matrices used for plotting

% Compute value of Gaussian pdf at each point in the grid
pdf =  exp(-(A(1,1)*(X-xmu).^2 + 2*A(1,2)*(X-xmu).*(Y-ymu) + A(2,2)*(Y-ymu).^2));
% Note! This pdf has no coefficient.


% compute Gaussian from basic equation
% coefficient = 1/(2*pi*sd*sd*sqrt(1-rho^2));
% z = coefficient * exp(-1*(((X-xmu).^2)./var + ((Y-ymu).^2)./var - (2*rho*(X-xmu).*(Y-ymu))./var )/(2*(1-rho^2)));
