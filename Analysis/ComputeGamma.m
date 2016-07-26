function gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu, sigma)
%%
% SpaceMin, SpaceMax, NPoints, cortical surface
% nx, number of basis functions
% sigma, covariance matrix of gaussian basis function
%% Generate data
% create a cortical surface
x = linspace(SpaceMin, SpaceMax, NPoints);
stepSize = x(2)-x(1);
[X, Y] = meshgrid(x, x);
%% define Gaussian basis functions
% define mu and sigma for gaussian basis functions
% numRow = sqrt(nx); % number of gaussians for each colomn
% numCol = nx / numRow; % number of columns
% 
% widthSpace = SpaceMax - SpaceMin;
% widthCentre = widthSpace / (numCol*2);
% 
% mu = zeros(nx, 2); % centres of each gaussian
% for m = 1 : numRow
%     for n = 1 : numCol
%         mu(n + numCol*(m-1), :) = [(SpaceMin - widthCentre + m*widthCentre*2) (SpaceMin - widthCentre + n*widthCentre*2)];
%     end
% end
% 
% % define gaussian basis functions
% gaussians = zeros(NPoints, NPoints, nx);
% for n = 1 : nx
%     gaussians(:,:, n) = Define2DGaussian_AnisotropicKernel(mu(n, 1), mu(n, 2), sigma, NPoints, SpaceMin, SpaceMax);
% end
%% plot
% figure, clf, shg; imagesc(squeeze(sum(gaussians, 3))), colorbar; title('Guassian basis functions in the field');
%% Compute gamma - analytic
gamma = zeros(nx, nx);

for m = 1 : nx
    for n = 1 : nx
        gamma(m, n) = InnerProductTwo2DGaussians(mu(m, :), mu(n, :), sigma, sigma);
    end
end

% plot
% figure, imagesc(gamma_analytic), colorbar, title('gamma matrix - analytic');
end