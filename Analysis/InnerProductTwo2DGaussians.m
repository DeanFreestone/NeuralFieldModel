function [product] = InnerProductTwo2DGaussians(mu_1, mu_2, sigma_1, sigma_2)
%% Inner product of two 2-D Gaussians
% To implement Appendix D Equation (D.7). Freestone et al., 2011, NeuroImage
% Miao Cao

% ~~~~~~~~~~~~~~~
% input parameter list:
% mu_1 - centres of Gaussian 1
% mu _2- centres of Gaussian 2
% sigma_1 - sigma of Gaussian 1
% sigma_2 - sigma of Gaussian 2

% output parameter list:
% product - inner product of two gaussians

%%
mu = (mu_1 - mu_2)';
% var_phi = sigma_phi; var_psi = sigma_psi;
var_1 = sigma_1(1,1); var_2 = sigma_2(1,1);

CovMat = (var_1 + var_2);
exponential = exp(-(mu'/CovMat*mu));

coefficient = (pi*var_1*var_2)/(var_1 + var_2);
product = coefficient*exponential;

end
