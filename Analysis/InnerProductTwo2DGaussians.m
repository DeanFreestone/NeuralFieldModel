function [product] = InnerProductTwo2DGaussians(mu_1, mu_2, sigma_1, sigma_2)
%% Inner product of two 2-D Gaussian basis functions
% To implement Appendix D Equation (D.7). Freestone et al., 2011, NeuroImage
% Miao Cao


% ~~~~~~~~~~~~~~~
% input parameter list:
% mu_1 - centre of Gaussian basis function 1
% mu_2- centre of Gaussian basis function 2
% sigma_1 - sigma/variance-covariance matrix of Gaussian basis function 1
% sigma_2 - sigma/variance-covariance matrix of Gaussian basis function 2

% output parameter list:
% product - inner product of two gaussians

%% implement the equation
% ~~~~~~~~~~~~~~~


mu = (mu_1 - mu_2)'; % centre of the resultant Gaussian

covMat_1 = sigma_1; covMat_2 = sigma_2; % variance-covariance matrix

CovMat = (covMat_1 + covMat_2);

exponential = exp(-(mu'/CovMat*mu)); % exponential part

% coefficient = (2*pi*sqrt(det(covMat_1))*sqrt(det(covMat_2))) / sqrt(det(CovMat)); % coefficient
coefficient = (pi*det(covMat_1)*det(covMat_2)) / sqrt(det(CovMat)); % coefficient

product = coefficient*exponential; % product of coefficient and exponential

end
