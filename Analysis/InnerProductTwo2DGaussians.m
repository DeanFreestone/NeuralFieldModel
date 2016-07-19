function [product] = InnerProductTwo2DGaussians(mu_phi, mu_psi, sigma_phi, sigma_psi)

mu = (mu_phi - mu_psi)';
var_phi = sigma_phi(1,1); var_psi = sigma_psi(1,1);

CovMat = (var_phi + var_psi);
exponential = exp(-(mu'/CovMat*mu));

coefficient = (pi*var_phi*var_psi)/(var_phi + var_psi);
product = coefficient*exponential;

end
