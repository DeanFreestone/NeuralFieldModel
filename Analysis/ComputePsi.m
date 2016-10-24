function Ts_invGamma_phi_psi = ComputePsi(X, Y, SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi)
%% Compute Psi
% Compute Psi in Equation (24), Freestone et al., 2011, NeuroImage
% Miao Cao

% ~~~~~~~~~~~~~~~
% input parameter list:
% parameters used to create a  2-D cortical surface
% SpaceMin - Edge of surface on the negative side
% SpaceMax - Edge of surface on the positive side
% NPoints - number of points along each dimension

% nTheta - number of connectivity kernel basis functions
% Ts - Time step
% nx - number of Gaussian basis functions
% mu_phi - centre of Gaussian basis function for phi
% sigma_phi - sigma of Gaussian basis function for phi
% mu_psi - a vector of centres of Gaussian basis functions of connectivity kernel
% vector_Sigma_Psi - a vector of sigma of Gaussian basis functions

% ~~~~~~~~~~~~~~~
% output parameter list:
% Psi - Psi matrix

%% parameters and variables are pre-defined here
% ~~~~~~~~~~~~~~~
% calc location of Gaussian basis functions. Uniformly distributed over 2-D surface.


numRow = sqrt(nx); % number of gaussians for each colomn
numCol = nx / numRow; % number of columns

widthSpace = SpaceMax - SpaceMin;
widthCentre = widthSpace / (numCol*2);

if isempty(mu_phi) || any(mu_phi) % if mu_phi is empty or only zeros
    mu_phi = zeros(nx, 2); % centres of each phi (gaussian)
    for m = 1 : numRow
        for n = 1 : numCol
            mu_phi(n + numCol*(m-1), :) = [(SpaceMin - widthCentre + m*widthCentre*2) (SpaceMin - widthCentre + n*widthCentre*2)];
        end
    end
end

covMat_phi = [sigma_phi 0; 0 sigma_phi]; % covariance matrix of phi

%% Compute gamma
% ~~~~~~~~~~~~~~~

Gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu_phi, covMat_phi); % compute gamma based on phi

%% Compute Psi - analytic
% ~~~~~~~~~~~~~~~
% now form the matrix
% these are the coefficients for the analytic convolution of psi and phi
% But, we haven't figure out covariance matrix here.


psi_phi_coefficient = zeros(length(vector_Sigma_Psi), 1);
for m = 1 : length(vector_Sigma_Psi)
    psi_phi_coefficient(m) = pi*vector_Sigma_Psi(m)*sigma_phi / (vector_Sigma_Psi(m)+sigma_phi);
end

% compute the convolution between phi and psi
psi_phi_basis = zeros(nTheta, nx, NPoints, NPoints); % nx * ntheta * fields
for m=1 : nTheta
    for n=1 : nx
        % these guys here are used with the LS algorithm for estimating
        % theta and xi
        
        mu = mu_phi(n, :) + mu_psi(m, :) + 2*mu_psi(m, :); % centre of Gaussian after convolution
        
        psi_phi = psi_phi_coefficient(m)*Define2DGaussian_AnisotropicKernel(mu(1), mu(2), [vector_Sigma_Psi(m) 0; 0 vector_Sigma_Psi(m)]+covMat_phi, NPoints, SpaceMin, SpaceMax);
        
        psi_phi_basis(m, n, :, :) = psi_phi(:, :);
        
        
        %         % cycle through every point on the cortical surface
        %         for p = 1 : NPoints
        %             for q = 1 : NPoints
        %
        %                 rPrime = [X(p, q) Y(p, q)]; % location r'
        %
        %                 mu = mu_phi(n, :) + mu_psi(m, :) + 2*mu_psi(m, :) + rPrime; % centre of Gaussian after convolution
        %
        %                 psi_phi = psi_phi_coefficient(m)*Define2DGaussian_AnisotropicKernel(mu(1), mu(2), [vector_Sigma_Psi(m) 0; 0 vector_Sigma_Psi(m)]+covMat_phi, NPoints, SpaceMin, SpaceMax);
        %
        %                 psi_phi_basis(m, n, :, :) = squeeze(psi_phi_basis(m, n, :, :)) + psi_phi(:, :);
        %
        %             end
        %         end
        %         theta_psi_phi_basis(nn,n,:) = theta(nn)*psi_phi_basis(nn,n,:);
    end
    
end

Ts_invGamma_phi_psi = zeros(nTheta, nx, NPoints, NPoints); % initialise the matrix of fields. nx * ntheta * fields

inv_Gamma = inv(Gamma);

for m = 1 : nTheta % cycle through each row of the matrix of fields
    
    fieldVector = squeeze(psi_phi_basis(m, :, :, :)); % a
    inv_Gamma_fieldVector = zeros(size(fieldVector));
    
    for p = 1 : nx
        for q = 1 : nx
            inv_Gamma_fieldVector(p, :, :) = inv_Gamma_fieldVector(p, :, :) + inv_Gamma(p, q) .* fieldVector(q, :, :);
        end
    end
    
    Ts_invGamma_phi_psi(m,:,:,:) = Ts * inv_Gamma_fieldVector;
end

end
