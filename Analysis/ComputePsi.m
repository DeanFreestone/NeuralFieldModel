function Ts_invGamma_phi_psi = ComputePsi(SpaceMin, SpaceMax, NPoints, nTheta, Ts, nx, mu_phi, sigma_phi, mu_psi, vector_Sigma_Psi)
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

%% parameters and centres of field basis functions are pre-defined here
% We align the centres of field basis function to either right on the
% discretisation points or the middle point of two adjacent discretisation
% points
% ~~~~~~~~~~~~~~~
% calc location of Gaussian basis functions. Uniformly distributed over 2-D surface.



x = linspace(SpaceMin, SpaceMax, NPoints*2-1); % coordinates of discretisation of the surface

numRow = sqrt(nx); % number of gaussians for each colomn
numCol = nx / numRow; % number of columns

widthSpace = SpaceMax - SpaceMin;
widthCentre = widthSpace / (numCol*2);

if isempty(mu_phi) || ~any(mu_phi(:)) % if mu_phi is empty or only zeros
    mu_phi = zeros(nx, 2); % centres of each phi (gaussian)
    for m = 1 : numRow
        for n = 1 : numCol
            r = [(SpaceMin - widthCentre + m*widthCentre*2) (SpaceMin - widthCentre + n*widthCentre*2)]; % location if uniform distribution
            rX = r(1); rY = r(2);
            [~, indX] = min(abs(x - rX)); [~, indY] = min(abs(x - rY));
            rX = x(indX); rY = x(indY); % aligned with discretisation
            mu_phi(n + numCol*(m-1), :) = [rX, rY];
        end
    end
end

covMat_phi = sigma_phi; % covariance matrix of phi

%% Compute gamma
% ~~~~~~~~~~~~~~~

Gamma = ComputeGamma(SpaceMin, SpaceMax, NPoints, nx, mu_phi, covMat_phi); % compute gamma based on phi

%% Compute Psi - analytic
% ~~~~~~~~~~~~~~~
% now form the matrix
% these are the coefficients for the analytic convolution of psi and phi
% But, we haven't figure out covariance matrix here.



% compute coeffcients
psi_phi_coefficient = zeros(length(vector_Sigma_Psi), 1);

for m = 1 : nTheta
    
    % covariance matrix of this basis function
    covMat_Psi = [vector_Sigma_Psi(m, 1) vector_Sigma_Psi(m, 2); vector_Sigma_Psi(m, 2) vector_Sigma_Psi(m, 1)];
    
    % coefficient of convolution of phi and psi basis functions
    %     psi_phi_coefficient(m) = pi*det(covMat_Psi)*det(covMat_phi) / det(covMat_Psi+covMat_phi); % coefficient of convolution of phi and psi basis functions
    psi_phi_coefficient(m) = 2*pi*sqrt(det(covMat_Psi))*sqrt(det(covMat_phi)) / sqrt(det(covMat_Psi+covMat_phi)); % coefficient of convolution of phi and psi basis functions
    
end


% compute the convolution of phi and psi
psi_phi_basis = zeros(nTheta, nx, NPoints, NPoints); % nx * ntheta * fields

for m=1 : nTheta % cycle through each connectivity basis function
    
    for n=1 : nx % cycle through each field basis function
        
        mu = mu_phi(n, :) + mu_psi(m, :) + 2*mu_psi(m, :); % centre of a Gaussian after convolution of phi and psi
        
        % covariance matrix of this basis function of connectivity kernel
        covMat_Psi = [vector_Sigma_Psi(m, 1) vector_Sigma_Psi(m, 2); vector_Sigma_Psi(m, 2) vector_Sigma_Psi(m, 1)];
        
        psi_phi = psi_phi_coefficient(m)*Define2DGaussian_AnisotropicKernel(mu(1), mu(2), covMat_phi + covMat_Psi, NPoints, SpaceMin, SpaceMax); % convolution of phi ans psi is another Gaussian
        
        psi_phi_basis(m, n, :, :) = psi_phi(:, :);
        
        
    end
    
end

Ts_invGamma_phi_psi = zeros(nTheta, nx, NPoints, NPoints); % initialise the matrix of fields. Psi matrix. Dimensions: nx * ntheta * fields

inv_Gamma = inv(Gamma); % inverse Gamma matrix

for m = 1 : nTheta % cycle through each row of the matrix of fields
    
    fieldVector = squeeze(psi_phi_basis(m, :, :, :));
    
    inv_Gamma_fieldVector = zeros(size(fieldVector));
    
    for p = 1 : nx % cycle through each field basis function
        
        for q = 1 : nx % cycle through each field basis function
            
            inv_Gamma_fieldVector(p, :, :) = inv_Gamma_fieldVector(p, :, :) + inv_Gamma(p, q) .* fieldVector(q, :, :); % calculate each element in the Psi matrix
            
        end
        
    end
    
    Ts_invGamma_phi_psi(m,:,:,:) = Ts * inv_Gamma_fieldVector; % Psi matrix
end

end
