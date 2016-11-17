function integral = spatialIntegral2Gaussians(X, Y, NPoints, mu1, sigma1, mu2, sigma2)
%% Integrate two 2-D gaussians over a space
% Miao Cao
% ~~~~~~~~~~~~~~~

% Parameter list:
% X - x coordinates of the field
% Y - y coordinates of the field
% NPoints - number of points in each row or column of the field
% mu1 - centre of the first Gaussian
% sigma1 - variance-covariance matrix of the first Gaussian
% mu2 - centre of the second Gaussian
% sigma2 - variance-covariance matrix of the second Gaussian


%% Compute the integration
% ~~~~~~~~~~~~~~~


stepSize = abs(X(1,1)-X(1,2)); % distance or step size of spatial discretisation

product = zeros(size(X)); % product of two Gaussians

for m = 1 : NPoints
    
    for n = 1 : NPoints
        
        r = [X(m, n); Y(m, n)];
        
        e1 = exp(-((r - mu1)'/sigma1*(r - mu1)));
        e2 = exp(-((r - mu2)'/sigma2*(r - mu2)));
        
        product(m, n) = e1 * e2;
        
    end
end

integral = sum(sum(product*stepSize^2, 2), 1); % integrate over the space

end
