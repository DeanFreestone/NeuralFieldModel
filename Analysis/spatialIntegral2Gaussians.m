function [integral] = spatialIntegral2Gaussians(X, Y, NPoints, mu1, sigma1, mu2, sigma2)
%%
stepSize = abs(X(1,1)-X(1,2));

product = zeros(size(X));
for m = 1 : NPoints
    for n = 1 : NPoints
        r = [X(m, n); Y(m, n)];
        
        e1 = exp(-((r - mu1)'/sigma1*(r - mu1)));
        e2 = exp(-((r - mu2)'/sigma2*(r - mu2)));
        
        product(m, n) = e1 * e2;
    end
end


% product = conv2(gaussian1, gaussian2, 'same');

integral = sum(sum(product*stepSize^2, 2), 1);
end