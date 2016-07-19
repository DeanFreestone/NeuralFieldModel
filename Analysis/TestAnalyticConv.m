%%
clear
clc
close all


rx= linspace(-1,1,100);
ry= linspace(-1,1,100);
r = meshgrid(rx, ry);

mu_i = [0.5 0.5];
mu_j = [-0.5 -0.5];

mu = mu_i + mu_j;

sigma_i = 0.1;
sigma_j = 0.1;

x = (pi*sigma_i^2*sigma_j^2) / (sigma_i^2 + sigma_j^2);
y = exp(-(r'*r)/(sigma_i^2+sigma_j^2));


imagesc(x*y)
