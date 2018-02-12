%
% plot ellipsoid representing the mixture of gaussians
% without plotting their weights
%

numGaussian = obj.NumComponents;
figure
hold on
for i = 1:numGaussian
    plot_gaussian_ellipsoid(obj.mu(i,:), obj.Sigma(:,:,i))
end