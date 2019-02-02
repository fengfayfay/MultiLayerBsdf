%
% plot ellipsoid representing the mixture of gaussians
% without plotting their weights
%
function plotGMM(obj, componentInterval)

    hold on
    numGaussian = obj.NumComponents;
    numVar = obj.NumVariables;
    if componentInterval == 0
        componentInterval = numGaussian;
    end
    plotCount = floor(numGaussian/componentInterval);
    disp(plotCount)
    for j = 0:plotCount-1
        figure
        for i = 1:componentInterval
            index = j * componentInterval + i;
            plot_gaussian_ellipsoid(obj.mu(index,:), obj.Sigma(:,:,index))
        end
        if numVar == 3
            view(129,36); set(gca,'proj','perspective'); grid on;
            grid on; axis equal; axis tight;
        end
    end
end