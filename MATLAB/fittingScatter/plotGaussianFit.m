function [fitresult, gof] = plotGaussianFit(polyDegree, anglevalues, target, targetName, targetFileName, resultDir, isMu)
%CREATEFIT(ANGLEVALUES,CV1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : anglevalues
%      Y Output: target
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Mar-2018 09:16:15


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( anglevalues, target );

% Set up fittype and options.
fitname = ['poly', num2str(polyDegree)];
ft = fittype( fitname );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figname = [targetName, ' fit']; 
myFig = figure( 'Name', 'Poly Fit');
h = plot( fitresult, xData, yData );
% Label axes
if isMu
    xlabel('\mu', 'FontSize', 20)
    filename = [resultDir, 'polyfit_mu_', targetFileName]
else
    xlabel('\theta_i', 'FontSize', 20)
    filename = [resultDir, 'polyfit_', targetFileName]
end
ylabel(targetName, 'FontSize', 20)
grid on

legend('data', ['fit RMSE ', num2str(gof.rmse)] );
title(['\fontsize{20}', 'Polynomial Fit of ', targetName]);
set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal')

saveas(gcf,[filename,'.jpeg']);

end
