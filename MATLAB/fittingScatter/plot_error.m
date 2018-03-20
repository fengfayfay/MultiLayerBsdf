function [fitresult, gof] = plot_error(anglevalues, errs, resultDir)
%CREATEFIT1(ANGLEVALUES,ERRS)
%  Create a fit.
%
%  Data for 'L2 error of gaussian fit' fit:
%      X Input : anglevalues
%      Y Output: errs
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Mar-2018 09:39:30


%% Fit: 'L2 error of gaussian fit'.
[xData, yData] = prepareCurveData( anglevalues, errs );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'L2 error of gaussian fit' );
h = plot( fitresult, xData, yData );
legend( h, 'errs vs. anglevalues', 'L2 error of gaussian fit', 'Location', 'NorthEast' );
% Label axes
xlabel anglevalues
ylabel errs
grid on

filename = [resultDir, 'gaussian_error'];
saveas(gcf,[filename,'.jpeg'])


