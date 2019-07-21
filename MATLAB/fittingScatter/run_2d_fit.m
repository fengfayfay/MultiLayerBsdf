function run_2d_fit()

%alphas = [0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20]
alphas = [0.90, 0.80, 0.70, 0.60, 0.50]
for i = 1:length(alphas)
    [mean_0, mean_1, cv_0, cv_1, energyRatios, anglevalues, errs, mus] =  OneLayerMu(1.0, 1.0, alphas(i), 3, 3, 3.5);
end
