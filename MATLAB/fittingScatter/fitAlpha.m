function [energy, mean0, mean1, cv0, cv1, mu, angleValue, alphaValues] = fitAlpha()
    %alphaValues = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    %alphaValues = [0.9 0.8 0.7 0.6 0.5];
    alphaValues = [.9 0.8 0.7 0.6 0.5 0.4];
    mean0 = [];
    mean1 = [];
    cv0 = [];
    cv1 = [];
    energy = [];
    mu = [];
    angleValue = [];
    for k = 1:length(alphaValues)
        alpha = alphaValues(k);
        [mean_0, mean_1, cv_0, cv_1, energyRatios, anglevalues, errs, mus] = OneLayerMu(alpha, 4, 3);
        mean0 = [mean0 mean_0];
        mean1 = [mean1 mean_1];
        cv0 = [cv0 cv_0];
        cv1 = [cv1 cv_1];
        energy = [energy energyRatios];
        mu = [mu mus];
        angleValue = [angleValue anglevalues];
    
        %figure out how to quickly do double array
        %also need to zero out values for energy, mean and cvs when energy is < 1e-3. Otherwise different 
        %alpha arrays would have different sizes which would cause problems with alpha plot     
    end
end 
