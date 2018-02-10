function [obj, W, M, R] = customized_em(train,numGaussian,maxiter,tol,softinit, oW, oM, oR)
    % use matlab gmcluster for em fitting
    tic
    options = statset('MaxIter',maxiter, 'Display','final','TolFun',tol);
    try
        if softinit
            S = struct('ComponentProportion',oW,'mu',oM,'Sigma',oR);
            obj = fitgmdist(train,numGaussian,'Options',options,'Start',S);
        else
            obj = fitgmdist(train,numGaussian,'Options',options);
        end
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message;
    end
    W = obj.ComponentProportion;
    M = obj.mu;
    R = obj.Sigma;
    fprintf('\nRuntime: %.2f seconds\n', toc);
end