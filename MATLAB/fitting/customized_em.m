function obj = customized_em(maxiter,tol,train,numGaussian)
    % use matlab gmcluster for em fitting
    tic
    options = statset('MaxIter',maxiter, 'Display','final','TolFun',tol);
    try
        obj = fitgmdist(train,numGaussian,'Options',options,'Start','customize');
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message;
    end
    fprintf('\nRuntime: %.2f seconds\n', toc);
end