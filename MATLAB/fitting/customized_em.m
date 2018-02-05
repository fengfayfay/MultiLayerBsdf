function obj = customized_em(train,numGaussian,maxiter,tol)
    % use matlab gmcluster for em fitting
    tic
    options = statset('MaxIter',maxiter, 'Display','final','TolFun',tol);
    try
%         obj = fitgmdist(train,numGaussian,'Options',options,'Start','customize');
        obj = fitgmdist(train,numGaussian,'Options',options);
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message;
    end
    fprintf('\nRuntime: %.2f seconds\n', toc);
end