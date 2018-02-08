function [obj, W, M, R, isigma] = accelerated_em(train,trainnum,numGaussian,maxiter,tol, softinit, oW, oM, oR, osigma)
    % use accelearted em for fitting
    tic;
    tree = buildtree(train, 0, 0, 3, ceil(trainnum/1000));
    [W,M,R,isigma,~,~,~,~,iter] = em(train,[],numGaussian,0,0,tree,maxiter,tol, softinit, oW, oM, oR, osigma);
    fprintf('\nRuntime: %.2f seconds\n', toc);
    fprintf('\niterations: %4d \n', iter);
    Rnew = reshape(R', size(train,2),size(train,2),numGaussian);
    obj = gmdistribution(M,Rnew,W');
end
