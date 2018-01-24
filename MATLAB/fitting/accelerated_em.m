function obj = accelerated_em(train,trainnum,numGaussian)
    % use accelearted em for fitting
    tic;
    tree = buildtree(train, 0, 0, 3, ceil(trainnum/1000));
    [W,M,R,~,~,~,~] = em(train,[],numGaussian,0,0,tree);
    fprintf('\nRuntime: %.2f seconds\n', toc);
    Rnew = reshape(R', size(train,2),size(train,2),numGaussian);
    obj = gmdistribution(M,Rnew,W');
end