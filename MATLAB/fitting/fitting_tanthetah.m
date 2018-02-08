function[train, result, obj, Y, predict] =  fitting_tanthetah(dir,fundir,alpha,angle,input,trainnum, ...
    generatenum, gaussiannumvec, tanthetanum,accelerated,maxiter,tol)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% Energy vs half vector from data
incident = [sin(angle*pi/180),0,cos(angle*pi/180)];
h = (input + incident)/2;
hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
h = h./hnorm;
tanthetavec = sqrt((1-h(:,3).^2)./h(:,3).^2);
tanthetarange = max(tanthetavec);
fprintf('tanthetah range is %4.4f\n',tanthetarange)
thetaunit = tanthetarange/tanthetanum;

train = zeros(trainnum,1);
for i = 1:trainnum
      train(i) = tanthetavec(i);
end

result = zeros(tanthetanum,1);
indexrange = trainnum+1:trainnum+generatenum;
for i = indexrange
    tancur = tanthetavec(i);
    index = ceil(tancur/thetaunit);
    result(index) =  result(index)+1;
end
result = result/sum(result(:));

figure
plot(linspace(0,tanthetarange,tanthetanum),result/sum(result),'linewidth',2)
xlabel('tanthetah')
hold on

distName = 'Loglogistic';
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using tan of half vector
    numGaussian = gaussiannumvec(j);
    
    if accelerated
        aemdir = [fundir,'accelerated_greedy_EM'];
        addpath(aemdir);
        obj = accelerated_em(train,trainnum,numGaussian,maxiter,tol);
    else
        %obj = customized_em(train,numGaussian,maxiter,tol);
        obj = fitdist(train, distName);
    end
  
    filename = [dir,'tanhalf',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    %% calculate energy
    predict = zeros(tanthetanum,1);
    count = 0;
    for i = 1:generatenum
        if Y(i)>=0 && Y(i)<tanthetarange
            predict(ceil(Y(i)/thetaunit)) =  predict(ceil(Y(i)/thetaunit))+1;
        else
            count = count+1;
        end
    end
    predict = predict/sum(predict);
    plot(linspace(0,tanthetarange,tanthetanum),predict/sum(predict),'linewidth',2)
    legend('input','fitting')
    
    %% calculate error
    
    [err,~,~,~] = relativel2err(result,predict);
    fprintf('err: %4.6f\n', err)
    
    errvec(j) = err;
    countvec(j) = count;
    
end

errvec_filename = [dir,'tanhalf',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'tanhalf',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')

end
