function fitting_tanhalfvector(dir,alpha,angle,x,y,z,weight,testafter, trainnum, ...
    generatenum, gaussiannumvec, incident, thetarange, thetanum)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% Energy vs half vector from data
thetaunit = thetarange/thetanum;
% tancur = zeros(length(x)-testafter,1);
thetaresult = zeros(thetanum,1);
for i = testafter+1:length(x)
    if z(i)>0
        h = ([x(i),y(i),z(i)] + incident)/2;
        h = h/norm(h);
%         tancur(i-testafter) = sqrt((1-h(3)^2)/h(3)^2);
        tancur = sqrt((1-h(3)^2)/h(3)^2);
        if ceil(tancur/thetaunit)==0
            index = 1;
        else
            index = ceil(tancur/thetaunit);
        end
        thetaresult(index) =  thetaresult(index)+weight(i);
    end
end

% thetarange = max(tancur);
% thetarange_filename = ['tanhalf',num2str(angle),'_angle_',num2str(alpha),'_thetarange.mat'];
% save(thetarange_filename,'thetarange')
% thetaunit = thetarange/thetanum;
% 
% thetaresult = zeros(thetanum,1);
% for i = 1:length(tancur)
%     if ceil(tancur(i)/thetaunit)==0
%         index = 1;
%     else
%         index = ceil(tancur(i)/thetaunit);
%     end
%     thetaresult(index) =  thetaresult(index)+weight(i);
% end

close all
figure
hold on
plot(linspace(0,thetarange,thetanum),thetaresult/generatenum)
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using tan of half vector
    numGaussian = gaussiannumvec(j);
    tantheta = zeros(trainnum,1);
    for i = 1:trainnum
        if z(i)>=0
            h = ([x(i),y(i),z(i)] + incident)/2;
            h = h/norm(h);
            tantheta(i) = sqrt((1-h(3)^2)/h(3)^2);
        end
    end
    train = tantheta;
    options = statset('MaxIter',500, 'Display','final');
    obj = gmdistribution.fit(train,numGaussian,'Options',options);
  
    filename = ['tanhalf',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    %% calculate energy
    predict = zeros(thetanum,1);
    count = 0;
    for i = 1:generatenum
        if Y(i)>=0 && Y(i)<thetarange
            predict(ceil(Y(i)/thetaunit)) =  predict(ceil(Y(i)/thetaunit))+1;
        else
            count = count+1;
        end
    end
    plot(linspace(0,thetarange,thetanum),predict/(generatenum-count))
    
    %% calculate error
    diff = predict/(generatenum-count)-thetaresult/generatenum;
    err = sqrt(sum(diff.*diff));
    
    errvec(j) = err;
    countvec(j) = count/generatenum;

    
end
title(['Energy generated using GMM on tan(htheta), alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
legendInfo{1} = ['heightfield'];
for i = 1:length(gaussiannumvec)
    legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
end
legend(legendInfo)
filename=['tanhalf',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
saveas(gcf,[filename,'.jpeg'])

errvec_filename = ['tanhalf',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = ['tanhalf',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')
end