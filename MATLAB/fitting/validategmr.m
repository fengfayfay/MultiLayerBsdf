%
%   plot gaussian mixture regression result and validate
%

generatenum=1e6;
for j = 1:9
    obj = gmdistribution(expData(1:2,(j-1)*10+1)',expSigma(:,:,(j-1)*10+1),1);
    Y = random(obj,generatenum);
    
    %    calculate predict matrix
    predict= zeros(xnum,ynum);
    count = 0;
    for i = 1:generatenum
        if abs(Y(i,1))<range && abs(Y(i,2))<range
            predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) = ...
                predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) + 1;
        else
            count = count +1;
        end
    end
    
    %     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/sum(sum(predict)))
    colorbar()
end