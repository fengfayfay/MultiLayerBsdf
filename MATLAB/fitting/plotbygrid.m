function [count, A] = plotbygrid(xnum,ynum,Y,...
    range,titlestring,filename)

% calculate predict matrix
    x_unit = 2.0 * range/xnum;
    y_unit = 2.0 * range/ynum;

    index_X0 = xnum * .5;
    index_Y0 = ynum * .5;
     
    A = zeros(xnum,ynum);
    count = 0;
    for i = 1:length(Y)
        xindex = ceil(Y(i, 1) / x_unit + .5) + index_X0;
        yindex = ceil(Y(i, 2) / y_unit + .5) + index_Y0;
        
        if xindex < xnum + 1 && yindex < ynum + 1
            A(xindex, yindex) = A(xindex, yindex) + 1;
        else
            count = count +1;
        end
    end

    [maxAR, maxRIndex] = max(A);
    [maxA, maxCIndex] = max(maxAR);
    
    maxX = ((maxRIndex(maxCIndex) - index_X0) + .5) * x_unit;
    maxY = ((maxCIndex - index_Y0) + .5)  * y_unit;

    disp([maxCIndex, maxRIndex(maxCIndex)]);

    %disp(maxA);
    %disp(maxIndex);
    disp([maxX, maxY]);
    disp(count);
     
    figure
    imagesc(A/sum(sum(A)))
    colorbar()
    ylabel('x/z')
    xlabel('y/z')
    title(titlestring)
    saveas(gcf,[filename,'.jpeg'])

end
