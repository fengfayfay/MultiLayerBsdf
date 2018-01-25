function [count, A] = plotbygrid(xnum,ynum,Y,...
    range,x_unit,y_unit,titlestring,filename)

% calculate predict matrix
A= zeros(xnum,ynum);
count = 0;
for i = 1:length(Y)
    if abs(Y(i,1))<range && abs(Y(i,2))<range
        A(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) = ...
            A(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) + 1;
    else
        count = count +1;
    end
end

figure
imagesc(A/sum(sum(A)))
colorbar()
ylabel('x/z')
xlabel('y/z')
title(titlestring)
saveas(gcf,[filename,'.jpeg'])

end