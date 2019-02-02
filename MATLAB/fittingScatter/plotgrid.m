function [count, A] = plotgrid(input,xnum, ynum, titlestring,filename, energyRatio)

xrange = [0 2*pi];
yrange = [0 1];

% calculate predict matrix
ynum = 100;
xnum = 400;

x_unit = (xrange(2) - xrange(1))/(xnum-1);
y_unit = (yrange(2) - yrange(1))/(ynum-1);
A= zeros(xnum,ynum);
count = 0;
for i = 1:length(input)
    phi = atan2(input(i,2), input(i,1));
    if phi < 0
        phi = phi+ 2 * pi;
    end

    x = ceil(phi/x_unit + .5);
    y = ceil(input(i, 3)/y_unit + .5);
    if (x <= xnum && y <=ynum) 
        A(x, y) = A(x, y) + 1;
    else
        count = count +1;
    end
end

B = A/sum(sum(A));
B = log(B * energyRatio);

figure
imagesc(B)
colorbar()
ylabel('phi')
xlabel('mu')
title(titlestring)
saveas(gcf,[filename,'.jpeg'])

end
