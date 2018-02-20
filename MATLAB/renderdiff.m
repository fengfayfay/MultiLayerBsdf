close all
alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.8];
figure
for i = 6

I1 = imread(['../build/killeroo-gauss_b+g.png']);
I2 = imread(['../build/killeroo-beck.png']);
I1 = im2double(I1);
I2 = im2double(I2);
diff = I1-I2;
disp(sum(sum(diff(:,:,1)>0))/sum(sum(diff(:,:,1)~=0)))
disp(sum(sum(diff(:,:,1)<0))/sum(sum(diff(:,:,1)~=0)))
disp(sum(sum(diff(:,:,1))))
% subplot(2,3,i);
imagesc(diff(:,:,1))
colorbar()

end