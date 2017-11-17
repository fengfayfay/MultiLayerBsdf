close all
cd('/Users/mandy/Github/MultiLayerBsdf/build');
alpha = 0.4;
filename = ['outputx_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C1 = textscan(fileID,'%f');
fclose(fileID);

filename = ['outputy_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C2 = textscan(fileID,'%f');
fclose(fileID);

filename = ['outputz_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C3 = textscan(fileID,'%f');
fclose(fileID);

filename = ['outputweight_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C4 = textscan(fileID,'%f');
fclose(fileID);

filename = ['outputdepth_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C5 = textscan(fileID,'%f');
fclose(fileID);


x = C1{1};
y = C2{1};
z = C3{1};
weight = C4{1};
depth = C5{1};

% figure
% scatter3(x,y,z)

% parameter for gaussian heightfield test
N = 4800;
len = 40;
observe = 6000;

% visualization parameters
munum = 100;
phinum = 400;

x1 = linspace(-len,len,N);
y1 = linspace(-len,len,N);
[X,Y] = meshgrid(x1,y1);
% hold on
% surf(X,Y,f)

x = x/observe;
y = y/observe;
z = z/observe;
mu_unit = 1/munum;
phi_unit = 2*pi/phinum;
result = zeros(phinum,munum);
result_d1 = zeros(phinum,munum);
result_d2 = zeros(phinum,munum);
epsilon = 1e-6;
for i = 1:length(x)
    if z(i) >=0
        if abs(x(i))<epsilon && y(i)>=0
            theta = pi/2;
        elseif (abs(x(i))<epsilon && y(i)<0)
            theta = 3*pi/2;
        else
            theta = atan(y(i)/x(i));
            if x(i)< 0 
                theta = theta + pi;
            else
                if y(i) < 0
                theta = theta + 2*pi;
                end
            end
        end
        result(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) = result(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) + weight(i);
        if depth(i)<=1
            result_d1(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) = result_d1(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) + weight(i);
        else
            result_d2(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) = result_d2(ceil(theta/phi_unit),ceil(z(i)/mu_unit)) + weight(i);
        end
    end
end
figure
imagesc(result)
title('reflection lobe')
xlabel('mu_o')
ylabel('phi_o')
colorbar

result2 = zeros(phinum,munum);
result2_d1 = zeros(phinum,munum);
result2_d2 = zeros(phinum,munum);
for i = 1:length(x)
    if z(i)<0
        if abs(x(i))<1e-6 && y(i)>=0
            theta = pi/2;
        elseif (abs(x(i))<1e-6 && y(i)<0)
            theta = 3*pi/2;
        else
            theta = atan(y(i)/x(i));
            if x(i)< 0 
                theta = theta + pi;
            else
                if y(i) < 0
                theta = theta + 2*pi;
                end
            end
        end
        result2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) = result2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) + weight(i);
        if depth(i)<=1
            result2_d1(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) = result2_d1(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) + weight(i);
        else
            result2_d2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) = result2_d2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) + weight(i);
        end
    end
end
figure
imagesc(result2)
title('transmission lobe')
xlabel('mu_o')
ylabel('phi_o')
colorbar

% cd('/Users/mandy/Github/pixar/ritest');
% fid = fopen('scatteredray.txt','w');
% fprintf(fid,'%6f\n',result)
% fclose(fid)

% cd('/Users/mandy/Github/pixar/ritest');
% fid = fopen('transmitray.txt','w');
% fprintf(fid,'%6f\n',result2)
% fclose(fid)

cd('/Users/mandy/Github/pixar/ritest');
filename = ['reflect_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result);
fclose(fid);

filename = ['transmit_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result2);
fclose(fid);

filename = ['reflect1_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result_d1);
fclose(fid);

filename = ['reflect2_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result_d2);
fclose(fid);

filename = ['transmit1_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result2_d1);
fclose(fid);

filename = ['transmit2_', num2str(alpha),'.txt'];
fid = fopen(filename,'w');
fprintf(fid,'%6f\n',result2_d2);
fclose(fid);

% reflect=[reflect result];
% transmit = [transmit result2];

% depth_r = depth(z>0);
% depth_t = depth(z<0);
% yt = 0:0.1:1;
% figure
% histogram(depth,'Normalization','cdf');
% title('accumulated energy percentage from reflectance & transmission, alpha = 0.1')
% set(gca, 'YTick', yt);
% figure
% histogram(depth_r,'Normalization','cdf');
% title('accumulated energy percentage from reflectance, alpha = 0.1')
% set(gca, 'YTick', yt);
% figure
% histogram(depth_t,'Normalization','cdf');
% title('accumulated energy percentage from transmission, alpha = 0.1')
% set(gca, 'YTick', yt);
