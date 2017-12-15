close all
dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror/pi:3/output';
cd(dir)

anglevec = [60];
trainnum = 1e6;
testafter = 1e6;

for j = 1:length(anglevec)
    angle = anglevec(j);
    alpha = 0.5;
    filename = [num2str(angle), 'outputx_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C1 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputy_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C2 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputz_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C3 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputweight_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C4 = textscan(fileID,'%f');
    fclose(fileID);
    
    x = C1{1};
    y = C2{1};
    z = C3{1};
    weight = C4{1};
    observe = 6000;
    x = x/observe;
    y = y/observe;
    z = z/observe;
    epsilon = 1e-6;
    
%     %% half vector
%     xynum = 200;
%     znum = 100;
%     % Energy value from data
%     xy_unit = 2*sqrt(2)/xynum;
%     z_unit = 2/znum;
%     result = zeros(xynum,znum);
%     for i = testafter+1:length(x)
%         if z(i)>=0
%             h = ([x(i),y(i),z(i)] + incident)/2;
%             h = h/norm(h);
%             
%             dim1 = sqrt(h(1)^2 + h(2)^2);
%             dim2 = h(3);
%             result(ceil((-dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) = result(ceil((-dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) + weight(i);
%             result(ceil((dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) = result(ceil((dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) + weight(i);
%         end
%     end
%     
%     figure
%     imagesc(result/generatenum)
%     colorbar()
    
    %% phi theta of outgoing direction
    thetanum = 400;
    phinum = 400;
    theta_unit = 2*pi/thetanum;
    phi_unit = 2*pi/phinum;
%     result = zeros(thetanum,phinum);
    result = zeros(thetanum,1);
    for i = 1:length(x)
        if abs(y(i)) < phi_unit/2
            theta = acos(z(i));
            if x(i) >0
                theta = -theta;
            end
            result(ceil((theta+pi)/theta_unit)) = result(ceil((theta+pi)/theta_unit)) + weight(i);
        end
%         if abs(x(i))<1e-6 && y(i)>=0
%             phi = pi/2;
%         elseif (abs(x(i))<1e-6 && y(i)<0)
%             phi = 3*pi/2;
%         else
%             phi = atan(y(i)/x(i));
%             if x(i)< 0 
%                 phi = phi + pi;
%             else
%                 if y(i) < 0
%                 phi = phi + 2*pi;
%                 end
%             end
%         end
%         theta = acos(z(i));
% %         if pi > pi/2 && theta < 3*pi/2
% %             theta = -theta;
% %         end
%         result(ceil(phi/phi_unit), ceil((theta+pi)/theta_unit)) = result(ceil(phi/phi_unit), ceil((theta+pi)/theta_unit))+ weight(i);
    end
    
%     figure
%     imagesc(result/sum(weight))
%     colorbar()
%     thetavec = linspace(-pi,pi,thetanum);
%     figure
%     plot(thetavec*180/pi,result(200,:)/sum(weight),'linewidth',2)
      
       figure
       plot(thetavec*180/pi,result/sum(weight),'linewidth',2)
    
end