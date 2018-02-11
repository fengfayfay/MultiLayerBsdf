function gm2brdf(obj,dim,alpha,extendratio)
%
% convert from slope domain probability to brdf*cos value
%

phinum = 360;
phi = 2*pi*((0:1:phinum-1)+0.5)/phinum;

munum = 100;
mu = ((0:1:munum-1)+0.5)/munum;


anglevec = 60;
for k = 1:length(anglevec)
    angle = anglevec(k);
    theta = angle*pi/180;
    wi = [sin(theta), 0, cos(theta)];
    brdfcos = zeros(phinum,munum);
    
    for i = 1:phinum
        for j = 1:munum
            sintheta = sqrt(1 - mu(j)^2);
            wo = [sintheta*cos(phi(i)), sintheta*sin(phi(i)), mu(j)];
            h1 = (wi + wo)/2;
            h = h1/norm(h1);
            
            if dim==2
                p = pdf(obj,[h(1)/h(3),h(2)/h(3)]);
            else
                p = pdf(obj,[h(1)/h(3),h(2)/h(3),theta]);
                %conditioned on this incident angle
                p = p/(1/(pi/2));
            end
            
            %Jacobian
            detJ = (mu(j) * wi(3) + wo(2)*wi(2) + wo(1)*wi(1) + 1)/(wi(3)+mu(j))^3;
            brdfcos(i,j) = p * detJ * (1 + extendratio*2);
        end
    end
    
    figure
    imagesc(brdfcos)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['brdf*cos angle=',num2str(angle), ' alpha=', num2str(alpha)])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:))*2*pi/(phinum*munum));

    
    
%     %% pbrt output
%     dir = '/Users/mandy/Github/MultiLayerBsdf/build/';
%     filename = [dir,num2str(angle),'brdfcos.txt'];
%     % filename = [dir,num2str(angle),'brdfcos_reflect.txt'];
%     fileID = fopen(filename);
%     C1 = textscan(fileID,'%f');
%     fclose(fileID);
%     
%     brdfcos1 = C1{1};
%     brdfcos1 = reshape(brdfcos1,phinum,munum);
%     
%     figure
%     imagesc(brdfcos1)
%     colorbar()
%     xlabel('mu')
%     ylabel('phi')
%     title(['brdf*cos angle=',num2str(angle), ' pbrt'])
%     
%     diff = abs(brdfcos - brdfcos1);
%     figure
%     imagesc(diff)
%     colorbar()
%     xlabel('mu')
%     ylabel('phi')
%     title(['brdf*cos diff, angle=',num2str(angle)])
end


% %% plot brdf*cos on phi(0,2*pi) theta(-pi,pi) domain 
% thetanum = 360;
% thetarange = 2*pi;
% thetavec = thetarange * ((0:1:thetanum-1)+0.5)/thetanum - thetarange/2;
% 
% % wi is incident
% angle = 89;
% theta_i = angle*pi/180;
% wi = [sin(theta_i), 0, cos(theta_i)];
% brdfcos = zeros(phinum,munum);
% 
% for i = 1:phinum
%     for j = 1:thetanum
%         wo = [sin(thetavec(j))*cos(phi(i)), sin(thetavec(j))*sin(phi(i)), cos(thetavec(j))];
%         h1 = (wi + wo)/2;
%         h = h1/norm(h1);
%         
%         if dim==2
%             p = pdf(obj,[h(1)/h(3),h(2)/h(3)]);
%         else
%             p = pdf(obj,[h(1)/h(3),h(2)/h(3),theta_i]);
%             % conditioned on this incident angle
%             p = p/(1/(pi/2));
%         end
%         
%         % Jacobian
%         detJ = (cos(thetavec(j)) * wi(3) + wo(2)*wi(2) + wo(1)*wi(1) + 1)/(wi(3)+cos(thetavec(j)))^3;
%         brdfcos(i,j)= (reflect+1)*p*detJ;
%     end
% end
% 
% figure
% imagesc(brdfcos)
% colorbar()
% xlabel('theta')
% ylabel('phi')
% title(['brdf*cos angle=',num2str(angle), ' alpha=', num2str(alpha)])
% 
% fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:))*2*pi/(phinum*munum));
end