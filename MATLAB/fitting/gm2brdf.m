

function gm2brdf(obj, dim, angle, alpha, extendratio, filename, energyRatio)
%
% convert from slope domain probability to brdf*cos value
%

phinum = 400;
phi = 2*pi*((0:1:phinum-1)+0.5)/phinum;

munum = 100;
mu = ((0:1:munum-1)+0.5)/munum;
scale = 2 * pi/(munum * phinum) * energyRatio;

anglevec = angle;
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
            detJ = computeJacobian(wi, wo);
            brdfcos(i,j) = p * detJ * (1 + extendratio*2) * scale;
        end
    end
   
    figure
    imagesc(brdfcos)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['brdf*cos angle=',num2str(angle), ' alpha=', num2str(alpha)])
    saveas(gcf,[filename,'.jpeg'])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:)));
end
end

function [jacobian] = computeJacobian(wo, wi)
    denom = (wo(3) + wi(3));
    denom = abs(denom * denom * denom);
    jacobian = abs(wo(3)* wi(3)  + wo(2)*wi(2) + wo(1)*wi(1) + 1)/denom;
end 
