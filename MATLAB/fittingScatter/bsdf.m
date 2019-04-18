
function bsdf(angle, alpha, filename, ni, no, energyRatio)
% convert from slope domain probability to brdf*cos value
%

phinum = 400;
phi = 2*pi*((0:1:phinum-1)+0.5)/phinum;

munum = 100;

if no ~= ni
    signmu = -1;
else
    signmu = 1;
end
mu = ((0:1:munum-1)+0.5)/munum * signmu;
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
            wo = [sintheta*cos(phi(i)), sintheta*sin(phi(i)), mu(j) * signmu];
            
            h1 = signmu* (ni * wi + no * wo);
            h = h1/norm(h1);
           
            D = BeckmannDistribution(h, alpha);
            G = GVCavity(wi, wo, h);
            J = abs(dot(wo, h) * dot (wi, h)) * no * no;
            denom = ni *  dot(wi, h) + no * dot(wo, h);
            denom = denom * denom * wi(3);
            brdfcos(i,j) = D * G * J/denom * scale;
        end
    end
   
    figure
    imagesc(brdfcos)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['Single Scattering BSDF,  angle=',num2str(angle), ' alpha=', num2str(alpha)])
    saveas(gcf,[filename,'.jpeg'])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:)));
end
end

function [D] = BeckmannDistribution(wh, alpha)
    cosTheta = wh(3);
    cos2Theta = cosTheta * cosTheta;
    sin2Theta = 1.0 - cos2Theta; 
    tan2Theta = sin2Theta/cos2Theta;

    cos4Theta = cos2Theta * cos2Theta;
    D =  exp(-tan2Theta / (alpha * alpha))  / (pi * alpha * alpha * cos4Theta);
    
end

function [G] = GVCavity(wo, wi, wh)
    GO = GVCavity1(wo, wh);
    GI = GVCavity1(wi, wh);
    
    G = min(GO, GI);
end

function [G] = GVCavity1(wo, wh)
    denom = dot(wo, wh);
    nom = 2 * wo(3) * wh(3);
    G = min(1, abs(nom/denom));
    G = max(0, G);
end
