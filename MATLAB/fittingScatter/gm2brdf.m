

function gm2brdf(obj, dim, angle, alpha, extendratio, filename, energyRatio, ni, no)
%
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
            %wo = [sintheta*cos(phi(i)), sintheta*sin(phi(i)), mu(j)];
            h1 = signmu* (ni * wi + no * wo);
            h = h1/norm(h1);

            x = h(1)/h(3);
            y = h(2)/h(3);
            
            if dim==2
                p = pdf(obj,[x,y]);
            else
                p = pdf(obj,[x, y,theta]);
                %conditioned on this incident angle
                p = p/(1/(pi/2));
            end
            
            %Jacobian
            detJ = computeJacobian(wi, wo, ni, no);
            %if ni > no 
            %    detJ = detJ * ni * ni;
            %end    
            brdfcos(i,j) = p * detJ * (1 + extendratio*2) * scale;
        end
    end
   
    figure
    imagesc(brdfcos)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['Reflectance from fitted Gaussian,  angle=',num2str(angle), ' alpha=', num2str(alpha)])
    saveas(gcf,[filename,'.jpeg'])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:)));
end
end

function [jacobian] = computeJacobian(wo, wi, no, ni)
    
    denom = (wo(3) + wi(3));
    denom = abs(denom * denom * denom);
    jacobian = abs(wo(3)* wi(3)  + wo(2)*wi(2) + wo(1)*wi(1) + 1)/denom;
     

    %{
    xo = wo(1);
    yo = wo(2);
    zo = wo(3);
    
    xi = wi(1);
    yi = wi(2);
    zi = wi(3);

    hx = ni * xi + no * xo;
    hy = ni * yi + no * yo;
    hz = ni * zi + no * zo;
    
    dxdxi = ni * hz * zi + ni * xi * hx;
    dydyi = ni * hz * zi + ni * yi * hy;
    denom = zi * hz * hz;
    dydxi = ni * hy * xi;
    dxdyi = ni * hx * yi;
    
    jacobian_xiyi = (dxdxi * dydyi - dxdyi * dydxi)/(denom * denom);
   
    %eta = ni/no;
    %dxydwi = (1.0/eta)  
    mu_t = sqrt(zi * zi + ni* ni- no * no)*ni/no; 
    %jacobian = abs(jacobian_xiyi * wi(3) * ni/no);
    jacobian = abs(jacobian_xiyi * wi(3));
    %}
end 
