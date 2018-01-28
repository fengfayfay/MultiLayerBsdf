% reciprocity test
dim = 3;
numg = 50;
reflect = false;
alpha = 0.5;

% 3d gm data
dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
filename = [dir,'3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numg),'_reflect_',num2str(reflect),'.mat'];
load(filename,'obj')

testnum = 1e3;
wi_theta = pi/2*rand(testnum,1);
wo_phi = 2 * pi * rand(testnum,1);
wo_mu = rand(testnum,1);
brdf = zeros(testnum,1);
brdf1 = zeros(testnum,1);
tol = 1e-10;
for i = 1:testnum
    %% wi as incident
    wi = [sin(wi_theta(i)); 0; cos(wi_theta(i))];
    
    wo_sintheta = sqrt(1-wo_mu(i)^2);
    wo = [wo_sintheta*cos(wo_phi(i)); wo_sintheta*sin(wo_phi(i)); wo_mu(i)];
    
    h1 = (wi + wo)/2;
    h = h1/norm(h1);
    
    p = pdf(obj,[h(1)/h(3),h(2)/h(3),wi_theta(i)]);
    % conditioned on this incident angle
    p = p/(1/(pi/2));
    
    % Jacobian
    detJ = (wo(3) * wi(3) + wo(2)*wi(2) + wo(1)*wi(1) + 1)/(wi(3)+wo(3))^3;
    brdf(i)= (reflect+1)*p*detJ / wo_mu(i);
    
    %% wo as incident
    % rotate wi and wo
    wo_rot = rotz(180/pi*-wo_phi(i)) * wo;
    wi_rot = rotz(180/pi*-wo_phi(i)) * wi;
    wo_rot = wo_rot/norm(wo_rot);
    wi_rot = wi_rot/norm(wi_rot);
    
    % check rotation
    cond1 = abs(wo(3)-wo_rot(3))<tol&&abs(wi(3)-wi_rot(3))<tol;
    cond2 = abs(wo_rot(2))<tol;
    assert(cond1,'incorrect rotation, z changes')
    assert(cond2,'incorrect rotation, y not zero')
    
    h1_new = (wi_rot + wo_rot)/2;
    h_new = h1/norm(h1_new);
    
    % check h preserves
    assert(abs(h(1)-h_new(1))<tol && abs(h(2)-h_new(2))<tol && abs(h(3)-h_new(3))<tol,'h doesn not preserve')
    
    wo_theta = acos(wo_rot(3));
    p1 = pdf(obj,[h_new(1)/h_new(3),h_new(2)/h(3),wo_theta]);
    p = p/(1/(pi/2));
    
    % Jacobian
    detJ1 = (wo_rot(3) * wi_rot(3) + wo_rot(2)*wi_rot(2) + wo_rot(1)*wi_rot(1) + 1)/(wi_rot(3)+wo_rot(3))^3;
    
    % check Jacobian preserves
    assert(abs(detJ-detJ1)<tol,'Jacobian doesn not preserve')
    
    brdf1(i)= (reflect+1)*p1*detJ1 / wi_rot(3);
end
figure
plot(abs(brdf-brdf1))