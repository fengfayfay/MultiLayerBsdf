%
%    Heigtfield generation
%    for each roughness(alpha) generates two random rough surfaces one is
%    around height 0, another is around height offset
%    
%    N is resolution
%    rL is length of surface
%    hvec is possible height value
%    lx is autocorrelation length
%    hvec/lx is roughness(alpha)
%    dir is saving directory
%    offset is second surface height center
% 


N = 4000;
rL =  20;
hvec = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09];
%hvec = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09];
lx = 0.1;

for i = 1:length(hvec)
    h = hvec(i) + 0.005;
    [f,x,y] = randomsurface(N,rL,h,lx);
end