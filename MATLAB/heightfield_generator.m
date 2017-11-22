N = 4000;
rL =  20;
hvec = [0.01 0.02 0.04 0.05];
lx = 0.1;

for i = 1:length(hvec)
    h = hvec(i);
    cd('/Users/mandy/Github/MultiLayerBsdf/MATLAB')
    [f,x,y] = randomsurface(N,rL,h,lx);
end