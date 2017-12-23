N = 4000;
rL =  20;
hvec = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09];
lx = 0.1;

for i = 1:length(hvec)
    h = hvec(i);
    [f,x,y] = randomsurface(N,rL,h,lx);
end