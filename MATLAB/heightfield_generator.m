N = 4000;
rL =  80;
hvec = [0.09];
lx = 0.1;

for i = 1:length(hvec)
    h = hvec(i);
    [f,x,y] = randomsurface(N,rL,h,lx);
end