
alphavec = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
W = [];
M = [];
R = [];
for i = 1:length(alphavec)
alpha = alphavec(i);
gaussiannum = 20;
extend = 0.5;

[obj, W, M, R] = OneLayer3D(alpha, gaussiannum, extend, W, M, R);
end