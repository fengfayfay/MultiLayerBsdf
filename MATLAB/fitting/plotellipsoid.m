figure
hold on
for i = 1:100
    plot_gaussian_ellipsoid(M(i,:), reshape(R(i,:),3,3))
end