function sim = gmpairsimilarity(dim, mu1,mu2,sigma1,sigma2)

sigdiv = sigma2\sigma1;
D = 1/2 * ( trace(sigdiv) - dim - log(det(sigma1)/det(sigma2)))...
     + 1/2 * (mu2 - mu1) * inv(sigma2) * (mu2-mu1)';
sim = 1/(1+D);

end