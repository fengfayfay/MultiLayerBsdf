function [SW, SWX, SWXX, LL, tree] = make_stats(tree, W, M, R, time)

% make_stats - function to compute statistics from 1 box
%
% SW - (kx1)
% SWX - (kxd)
% SWXX - (kxd^2)
% LL - (1x1) average loglikelihood of box under entire mixture
% tree - node
% W, M, R - parameters of mixture
% time - vector with time index for each component
%
% Jan Nunnink, 2003

[k,d] = size(M);
   
% check on cached loglikelihoods
a = size(tree.ll,2);
b = length(time);
if a<b
   tree.ll(:,(a+1):b)=0;        % make size of ll equal to time
end
ff = (time==tree.ll(2,:)); 					% check which cached ll's are old and new
old = find(~ff);
new = find(ff);

% compute loglikelihoods
L(old,:) = opt_approx(tree, M(old,:), R(old,:));		% optimal shared loglikelihood
L(new,:) = tree.ll(1,new)';					% get rest from tree

% update tree
tree.ll(1,old) = L(old)';
tree.ll(2,old) = time(old);

ff = exp(L)'*W;
ff(find(ff<realmin))=realmin;
LL = log(ff)*tree.numpoints;		% approx total log-likelihood of mixture
                                                % = free energy of mixture under optimal P 

L = L - max(L);                     % to get better values. this has no effect on value of P
L = exp(L);

% compute posteriors
L = L.*W;										   
P = L/sum(L);
SW = P*tree.numpoints;
SWX = SW*tree.centroid;

for i=1:k
   tmp = SW(i)*tree.cov;
   tmp = tmp(:)';
   SWXX(i,:) = tmp;
end

