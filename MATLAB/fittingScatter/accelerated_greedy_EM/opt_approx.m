function L = opt_approx(tree, M, R)

% function to compute the average log-likelihood of a box
%
% L = opt_approx(tree, M, R};
%
% L - (kx1) average loglikelihood for each component
% tree - box
% M, R - parameters of the mixture
%
% Jan Nunnink, 2003

[k, d] = size(M);

for i=1:k
    
   S = reshape(R(i,:),d,d);
   if cond(S)>1e10
      L(i,1)=0;
   else
      
    Si = inv(S);
    
    L(i,1) = ...
        -(d/2)*log(2*pi) -.5*log(det(S)) ...
        -.5*( trace(Si*tree.cov) + M(i,:)*Si*M(i,:)' - 2*M(i,:)*Si*tree.centroid' );
	end    
end
