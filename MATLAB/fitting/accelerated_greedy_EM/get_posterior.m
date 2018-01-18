function [P,sLw,n] = get_posterior(tree, W)

% P     - (1xk) vector of posteriors of each component k for this box
% sLw   - (1x1) total likelihood of mixture for this box
% n     - (1x1) number of points in this box
%
% Jan Nunnink, 2003

Lw = exp(tree.ll(1,:)).*W';         % for comparison only so keep it simple
sLw = sum(Lw);
if sLw==0
   P=Lw;
else
   P = Lw/sLw;
end
n = tree.numpoints;
