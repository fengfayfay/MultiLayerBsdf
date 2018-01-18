function plotmix(X, W, M, R);

% plot data and mixture
%
% plotmix(X, W, M, R)
%
% X - datapoints
% W, M, R - parameters of mixture
%
% Jan Nunnink, 2003

[n,d] = size(X);

if d~=2
   return;
end

% plot data
if isempty(get(gcf,'children')) % no axis defined yet
	plot(X(:,1), X(:,2), 'g+');
   hold on;
else
   h=get(gca,'children');
   delete(h(1:end-1));
end

% plot mixture
for j = 1:length(W)
   S = reshape(R(j,:),d,d);
   [U,L,V] = svd(S); 
   l = diag(L);
   phi = acos(V(1,1));
   if V(2,1) < 0
      phi = 2*pi - phi;
   end
   plot(M(j,1),M(j,2),'r.',M(j,1),M(j,2),'r+');
   ellipse(2*sqrt(l(1)),2*sqrt(l(2)),phi,M(j,1),M(j,2),[W(j)/max(W) 0 0]);
end






