function Node = buildtree(X, leafbound, plot, stop, opt)

% builds a kd-tree for a data set.
%
% tree = buildtree(X, leafbound, plot, stop, opt)
%
% tree      -  resulting tree
%
% X         -  (nxd) datapoints
% leafbound -  bound on size of bounding box, not used at the moment
% plot      -  plot all leaf boxes or not
% stop      -  criterion to stop building:
%                0 - never stop (build full tree)
%                1 - build only this node and not it's children
%                2 - go only 'opt' levels deep
%                3 - don't split node when it contains less than 'opt' points 
% opt       -  (optional) extra parameters
%
% each node of the tree contains the following fields:
% left      -  (tree) left tree
% right     -  (tree) right tree
% numpoints -  number of points in node
% centroid  -  (1xd) mean of points in node
% cov       -  (dxd) 'covariance' of the points w.r.t. the axes centre
% hyperrect -  (2xd) hyperrectangle bounding the points
% type      -  'leaf' = node is leaf
%              'node' = node has 2 children
% ll        -  (2xk) matrix used to store previously computed average log-likelihoods
%
% Jan Nunnink, 2003

if nargin==4
    opt=0;
end

[n,d] = size(X);
if ~(d==2) plot=0; end

if n==0 fprintf('Error: 0 points in node, causing endless loop, press ctrl-C.\n'); end

Node.ll=zeros(2,0);
Node.type='node'; % default

% if there is 1 datapoint left, make a leaf
if n==1
    Node.left=[];
    Node.right=[];
    Node.centroid=X;
    Node.cov=X'*X;
    Node.hyperrect=[X;X];
    Node.type='leaf';
    Node.numpoints=1;
    if plot==1	plotbox(Node); end
    return;
end

% at least 2 points, calculate the node's easy values
Node.numpoints = n;
a = min(X); b = max(X);
Node.hyperrect = [a; b];

% no box so leaf
if a==b Node.type='leaf'; end

% if we only build this node then stop here and leave children undefined.
if stop==1
    Node.type = 'leaf';
end

% if we only build the tree up to a certain level check if that level is reached.
if stop==2
    if opt==0
        Node.type = 'leaf';
    end
    opt=opt-1;   
end

% check if number of points is below treshold
if stop==3 & n<=opt
    Node.type = 'leaf';
end

% recursively build rest of tree

if Node.type=='leaf'
    Node.centroid = mean(X);
    Node.cov=zeros(d,d);
    Node.cov=(X'*X)/n;
    
    if plot==1 plotbox(Node); end
else
    % split box
    [maxwidth, Node.splitdim] = max(b-a);
    Node.splitval = (a(Node.splitdim)+b(Node.splitdim))/2;
    i=X(:, Node.splitdim)<Node.splitval;
    Node.left = buildtree(X(find(i), :), leafbound, plot, stop, opt);
    Node.right = buildtree(X(find(~i), :), leafbound, plot, stop, opt);
    
    % calc node's values from children's
    Node.centroid = (Node.left.centroid*Node.left.numpoints + ...
                     Node.right.centroid*Node.right.numpoints)/n;
    Node.cov = (Node.left.cov*Node.left.numpoints + ...
                Node.right.cov*Node.right.numpoints)/n;
end



