function plotbox(node)

% plot the bounding box of the points in the node
%
% node     -  the node
%
% plots a rectangle when possible (2 or more points in the node).
% otherwise plots a circle.
% it plots leaves in black and open nodes in blue
%
% Jan Nunnink, 2003

hold on;
if node.numpoints==1
   ellipse(0.3,0.3,0,node.centroid(1),node.centroid(2),'k');
else
   a=node.hyperrect(1,:);
   b=node.hyperrect(2,:);
   c=[a(1) a(1) b(1) b(1);a(1) b(1) b(1) a(1)];
   d=[a(2) b(2) b(2) a(2);b(2) b(2) a(2) a(2)];
   if node.type=='leaf'
      line(c,d, 'color', 'k');
   else
      line(c,d, 'color', 'b');
   end
end
