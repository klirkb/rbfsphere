plotSphNodes <- function(x){
#PLOTSPHNODES Makes a simple plot of a given set of nodes on the surface of
#   the unit sphere.
#
#   plotSphNodes(X) plots the nodes contained in X on the surface of the
#   unit sphere.  Here X is assumed to be an N-by-3 array with each row
#   consisting of the (x,y,z) Cartesian coordinates for a node.
#
#   cax = plotSphNodes(X) returns a handle to Axes for the plot.
#
#   Example:
#       x = 2*rand(101,3)-1;
#       x = bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));  # Project to the sphere.
#       plotSphNodes(x);
#

# Author: Grady Wright, 2014

# Make sure the nodes are on the unit sphere.
#x <- bsxfun(@rdivide,x,sqrt(sum(x.^2,2)))
x <- apply(x,2,function(x) x/sqrt(apply(apply(x,1,function(x) x^2),1,sum)))
#
# Generate a unit sphere.
#
#[xx,yy,zz] <- sphere(101)
# Color of the sphere will be yellow:
#clr <- [255 255 102]/255
# Plot the sphere
#cax <- newplot
#hold on
persp(xx,yy,zz) #,1+0*xx,'EdgeColor','None','FaceColor',clr)

#
# Add the nodes as black solid circles.
#
cloud(x[,3]~x[,1]*x[,2])
#plot3(x(:,1),x(:,2),x(:,3),'k.','MarkerSize',6)
#hold off
#daspect([1 1 1])
#view(3)

#if (nargout == 1){
#    varargout <- cax
#}

}



