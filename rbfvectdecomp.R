rbfvectdecomp <- function(x,ux,xi,ep){
#RBFVECDECOMP Computes the Helmholtz decomposition of a vector field
#   tangent to the sphere from scattered samples using vector RBF
#   interpolation.
#
#   [UDF,UCF,SF,VP] = rbfvecdecomp(X,UX,XI,EP) returns an approximation to
#   the divergence-free (UDF) and curl-free (UCP) parts of the field UX
#   sampled at the (scattered) nodes X. These values are returned at the
#   node locations specified in XI. The field UDF and UCF are the
#   components of the Helmholtz decomposition of UX.  The field UDF+UCF
#   interpolates the field UX at the XI. SF is a stream function for the
#   UDF field and VP is a velocity potential for the UCF field.
#
#   X is assumed to be an array of size N-by-2 or N-by-3.  In the first
#   case the columns of X are assumed to correspond to the longitude and
#   latitude coordinates of the sample points (LAM,TH), where LAM is the
#   azimuthal angle and TH the elevation (measured from the equator), both
#   in radians.  If X is N-by-3 then the columns of X are assumed to
#   correspond to the Cartesian coordinates of the sample points (X,Y,Z).
#
#   XI should be arranged similarly to X and contain the point locations
#   where the fields should be evaluated.
#
#   UX can either be expressed with respect to spherical or Cartesian
#   coordinates.  In the former case it should be an N-by-2 array, in the
#   latter a N-by-3 array, with column elements corresponding to the column
#   elements of X.
#
#   The code is based on the vector interpolation/decomposition method
#   described in
#      E.J. Fuselier and G.B. Wright. Stability and error estimates for
#      vector field interpolation and decomposition on the sphere with
#      RBFs. SIAM J. Numer. Anal., 47 (2009), 3213-3239.
#   Presently, the code uses the Gaussian RBF in the construction and
#   requires the user input a shape parameter EP.
#
#   The present implementation is a bit memory intensive.

#   Author: Grady Wright 2014.

nd <- size(x,1)

# Determine whether the inputs are with respect to Cartesian or spherical
# coordinates.
if (size(x,2) == 3){
    [lam,th] <- cart2sphm(x)
} else {
    lam <- x(:,1)
    th <- x(:,2)
    x <- [x ones(nd,1)]
    [x(:,1),x(:,2),x(:,3)] <- sph2cart(lam,th,1)
}

if (size(ux,2) == 3){
    t <- zeros(nd,2)
    [!,!,t(:,1),t(:,2)] <- cartv2sphv(x(:,1),x(:,2),x(:,3),ux(:,1),ux(:,2),ux(:,3))
} else {
    t <- ux
}

f <- [-sin(lam) cos(lam) zeros(size(lam))]
e <- [-cos(lam).*sin(th) -sin(lam).*sin(th) cos(th)]

d <- [t(:,2)';t(:,1)']
d <- d(:)

# Distance matrix and basis vectors:
rd2 <- zeros(nd); fjfk <- rd2; fjek <- rd2; ejek <- rd2; fjxk <- rd2; ejxk <- rd2
for (j in 1:3){
   [xd1,xd2] <- meshgrid(x(:,j))
   rd2 <- rd2 + (xd1-xd2).^2
   [fk1,fj1] <- meshgrid(f(:,j))
   fjfk <- fjfk + fj1.*fk1
   [ek1,fj1] <- meshgrid(e(:,j),f(:,j))
   fjek <- fjek + fj1.*ek1
   [ek1,ej1] <- meshgrid(e(:,j))
   ejek <- ejek + ej1.*ek1
   [xk1,fj1] <- meshgrid(x(:,j),f(:,j))
   fjxk <- fjxk + fj1.*xk1
   [xk1,ej1] <- meshgrid(x(:,j),e(:,j))
   ejxk <- ejxk + ej1.*xk1
}
FM <- rbf(rd2,ep,2)
GM <- rbf(rd2,ep,3)
fjfk <- FM.*fjfk
fjek <- FM.*fjek
ejek <- FM.*ejek
clear FM

B1df <- reshape(shiftdim(reshape([-fjfk;fjek],nd,nd,2),1),2*nd,nd)'){
B2df <- reshape(shiftdim(reshape([fjek';-ejek],nd,nd,2),1),2*nd,nd)'){
B1cf <- reshape(shiftdim(reshape([-ejek;-fjek'],nd,nd,2),1),2*nd,nd)'){
B2cf <- reshape(shiftdim(reshape([-fjek;-fjfk],nd,nd,2),1),2*nd,nd)'){
clear fjfk fjek ejek

Bdf <- reshape(shiftdim(reshape([B1df B2df],2*nd,nd,2),2),2*nd,2*nd)){
Bcf <- reshape(shiftdim(reshape([B1cf B2cf],2*nd,nd,2),2),2*nd,2*nd)){
clear B1df B2df B1cf B2cf

t1 <- GM.*fjxk.*fjxk'
t2 <- -GM.*fjxk.*ejxk'
t3 <-  GM.*ejxk.*ejxk'
clear GM fjxk ejxk

C1df <- reshape(shiftdim(reshape([t1;t2],nd,nd,2),1),2*nd,nd)'){
C2df <- reshape(shiftdim(reshape([t2';t3],nd,nd,2),1),2*nd,nd)'){
C1cf <- reshape(shiftdim(reshape([t3;-t2'],nd,nd,2),1),2*nd,nd)'){
C2cf <- reshape(shiftdim(reshape([-t2;t1],nd,nd,2),1),2*nd,nd)'){
clear t1 t2 t3
Cdf <- reshape(shiftdim(reshape([C1df C2df],2*nd,nd,2),2),2*nd,2*nd)){
Ccf <- reshape(shiftdim(reshape([C1cf C2cf],2*nd,nd,2),2),2*nd,2*nd)){
clear C1df C2df C1cf C2cf

Adf <- Bdf + Cdf
Acf <- Bcf + Ccf
clear Bdf Cdf Bcf Ccf

condest(Adf+Acf)
pause

c <- (Adf + Acf)\d
clear Adf Acf

alpha <- c(1:2:length(c));  beta <- c(2:2:length(c))
clear c d

st <- -repmat(beta,1,3).*e + repmat(alpha,1,3).*f
s <- repmat(alpha,1,3).*e + repmat(beta,1,3).*f

clear alpha beta

#
# Evaluation.
#

nv <- size(xi,1)

# Determine whether the inputs are with respect to Cartesian or spherical
# coordinates.
if (size(xi,2) == 3){
    [li,ti] <- cart2sphm(xi)
    cartcoords <- true
} else {
    li <- xi(:,1)
    ti <- xi(:,2)
    xi <- [xi ones(nv,1)]
    [xi(:,1),xi(:,2),xi(:,3)] <- sph2cart(li,ti,1)
    cartcoords <- false
}

fv <- [-sin(li) cos(li) zeros(size(li))]
ev <- [-cos(li).*sin(ti) -sin(li).*sin(ti) cos(ti)]

udf <- zeros(nv,2)
ucf <- zeros(nv,2)
sf <- zeros(nv,1)
vp <- zeros(nv,1)

for (j in 1:nv){
   # Compute the square of the distances between the node points and evaluation
   # points.
   rv2 <- (2 - 2*dot(repmat(xi(j,:),nd,1),x,2))'
   F <- rbf(rv2,ep,2)
   G <- rbf(rv2,ep,3)
   udf(j,2) <- -F*dot(repmat(fv(j,:),nd,1),st,2) + G*(dot(repmat(fv(j,:),nd,1),x,2).*dot(repmat(xi(j,:),nd,1),st,2))
   udf(j,1) <-  F*dot(repmat(ev(j,:),nd,1),st,2) - G*(dot(repmat(ev(j,:),nd,1),x,2).*dot(repmat(xi(j,:),nd,1),st,2))
   ucf(j,2) <- -F*dot(repmat(ev(j,:),nd,1),s,2) + G*(dot(repmat(ev(j,:),nd,1),x,2).*dot(repmat(xi(j,:),nd,1),s,2))
   ucf(j,1) <- -F*dot(repmat(fv(j,:),nd,1),s,2) + G*(dot(repmat(fv(j,:),nd,1),x,2).*dot(repmat(xi(j,:),nd,1),s,2))
   vp(j,1) <- -F*dot(repmat(xi(j,:),nd,1)-x,(s-x.*repmat(dot(s,x,2),1,3)),2)
   sf(j,1) <- F*dot(repmat(xi(j,:),nd,1)-x,cross(x,s,2),2)
}

# Return the vector field with respect to Cartesian coordinates
if (cartcoords){
#     [~,~,~,udf(:,1),udf(:,2),udf(:,3)] = sphv2cartv(li,ti,udf(:,1),udf(:,2));
    udf <- repmat(udf(:,2),1,3).*ev + repmat(udf(:,1),1,3).*fv
    ucf <- repmat(ucf(:,2),1,3).*ev + repmat(ucf(:,1),1,3).*fv
}


# Function for computing the derivatives of the Gaussian kernel that are
# used in the code above.
rbf <- function(r2,ep,k){

switch k
   case 1
      phi <- exp(-ep^2*r2)
   case 2
      phi <- -2*exp(-ep^2*r2)
   case 3
      phi <- 4*ep^2*exp(-ep^2*r2)
   otherwise
      error('not implemented')
}

}

}
