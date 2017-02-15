setupSBRot <- function(x,y,z,alpha){
#SETUPSBROT Computes the velocity field and initial condition for the solid
#   body rotation test case from Williamson et. al. (1992) on the UNIT
#   sphere.
#
#   [U,V,W,H] = setupSBRot(X,Y,Z,ALPHA) Computes the velocity field and
#   initial condition for the solid body rotation test case from the
#   standard test suite for the shallow water wave equtions form
#   Williamson et. al. (1992).
#
#   X,Y,Z are represent the Cartesian coordinates for where the velocity
#   fied and initial condition are to be computed.  ALPHA is the the angle
#   of the rotation axis as measured from the north pole.
#
#   U,V,W are the velocity field for solid body rotation expressed with
#   resepct to the Cartesian coordinate system.
#
#   H is the initial condition, which corresponds a compactly supported
#   cosine bell centered at (1,0,0) on the sphere.
#
#   For more details on this test see
#   D. L. Williamson, J. B. Drake, J. J. Hack, R. Jakob, P. N.
#   Swarztrauber, A standard test set for numerical approximations to the
#   shallow water equations in spherical geometry, J. Comput. Phys. 102
#   (1992) 211?224.

# Author: Grady Wright 2014

u <- y*cos(alpha)-z*sin(alpha)
v <- -x*cos(alpha)
w <- x*sin(alpha)

r <- acos(x)

#Cosine Bell
h <- 1/2*(1+cos(3*pi*r))
h(r >= 1/3) <- 0

}
