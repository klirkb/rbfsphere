computeCyclogenesis <- function(x,y,z,t){
#COMPUTECYCLOGENESIS Computes the solution to the idealized cyclogenesis
#   test case from Nair et. al. (199i) on the UNIT sphere.
#
#   H = computeCyclogenesis(X,Y,Z,T) returns the solution for the idealized
#   cyclogenesis test case from Nair et. al. (199i) on the UNIT sphere at
#   the locations (X,Y,Z) and at time T.
#
#   For more details on this test see
#   R. D. Nair, J. C\^ot\'e, A. Staniforth, Cascade interpolation for
#   semi-Lagrangian advection over the sphere, Quart. J. Roy. Meteor. Soc.
#   125 (1999) 1445?1468.

# Author: Grady Wright 2014

#
# Parameters for this case
#
rho0 <- 3
gamma <- 5

[lam,th] <- cart2sph(x,y,z)

# The angular velocity of the vorticies.
rho <- rho0*cos(th)
w <- (3*sqrt(3)/2*sech(rho).^2.*tanh(rho))./rho
w(abs(rho) < 2*eps) <- 0

# Solution
h <- 1 - tanh(rho/gamma.*sin(lam - w*t))

}
