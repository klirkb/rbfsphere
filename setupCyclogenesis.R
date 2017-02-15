setupCyclogenesis <- function(x,y,z){
#SETUPCYCLORGENESIS Computes the velocity field and initial condition for
#   the idealized cyclogenesis test case from Nair et. al. (199i) on the UNIT
#   sphere.
#
#   [U,V,W,H] = setupCyclogenesis(X,Y,Z) returns the velocity field and
#   initial condition for the idealized cyclogenesis test case from Nair
#   et. al. (199i) on the UNIT sphere at the locations (X,Y,Z).
#
#   U,V,W are the velocity field for expressed with resepct to the
#   Cartesian coordinate system.  H is the initial condition, which
#   corresponds a compactly supported cosine bell centered at (1,0,0) on
#   the sphere.
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

us <- w.*cos(th)
vs <- 0*us

[!,!,!,u,v,w] <- sphv2cartv(lam,th,us,vs)

# Initial condition
h <- computeCyclogenesis(x,y,z,0)

}
