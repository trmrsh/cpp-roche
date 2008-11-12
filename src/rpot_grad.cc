#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_vec3.h"
#include "trm_roche.h"

/**
 * rpot_grad computes the gradient in phi, lambda space of the Roche potential.
 * phi refers to the orbital phase, lambda to a multiplier that specified the position
 * of a point from an origin plus the multiplier time lambda.
 * \param q mass ratio  = M2/M1
 * \param iangle  the orbital inclination, degrees. 90 = edge on.
 * \param p position of origin (units of separation)
 * \param phi phase
 * \param lam multiplier
 * \param dphi first derivative of Roche potential wrt phi
 * \param dlam first derivative of Roche potential wrt lamda
 */

void Roche::rpot_grad(double q, double iangle, const Subs::Vec3& p, double phi, double lam, double& dphi, double& dlam){

  if(q <= 0.) 
    throw Roche_Error("Roche::rpot_grad(double, double, const Subs::Vec3&, double, double, double& double&, double&): q = " + Subs::str(q) + " <= 0.");

  // Compute cosine and sine of inclination if need be:
  static double iangle_old = -100., sini, cosi;
  if(iangle != iangle_old){
    iangle_old = iangle;
    sini = sin(Constants::PI*iangle/180.);
    cosi = cos(Constants::PI*iangle/180.);
  }


  static double phi_old = -1.e30, sphi, cphi;
  if(phi != phi_old){
    phi_old = phi;
    sphi = sin(Constants::TWOPI*phi);
    cphi = cos(Constants::TWOPI*phi);
  }

  Subs::Vec3 e(  sini*cphi, -sini*sphi, cosi);
  Subs::Vec3 ed(-sini*sphi, -sini*cphi, 0.);

  Subs::Vec3 r = p + lam*e;
  Subs::Vec3 d = Roche::drpot(q, r);

  // derivative wrt phi
  dphi  = Constants::TWOPI*lam*Subs::dot(d, ed);

  // derivative wrt lambda
  dlam  = Subs::dot(d, e);

}
