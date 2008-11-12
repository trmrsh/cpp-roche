#include <cmath>
#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_vec3.h"
#include "trm_roche.h"

/**
 * rpot computes the Roche potential at a given point.
 * \param q mass ratio = M2/M1
 * \param p the point in question (units scaled by separation)
 * \return the Roche potential.
 */

double Roche::rpot(double q, const Subs::Vec3& p){

  double mu, comp, x2y2, z2, r1, r2, r1sq;

  if(q <= 0.) throw Roche_Error("q = " + Subs::str(q) + "(<= 0.) in rpot");

  mu   = q/(1+q);
  comp = 1.-mu;
  x2y2 = Subs::sqr(p.x()) + Subs::sqr(p.y());
  z2   = Subs::sqr(p.z());
  r1   = sqrt(r1sq = x2y2+z2);
  r2   = sqrt(r1sq + 1. - 2.*p.x());
  return (-comp/r1-mu/r2-(x2y2+mu*(mu-2.*p.x()))/2);
}
