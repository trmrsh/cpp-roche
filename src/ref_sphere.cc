#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "trm_subs.h"
#include "trm_vec3.h"
#include "trm_constants.h"
#include "trm_roche.h"

using Subs::operator+;

/**
 * ref_sphere computes the radius of a reference sphere just touching a Roche distorted star
 * along the line of centres and centred upon its centre of mass. This sphere, which is guaranteed
 * to enclose the equipotential inquestion can then be used to define regions for searching for equipotential
 * crossing when computing eclipses. The Roche-distorted star is defined by the mass ratio
 * and the (linear) filling factor defined as the distance from the centre of mass of the star to its
 * surface in the direction of the L1 point divided by the distance to the L1 point. A filling factor
 * = 1 is Roche filling.
 *
 * \param q    the mass ratio = M2/M1
 * \param ffac linear filling factor. 
 * \param star specifies which star, primary or secondary is under consideration
 * \param rref radius of the reference sphere
 * \param pref reference potential. Roche potential on surface of distorted star.
 */

void Roche::ref_sphere(double q, double ffac, STAR star, double& rref, double& pref){
  if(star == PRIMARY){
    rref = ffac*xl1(q);
    pref = rpot(q, Subs::Vec3(rref,0.,0.));
  }else if(star == SECONDARY){
    rref = ffac*(1.-xl1(q));
    pref = rpot(q, Subs::Vec3(1.-rref,0.,0.));
  } 
}





