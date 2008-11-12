#include <cstdlib>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_vec3.h"
#include "trm_roche.h"

/** 
 * set_earth sets the earth vector given an inclination and orbital phase. For speed
 * it stores the cosine and sine of the various angles and checks on next entry whether
 * they need changing. 
 * \param  iangle orbital inclination
 * \param  phase  orbital phase (1 = one orbit)
 * \param  earth  the earth vector
 */

void Roche::set_earth(const double iangle, const double phase, Subs::Vec3& earth){

  static double iangle_old = -1.e30, sini, cosi;
  if(iangle != iangle_old){
    iangle_old = iangle;
    cosi       = cos(Constants::TWOPI*iangle/360.);
    sini       = sin(Constants::TWOPI*iangle/360.);
  }

  static double phase_old = -1.e30, sinp, cosp;
  if(phase != phase_old){
    phase_old = phase;
    cosp      = cos(Constants::TWOPI*phase);
    sinp      = sin(Constants::TWOPI*phase);
  }

  earth.set(sini*cosp, -sini*sinp, cosi);
}


