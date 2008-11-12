#include "trm_subs.h"
#include "trm_vec3.h"
#include "trm_roche.h"

/**
 * Trivial routine to compute the rate of change of radius from the accretor
 * of a particle of given velocity and position.
 * \param r position
 * \param v velocity
 */
double Roche::rdot(const Subs::Vec3& r, const Subs::Vec3& v){
  return Subs::dot(r, v)/r.length();
}


