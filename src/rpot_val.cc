#include <cstdlib>
#include <cmath>
#include <iostream>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_vec3.h"
#include "trm_roche.h"

/**
 * rpot_val computes the value of the Roche potential for a specific value of phi & lambda.
 * phi refers to the orbital phase, lambda to a multiplier that specified the position
 * of a point from an origin plus the multiplier time lambda.
 * \param q mass ratio  = M2/M1
 * \param iangle  the orbital inclination, degrees. 90 = edge on.
 * \param p position of origin (units of separation)
 * \param phi phase
 * \param lam multiplier
 * \return Roche potential at point.
 */

double Roche::rpot_val(double q, double iangle, const Subs::Vec3& p, double phi, double lam){

    if(q <= 0.) 
	throw Roche_Error("Roche::rpot_val(double, double, const Subs::Vec3&, double, double): q = " + Subs::str(q) + " <= 0.");
    
    Subs::Vec3 earth;
    set_earth(iangle, phi, earth); 
    
    Subs::Vec3 r = p + lam*earth;
    
    return rpot(q, r);
}
