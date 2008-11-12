#ifndef TRM_ROCHE
#define TRM_ROCHE

//! \file

#include <string>
#include "trm_subs.h"
#include "trm_vec3.h"

//! Roche namespace

/**
 * This namespace wraps up all the functions associated with the roche package. The roche package
 * deals with gas streams, Roche lobes and the like. The accretor
 * is defined as star 1 or the "primary star". The mass ratio is always q = M2/M1. Distances
 * are normalised by the binary separation. The origin is centred on the primary star, with the
 * X axis pointing from the primary to the secondary star.
 */

namespace Roche {

  //! Enum specifying which star in a binary
  enum STAR {PRIMARY, SECONDARY};

  //! Default directory for command defaults
  const char ROCHE_DIR[] = ".roche";

  //! Environment variable for switching directory for command defaults
  const char ROCHE_ENV[] = "ROCHE_ENV";

  //! An exception class.

  /** Roche::Roche_Error is the error class for the Roche programs.
   * It is inherited from the standard string class.
   */
  class Roche_Error : public std::string {
  public:

    //! Default constructor
    Roche_Error() : std::string() {}

    //! Constructor storing a message
    Roche_Error(const std::string& err) : std::string(err) {} 
  };

  //! Computes L1 point distance from primary
  double xl1(double q);
  
  //! Computes L2 point distance from primary
  double xl2(double q);
  
  //! Computes L3 point distance from primary
  double xl3(double q);
  
  //! Computes Roche potential
  double rpot(double q, const Subs::Vec3& p);
  
  //! Computes derivative of Roche potential
  Subs::Vec3 drpot(double q, const Subs::Vec3& p);
  
  //! Calculates primary star's Roche lobe in orbital plane
  void lobe1(double q, float *x, float *y, int n);

  //! Calculates arbitrary Roche equipotential around the Primary in the orbital plane
  void flobe1(double q, float *x, float *y, int n, double pot);
  
  //! Calculates primary star's Roche lobe in orbital plane, velocity Subs::Vec3s
  void vlobe1(double q, float *vx, float *vy, int n);
  
  //! Calculates secondary star's Roche lobe in orbital plane
  void lobe2(double q, float *x, float *y, int n);
  
  //! Calculates secondary star's Roche lobe in orbital plane, velocity coords
  void vlobe2(double q, float *vx, float *vy, int n);
  
  //! Initialises gas stream
  void strinit(double q, Subs::Vec3 &r, Subs::Vec3 &v);
  
  //! Integrates gas stream
  void gsint(double q, Subs::Vec3 &r, Subs::Vec3 &v, double ttry, 
	     double &tdid, double &tnext, double &time, double eps);
  
  //! Works out Jacobi constant
  double jacobi(double q, Subs::Vec3 r, Subs::Vec3 v);
  
  //! Works out whether gas stream hits a given equi-potential around the primary star
  bool hits(double q, double pot, double& x, double& y);

  //! Works out whether gas stream hits a given equi-potential around the primary star
  bool hits(double q, double pot, double& x, double& y, double& vx, double& vy);

  //! Sets the earth vector for a given inclination angle and orbital phase
  void set_earth(const double iangle, const double phase, Subs::Vec3& earth);
 
    //! Works out acceleration in Roche potential
  Subs::Vec3 rocacc(double q, const Subs::Vec3& r, const Subs::Vec3& v);
  
  //! Calculates a gas stream
  void stream(double q, double step, float *x, float *y, int n);

  //! Calculates a gas stream
  void streamr(double q, double rad, float *x, float *y, int n);
  
  //! Calculates a gas stream in velocity coords
  void vtrans(double q, int type, double x, double y, 
	      double vx, double vy, double &tvx, double &tvy);
  
  //! Calculates a gas stream in velocity coords
  void vstream(double q, double step, float *vx, float *vy, 
	       float *rad, int n, int type);
  
  //! Advances a free particle orbit in a Roch potential
  void stradv(double q, Subs::Vec3 &r, Subs::Vec3 &v, double rad, 
	      double acc, double smax);
  
  //! Gas stream at regular radius steps in velocity coords
  void vstrreg(double q, double step, float *vx, float *vy, int n, int type);
  
  //! Determines whether a given point is eclipsed by a Roche-lobe filling secondary star
  bool blink(double q, const Subs::Vec3& r, const Subs::Vec3& e, double acc);

  //! Determines third contact phase
  bool third_contact(double q, double iangle, double rw, int ntheta, double acc, double& contact);

  //! Determines fourth contact phase
  bool fourth_contact(double q, double iangle, double rw, int ntheta, double acc, double& contact);

  //! Compute probability of eclipse by a Roche-lobe filling star
  double eprob(double q);
  
  //! Locate next point at minimum or maximum distanc from the accretor
  void strmnx(double q, Subs::Vec3 &r, Subs::Vec3 &v, double acc);
  
  //! Computes Eggleton's volume-integrated Roche lobe radius
  double rlobe_eggleton(double q);

  //! Computes d log(R/a)/d log(M2) for Eggleton's volume-integrated Roche lobe radius
  double zeta_rlobe_eggleton(double q);

  //! Computes derivative wrt to q of d log(R/a)/d log(M2) for Eggleton's volume-integrated Roche lobe radius
  double dzetadq_rlobe_eggleton(double q);
  
  //! Computes rate of change of radius
  double rdot(const Subs::Vec3& r, const Subs::Vec3& v);
  
  //! Works out whether donor will be irradiated
  bool irrad(double q, double x, double y);
  
  //! Computes data for eclipse by sphere, phase & lambda
  bool sphere_eclipse(double iangle, const Subs::Vec3& r, const Subs::Vec3& c, double rsphere, 
		      double& phi1, double& phi2, double& lam1, double& lam2);

  //! Computes data for eclipse by sphere, lambda only
  bool sphere_eclipse(double iangle, const Subs::Vec3& r, const Subs::Vec3& c, double rsphere, double phase, 
		      double& lam1, double& lam2);

  //! Computes information for eclipse by a disc
  std::vector<std::pair<double,double> > disc_eclipse(double iangle, double rdisc1, double rdisc2, double beta, double height, const Subs::Vec3& r);

  //! Computes radius of reference sphere and Roche potential for Roche-distorted star
  void ref_sphere(double q, double ffac, STAR star, double& rref, double& pref);
  
  //! Computes potential, first and second derivatives wrt coordinates at a point
  void rpot_derivs(double q, const Subs::Vec3& r, double& rpot, Subs::Vec3& fd, Subs::Vec3& sdx, Subs::Vec3& sdy, Subs::Vec3& sdz);

  //! Computes gradients of Roche potential wrt phi and lambda
  void rpot_grad(double q, double iangle, const Subs::Vec3& p, double phi, double lam, double& dphi, double& dlam);

  //! Computes value and gradients of Roche potential wrt phi and lambda
  void rpot_val_grad(double q, double iangle, const Subs::Vec3& p, double phi, double lam, double& rpot, double& dphi, double& dlam);

  //! Computes value of Roche potential as function of phi, lambda
  double rpot_val(double q, double iangle, const Subs::Vec3& p, double phi, double lam);

  //! Computes shadows of Roche lobe in equatorial plane
  void roche_shadow(double q, double iangle, double phi, double dist, double acc, float x[], float y[], bool shade[], int n);

  //! Minmisation of Roche potential
  bool pot_min(double q, double iangle, const Subs::Vec3& p, double phi1, double phi2, double lam1, double lam2, 
	       double rref, double pref, double ftol, double& phi, double& lam);

  //! minimisation along a line accounting for boundaries.
  void linmin(double q, double iangle, const Subs::Vec3& p, double& phi, double& lam, double dphi, double dlam, 
		   double phi1, double phi2, double lam1, double lam2, double pref, double acc, double& pmin, bool& jammed);

  //! Is a point in eclipse or not?
  bool fblink(double q, double iangle, double phi, const Subs::Vec3& p, STAR star, double ffac, double acc);

  //! Computes ingress & egress phases in Roche geometry
  bool ingress_egress(double q, double ffac, double iangle, const Subs::Vec3& r, double& ingress, double& egress, double delta, STAR star);

  //! Computes the position of a point on the Roche distorted surface
  void face(double q, const Subs::Vec3& dirn, double rref, double pref, double acc, STAR star, Subs::Vec3& pvec, Subs::Vec3& dvec, double& r, double& g);

  //! Function object for carrying out line minimisation in phi lambda space
  class Rpot : public Subs::Sfunc {
    
  public:
    Rpot(double q, double iangle, const Subs::Vec3& p, double phi, double dphi, double lam, double dlam) 
      : q_(q), iangle_(iangle), p_(p), phi_(phi), dphi_(dphi), lam_(lam), dlam_(dlam) {}
    
    void reset(double q, double iangle, const Subs::Vec3& p, double phi, double dphi, double lam, double dlam) {
      q_      = q;
      iangle_ = iangle;
      p_      = p;
      phi_    = phi;
      dphi_   = dphi;
      lam_    = lam;
      dlam_   = dlam;
    }

    double operator()(double x) {
      return rpot_val(q_, iangle_, p_, phi_+dphi_*x, lam_+dlam_*x);
    }
    
  private:
    double q_, iangle_;
    Subs::Vec3 p_;
    double phi_, dphi_, lam_, dlam_;
  };
  
  //! Function object for carrying out line minimisation in phi lambda space (gradient)
  class Drpot : public Subs::Sfunc {
    
  public:

    Drpot(double q, double iangle, const Subs::Vec3& p, double phi, double dphi, double lam, double dlam) 
      : q_(q), iangle_(iangle), p_(p), phi_(phi), dphi_(dphi), lam_(lam), dlam_(dlam) {}
    
    void reset(double q, double iangle, const Subs::Vec3& p, double phi, double dphi, double lam, double dlam) {
      q_      = q;
      iangle_ = iangle;
      p_      = p;
      phi_    = phi;
      dphi_   = dphi;
      lam_    = lam;
      dlam_   = dlam;
    }

    double operator()(double x) {
      double dp, dl;
      rpot_grad(q_, iangle_, p_, phi_+dphi_*x, lam_+dlam_*x, dp, dl);
      return (dp*dphi_+dl*dlam_);
    }
    
  private:
    double q_, iangle_;
    Subs::Vec3 p_;
    double phi_, dphi_, lam_, dlam_;
  };


  //! Function object for carrying out line minimisation in lambda space
  class Rlpot : public Subs::Sfunc {
    
  public:
    Rlpot(double q, double iangle, const Subs::Vec3& p, double phi) 
      : q_(q), iangle_(iangle), p_(p), phi_(phi) {}
    
    double operator()(double x) {
      return rpot_val(q_, iangle_, p_, phi_, x);
    }
    
  private:
    double q_, iangle_;
    Subs::Vec3 p_;
    double phi_;
  };
  
  //! Function object for carrying out line minimisation in lambda space (gradient)
  class Dlrpot : public Subs::Sfunc {
    
  public:
    
    Dlrpot(double q, double iangle, const Subs::Vec3& p, double phi) 
      : q_(q), iangle_(iangle), p_(p), phi_(phi) {}
    
    double operator()(double x) {
      double dp, dl;
      rpot_grad(q_, iangle_, p_, phi_, x, dp, dl);
      return dl;
    }
    
  private:
    double q_, iangle_;
    Subs::Vec3 p_;
    double phi_;
  };


  //! Function object to compute Roche potential and its derivative
  /**
   * This is needed for 'rtsafe' inside 'lobe1' and 'flobe1'
   */

  class Lfunc1 : public Subs::RTfunc {

    double dx, dy, qp, cpot;

  public:

    //! Constructor storing fixed data
    Lfunc1(double dxi, double dyi, double qpi, double cpoti) : 
      dx(dxi), dy(dyi), qp(qpi), cpot(cpoti) {}

    //! Function operator
    void operator()(double t, double& f, double& d) const { 
      Subs::Vec3 p, dp;
      p.x() = t*dx;
      p.y() = t*dy;
      p.z() = 0.0;
      f   = rpot(qp, p)-cpot;
      dp  = drpot(qp, p);
      d   = dx*dp.x()+dy*dp.y();
    }

  };

  //! Function object to compute Roche potential and its derivative
  /**
   * This is needed for 'rtsafe' inside 'lobe2'
   */

  class Lfunc2 : public Subs::RTfunc {

    double dx, dy, qp, cpot;

  public:

    //! Constructor storing fixed data
    Lfunc2(double dxi, double dyi, double qpi, double cpoti) : 
      dx(dxi), dy(dyi), qp(qpi), cpot(cpoti) {}

    //! Function operator
    void operator()(double t, double &f, double &d) const { 
      Subs::Vec3 p, dp;
      p.x() = 1.0+t*dx;
      p.y() = t*dy;
      p.z() = 0.0;
      f   = rpot(qp, p)-cpot;
      dp  = drpot(qp, p);
      d   = dx*dp.x()+dy*dp.y();
    }
  };

}
  
#endif



