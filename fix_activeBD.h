#ifdef FIX_CLASS

FixStyle(activebd,FixactiveBD)

#else

#ifndef LMP_FIX_activeBD_H
#define LMP_FIX_activeBD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixactiveBD : public Fix {
 public:
  FixactiveBD(class LAMMPS *, int, char **);
  virtual ~FixactiveBD();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);


 private: 
 double dt, sqrtdt;

 protected:
  class RanMars *random;
  int seed;
  double t_start,t_stop,t_period,t_target,tsqrt;
  double gamma1,gamma2,gamma3;
  double D,cosphi,sinphi;
  double v_active;
  char *id_temp;
  void compute_target();

};

}

#endif
#endif