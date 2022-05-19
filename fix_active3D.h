#ifdef FIX_CLASS

FixStyle(active3d,Fixactive3D)

#else

#ifndef LMP_FIX_active3D_H
#define LMP_FIX_active3D_H

#include "fix.h"

namespace LAMMPS_NS {

class Fixactive3D : public Fix {
 public:
  Fixactive3D(class LAMMPS *, int, char **);
  virtual ~Fixactive3D();
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