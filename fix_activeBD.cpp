#include "fix_activeBD.h"

#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_langevin.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixactiveBD::FixactiveBD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{
  if (narg != 8) error->all(FLERR,"Illegal fix active_bd command");


  t_start = utils::numeric(FLERR,arg[3],false,lmp);
  t_target = t_start;
  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  D = utils::numeric(FLERR,arg[5],false,lmp);
  if (D <= 0.0) error->all(FLERR,"Fix bd diffusion coefficient must be > 0.0");
  v_active = utils::numeric(FLERR,arg[6],false,lmp);
  seed = utils::numeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix active bd command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}
/* ---------------------------------------------------------------------- */

FixactiveBD::~FixactiveBD()
{

  delete random;

}

/* ---------------------------------------------------------------------- */


int FixactiveBD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixactiveBD::init()
{
compute_target();


gamma1 = D / force->boltz;
gamma2 = sqrt(2*D);
gamma3 = sqrt(2*3*D);


}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixactiveBD::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}

void FixactiveBD::initial_integrate(int vflag)
{
 
  // update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  // cheating here....
  double *phi = atom->radius;
  double xtmp,ytmp,ztmp;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  double dtinv = 1.0/dt;
  sqrtdt = sqrt(dt);

  // set square root of temperature
  compute_target();

 if (step <= 1) {
    for (int i = 0; i < nlocal; i++)       
     	 
	 phi[i] = 2*3.14159265359*random->uniform();
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      
      phi[i] += sqrtdt * gamma3 * random->gaussian();
      cosphi = cos(phi[i]);
      sinphi = sin(phi[i]);

      xtmp =  v_active*cosphi*dt + dt*gamma1*f[i][0]/t_target + sqrtdt*gamma2*(random->gaussian());
      ytmp =  v_active*sinphi*dt + dt*gamma1*f[i][1]/t_target + sqrtdt*gamma2*(random->gaussian());
      ztmp = dt*gamma1*f[i][2]/t_target + sqrtdt*gamma2*(random->gaussian());

      x[i][0] +=  xtmp;
      x[i][1] +=  ytmp;
      x[i][2] +=  ztmp;

      v[i][0] = xtmp*dtinv;
      v[i][1] = ytmp*dtinv;
      v[i][2] = ztmp*dtinv;

      }
}

/* ---------------------------------------------------------------------- */

