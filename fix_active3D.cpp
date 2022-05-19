#include "fix_active3D.h"
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math_extra.h"
#include "atom.h"
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


void crossProduct(double v_A[], double v_B[], double c_P[]) {
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

/* ---------------------------------------------------------------------- */

Fixactive3D::Fixactive3D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
  
{
  if (narg != 8) error->all(FLERR,"Illegal fix active_3d command");


  t_start = utils::numeric(FLERR,arg[3],false,lmp);
  t_target = t_start;
  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  D = utils::numeric(FLERR,arg[5],false,lmp);
  if (D <= 0.0) error->all(FLERR,"Fix bd diffusion coefficient must be > 0.0");
  v_active = utils::numeric(FLERR,arg[6],false,lmp);
  seed = utils::numeric(FLERR,arg[7],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix active 3d command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}
/* ---------------------------------------------------------------------- */

Fixactive3D::~Fixactive3D()
{

  delete random;

}

/* ---------------------------------------------------------------------- */


int Fixactive3D::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void Fixactive3D::init()
{
compute_target();


gamma1 = D / force->boltz;
gamma2 = sqrt(2*D);
gamma3 = sqrt(2*3*D);


}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void Fixactive3D::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}

void Fixactive3D::initial_integrate(int vflag)
{
 
  // function to update x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  double mag;


  // Using 'radius' to store particle orentiation angle
  //double *phi = atom->radius;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int step = update->ntimestep;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);

  // set square root of temperature
  compute_target();

  // Set initial particle orientation
 if (step <= 1) {
    for (int i = 0; i < nlocal; i++){       
       // Initialise Active vector with random direction at beginining of simulation
       mu[i][0] = random->gaussian();       
       mu[i][1] = random->gaussian(); 
       mu[i][2] = random->gaussian();

       // Normailise active vector
       mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2)+pow(mu[i][2],2));
       double mag_inv = 1.0 / mag;
       mu[i][0] *= mag_inv;
       mu[i][1] *= mag_inv;
       mu[i][2] *= mag_inv;
    }
  }
  // Integrator 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      
      // Initialise and normalise Gaussian white-noise vector
      double gauss_vector[] = {random->gaussian(),random->gaussian(),random->gaussian()};
      mag = sqrt(pow(gauss_vector[0],2)+pow(gauss_vector[1],2)+pow(gauss_vector[2],2));
      double gauss_mag_inv = 1.0 / mag;
      gauss_vector[0]*=gauss_mag_inv;
      gauss_vector[1]*=gauss_mag_inv;
      gauss_vector[2]*=gauss_mag_inv;

      // Cross product Active vector with Noise
      double c_P[] = {0.0,0.0,0.0};
      crossProduct(mu[i],gauss_vector,c_P);

      // Scale with rotation diffusion coefficient and update active direction
      mu[i][0] += c_P[0]*sqrtdt*gamma3;
      mu[i][1] += c_P[1]*sqrtdt*gamma3;
      mu[i][2] += c_P[2]*sqrtdt*gamma3;

      // Normalise updated Active vector
      mag = sqrt(pow(mu[i][0],2)+pow(mu[i][1],2)+pow(mu[i][2],2));
      double mu_mag_inv = 1.0/mag;
      mu[i][0] *= mu_mag_inv;
      mu[i][1] *= mu_mag_inv;
      mu[i][2] *= mu_mag_inv;

      // Update positions
      x[i][0] +=  v_active*mu[i][0]*dt + dt*gamma1*f[i][0]/t_target + sqrtdt*gamma2*(random->gaussian());
      x[i][1] +=  v_active*mu[i][1]*dt + dt*gamma1*f[i][1]/t_target + sqrtdt*gamma2*(random->gaussian());
      x[i][2] +=  v_active*mu[i][2]*dt + dt*gamma1*f[i][2]/t_target + sqrtdt*gamma2*(random->gaussian());                                                                                 
      }
}

/* ---------------------------------------------------------------------- */
