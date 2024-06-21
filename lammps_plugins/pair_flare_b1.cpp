#include "pair_flare_b1.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <sys/time.h>

// flare++ modules
#include "cutoffs.h"
#include "lammps_descriptor.h"
#include "radial.h"
#include "y_grad.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

typedef unsigned long long timestamp_t;

static timestamp_t
get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}



/* ---------------------------------------------------------------------- */

PairFLAREB1::PairFLAREB1(LAMMPS *lmp) : Pair(lmp) {
  restartinfo = 0;
  manybody_flag = 1;

  beta = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairFLAREB1::~PairFLAREB1() {
  if (copymode)
    return;

  memory->destroy(beta);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairFLAREB1::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype, n_inner, n_count;
  double evdwl, delx, dely, delz, xtmp, ytmp, ztmp, rsq;
  double *coeff;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int beta_init, beta_counter;
  double norm_squared;
  // New vars: beta_p, partial forces
  Eigen::VectorXd single_bond_vals, vals, env_dot, u, beta_p, partial_forces;
  Eigen::MatrixXd single_bond_env_dervs, env_dervs;
  double empty_thresh = 1e-8;

  for (ii = 0; ii < inum; ii++) {
    i = list->ilist[ii];
    itype = type[i];
    jnum = numneigh[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];

    // Count the atoms inside the cutoff.
    n_inner = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      int s = type[j] - 1;
      double cutoff_val = cutoff_matrix(itype-1, s);

      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < (cutoff_val * cutoff_val))
        n_inner++;
    }

    // Compute covariant descriptors.
    double secs;
    single_bond_multiple_cutoffs(x, type, jnum, n_inner, i, xtmp, ytmp, ztmp,
                                 jlist, basis_function, cutoff_function,
                                 n_species, n_max, l_max, radial_hyps,
                                 cutoff_hyps, single_bond_vals,
                                 single_bond_env_dervs, cutoff_matrix);

    // Compute invariant descriptors.
    // Study B1 vs B2
    B1_descriptor(vals, env_dervs, norm_squared, env_dot,
                  single_bond_vals, single_bond_env_dervs, 
                  n_species, n_max, l_max);

    // Continue if the environment is empty.
    if (norm_squared < empty_thresh)
      continue;

    // why do we compute here? Compare to b2 -> Investigaet where power and beta_matrices is pulled from, understand the calculations
    // Compute local energy and partial forces.
    if (power == 1) {
      double B1_norm = pow(norm_squared, 0.5);
      evdwl = beta_matrices[itype - 1].col(0).dot(vals) / B1_norm;
      partial_forces = - env_dervs * beta_matrices[itype - 1].col(0) / B1_norm + evdwl * env_dot / norm_squared;
    } else if (power == 2) {
      beta_p = beta_matrices[itype - 1] * vals;
      evdwl = vals.dot(beta_p) / norm_squared;
      partial_forces =
        2 * (- env_dervs * beta_p + evdwl * env_dot) / norm_squared;
    }

    // Update energy, force and stress arrays.
    n_count = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
  
      if (rsq < (cutoff * cutoff)) {
        double fx = -partial_forces(n_count * 3);
        double fy = -partial_forces(n_count * 3 + 1);
        double fz = -partial_forces(n_count * 3 + 2);
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
  
        if (vflag) {
          ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, delx,
                       dely, delz);
        }
        n_count++;
      }
    }
    // Compute local energy.
    if (eflag)
      ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
  }

  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairFLAREB1::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");

  
  // Set the entire setflag to 1 (otherwise pair.cpp will throw an error)
  // off-diagonal needed for hybrid/overlay (note: we only support pair_coeff * *)
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      setflag[i][j] = 1;
    }
  }

  // Create cutsq array (used in pair.cpp)
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairFLAREB1::settings(int narg, char ** /*arg*/) {
  // "flare" should be the only word after "pair_style" in the input file.
  if (narg > 0)
    error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read DYNAMO funcfl file
------------------------------------------------------------------------- */

void PairFLAREB1::coeff(int narg, char **arg) {
  if (!allocated)
    allocate();

  // Should be exactly 3 arguments following "pair_coeff" in the input file.
  if (narg != 3)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // Ensure I,J args are "* *".
  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  read_file(arg[2]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairFLAREB1::init_style() {
  // Require newton on.
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style requires newton pair on");

  // Request a full neighbor list.
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairFLAREB1::init_one(int i, int j) {
  // init_one is called for each i, j pair in pair.cpp after calling init_style.

  return cutoff;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairFLAREB1::read_file(char *filename) {
  int me = comm->me;
  // diverges cuz no need to have kernel_string here, because B1
  char line[MAXLINE], radial_string[MAXLINE], cutoff_string[MAXLINE], body_order_string[MAXLINE];
  int radial_string_length, cutoff_string_length;
  FILE *fptr;

  // Check that the potential file can be opened.
  if (me == 0) {
    fptr = utils::open_potential(filename,lmp,nullptr);
    if (fptr == NULL) {
      char str[128];
      snprintf(str, 128, "Cannot open potential file %s", filename);
      error->one(FLERR, str);
    }
  }

  int tmp, nwords;
  if (me == 0) {
    fgets(line, MAXLINE, fptr); // Date and contributor

    fgets(line, MAXLINE, fptr); // Power, use integer instead of double for simplicity
    sscanf(line, "%i", &power);

    fgets(line, MAXLINE, fptr); // Body order, B1/2/3
    sscanf(line, "%s", body_order_string);
    if (strcmp(body_order_string, "B1")) {
      error->all(FLERR, "Potential has to be B1");
    }

    fgets(line, MAXLINE, fptr);
    sscanf(line, "%s", radial_string); // Radial basis set
    radial_string_length = strlen(radial_string);

    fgets(line, MAXLINE, fptr);
    sscanf(line, "%i %i %i %i", &n_species, &n_max, &l_max, &beta_size);

    fgets(line, MAXLINE, fptr);
    sscanf(line, "%s", cutoff_string); // Cutoff function
    cutoff_string_length = strlen(cutoff_string);
  }

  MPI_Bcast(&power, 1, MPI_INT, 0, world);
  MPI_Bcast(&n_species, 1, MPI_INT, 0, world);
  MPI_Bcast(&n_max, 1, MPI_INT, 0, world);
  MPI_Bcast(&l_max, 1, MPI_INT, 0, world);
  MPI_Bcast(&beta_size, 1, MPI_INT, 0, world);
  MPI_Bcast(&cutoff, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&radial_string_length, 1, MPI_INT, 0, world);
  MPI_Bcast(&cutoff_string_length, 1, MPI_INT, 0, world);
  MPI_Bcast(radial_string, radial_string_length + 1, MPI_CHAR, 0, world);
  MPI_Bcast(cutoff_string, cutoff_string_length + 1, MPI_CHAR, 0, world);

  // Parse the cutoffs.
  int n_cutoffs = n_species * n_species;
  memory->create(cutoffs, n_cutoffs, "pair:cutoffs");
  if (me == 0)
    grab(fptr, n_cutoffs, cutoffs);
  MPI_Bcast(cutoffs, n_cutoffs, MPI_DOUBLE, 0, world);

  // Fill in the cutoff matrix.
  cutoff = -1;
  cutoff_matrix = Eigen::MatrixXd::Zero(n_species, n_species);
  int cutoff_count = 0;
  for (int i = 0; i < n_species; i++){
    for (int j = 0; j < n_species; j++){
      double cutoff_val = cutoffs[cutoff_count];
      cutoff_matrix(i, j) = cutoff_val;
      if (cutoff_val > cutoff) cutoff = cutoff_val;
      cutoff_count ++;
    }
  }

  // Set number of descriptors.
  int n_radial = n_max * n_species;
  n_descriptors = n_radial;

  // Check the relationship between the power spectrum and beta.
  int beta_check;
  if (power == 1) {
    beta_check = n_descriptors;
  } else if (power == 2) {
    beta_check = n_descriptors * (n_descriptors + 1) / 2;
  } else {
    error->all(FLERR, "Power should be 1 or 2.");
  }
  if (beta_check != beta_size)
    error->all(FLERR, "Beta size doesn't match the number of descriptors.");

  // Set the radial basis.
  if (!strcmp(radial_string, "chebyshev")) {
    basis_function = chebyshev;
    radial_hyps = std::vector<double>{0, cutoff};
  }

  // Set the cutoff function.
  if (!strcmp(cutoff_string, "quadratic"))
    cutoff_function = quadratic_cutoff;
  else if (!strcmp(cutoff_string, "cosine"))
    cutoff_function = cos_cutoff;

  // Parse the beta vectors.
  memory->create(beta, beta_size * n_species, "pair:beta");
  if (me == 0)
    grab(fptr, beta_size * n_species, beta);
  MPI_Bcast(beta, beta_size * n_species, MPI_DOUBLE, 0, world);

  // Fill in the beta matrix.
  // TODO: Remove factor of 2 from beta.
  Eigen::MatrixXd beta_matrix;
  int beta_count = 0;
  double beta_val;

  if (power == 1) {
    for (int k = 0; k < n_species; k++) {
      beta_matrix = Eigen::MatrixXd::Zero(n_descriptors, 1);
      for (int i = 0; i < n_descriptors; i++) {
        beta_matrix(i, 0) = beta[beta_count];
        beta_count++;
      }
      beta_matrices.push_back(beta_matrix);
    }
  } else if (power == 2) {
    for (int k = 0; k < n_species; k++) {
      beta_matrix = Eigen::MatrixXd::Zero(n_descriptors, n_descriptors);
      for (int i = 0; i < n_descriptors; i++) {
        for (int j = i; j < n_descriptors; j++) {
          if (i == j)
            beta_matrix(i, j) = beta[beta_count];
          else if (i != j) {
            beta_val = beta[beta_count] / 2;
            beta_matrix(i, j) = beta_val;
            beta_matrix(j, i) = beta_val;
          }
          beta_count++;
        }
      }
      beta_matrices.push_back(beta_matrix);
    }
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairFLAREB1::grab(FILE *fptr, int n, double *list) {
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line, MAXLINE, fptr);
    ptr = strtok(line, " \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL, " \t\n\r\f")))
      list[i++] = atof(ptr);
  }
}
