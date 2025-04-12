/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "mpi.h"

// Prototypes
void chem_solver(double, double*, double*, double*,
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // initialize MPI
  int ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
     std::cerr << "Error in calling MPI_Init\n";
     return 1;
  }

  int numprocs, myid;
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // 1. set solver input parameters
  const int maxit = 1000000;
  const double lam = 1.e-2;
  const double eps = 1.e-10;

  double Pbuf[4];
  double Sbuf[3];
  MPI_Status status;

  if (myid == 0) {

    // 2. input the number of intervals
    int n;
    std::cout << "Enter the number of intervals (0 quits):\n";
    std::cin >> n;
    if (n < 1) {
      return 1;
    }

    // 3. allocate temperature and solution arrays
    double *T = new double[n];
    double *u = new double[n];
    double *v = new double[n];
    double *w = new double[n];

    // 4. set random temperature field, initial guesses at chemical densities
    for (int i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
    for (int i=0; i<n; i++)  u[i] = 0.35;
    for (int i=0; i<n; i++)  v[i] = 0.1;
    for (int i=0; i<n; i++)  w[i] = 0.5;

    // 5. start timer
    double stime = MPI_Wtime();

    std::cout << "Starting block 6..." << std::endl;

    // 6. call solver over n intervals
    int numsent = 0;
    for (int i = 1; i < numprocs && numsent < n; i++) {
      Pbuf[0] = T[numsent];
      Pbuf[1] = u[numsent];
      Pbuf[2] = v[numsent];
      Pbuf[3] = w[numsent];
      MPI_Send(Pbuf, 4, MPI_DOUBLE, i, numsent + 1, MPI_COMM_WORLD);
      numsent++;
    }

    int recv = 0;
    while (recv < n) {
      MPI_Recv(Sbuf, 3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      recv++;

      int i = status.MPI_TAG;

      int its = Sbuf[0];
      double res = Sbuf[1];

      if (res < eps) {
        std::cout << "    i = " << i << "  its = " << its << std::endl;
      } else {
        std::cout << "    error: i=" << i << ", its=" << its << ", res=" << res
                  << ", u=" << u[i] << ", v=" << v[i] << ", w=" << w[i] << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      if (numsent < n) {
        Pbuf[0] = T[numsent];
        Pbuf[1] = u[numsent];
        Pbuf[2] = v[numsent];
        Pbuf[3] = w[numsent];
        MPI_Send(Pbuf, 4, MPI_DOUBLE, status.MPI_SOURCE, numsent + 1, MPI_COMM_WORLD);
        numsent++;
      } else {
        MPI_Send(Pbuf, 0, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
      }
    }

    // 7. stop timer
    double ftime = MPI_Wtime();
    double runtime = ftime - stime;

    // 8. output solution time
    std::cout << "     runtime = " << runtime << std::endl;

    // 9. free temperature and solution arrays
    delete[] T;
    delete[] u;
    delete[] v;
    delete[] w;

    // finalize MPI
    ierr = MPI_Finalize();
  } else {
    std::cout << "worker!" << std::endl;

    int its;
    double res;
    int tag;

    while (true) {
      ierr = MPI_Recv(Pbuf, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      tag = status.MPI_TAG;
      if (tag == 0) {
        break;
      }
      chem_solver(Pbuf[0], &(Pbuf[1]), &(Pbuf[2]), &(Pbuf[3]), lam, eps, maxit, &its, &res);
      if (res >= eps) {
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
      }
      Sbuf[0] = its;
      Sbuf[1] = res;
      Sbuf[2] = Pbuf[1];
      ierr = MPI_Send(Sbuf, 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }


  }

  return 0;
} // end main
