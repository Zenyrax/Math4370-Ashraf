/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>

// Prototypes
void vector_sum(int, int, int, double*, double*, double*);
void vector_scale(int, int, int, double, double*, double*);
double vector_infnorm(int, int, int, double*);
double vector_dotprod(int, int, int, double*, double*);


/* Allocates and fills the multi-dimensional vectors x, y, z
   with random numbers.
   Calls a set of vector routines.
   Writes the total time taken for all vector routines. */
int main(int argc, char* argv[]) {

  // set problem parameters
  const int l = 100;   // length of first index
  const int m = 100;   // length of second index
  const int n = 100;   // length of third index

  // allocate 3D data arrays
  int i, j, k;
  double *x = new double[l*m*n];
  double *y = new double[l*m*n];
  double *z = new double[l*m*n];

  // fill in vectors x and y
  for (i=0; i<l*m*n; i++)
	x[i] = random() / (pow(2.0,31.0) - 1.0);
  for (i=0; i<l*m*n; i++)
	y[i] = random() / (pow(2.0,31.0) - 1.0);

  // start timer
  clock_t stime = clock();

  // call the vector routines a number of times
  double a;
  for (i=0; i<100; i++) {
    vector_sum(l, m, n, x, y, z);
    vector_scale(l, m, n, 2.0, x, z);
    a = vector_infnorm(l, m, n, z);
    a = vector_dotprod(l, m, n, x, y);
  }

  // stop timer
  clock_t ftime = clock();
  double runtime = ((double) (ftime - stime))/CLOCKS_PER_SEC;

  // output solution time
  std::cout << " Result from computation = " << a << std::endl;
  std::cout << " Total run time = " << runtime << std::endl;

  // free data arrays
  delete[] x;
  delete[] y;
  delete[] z;

} // end main
