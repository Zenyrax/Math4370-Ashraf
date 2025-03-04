/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <iostream>
#include "vec1d.hpp"


// tolerance for linear independence checks
#define TOL 1e-12

// Gram-Schmidt process for orthonormalizing a set of vectors
int GramSchmidt1d(vec1d X[], int numvectors) {

  // check that there is work to do
  if (numvectors < 1)  return 0;

  // normalize first vector
  double Xnorm = TwoNorm(X[0]);
  if (std::abs(Xnorm) < TOL) {
    std::cerr << "GramSchmidt error: vectors are linearly-dependent!" << std::endl;
    return 1;
  } else {
    X[0].Scale(1.0/Xnorm);
  }

  // iterate over remaining vectors, performing Gram-Schmidt process
  for (int i=1; i<numvectors; i++) {

    // subtract off portions in directions of existing basis vectors
    for (int j=0; j<i; j++)
      X[i].LinearSum(1.0, X[i], -Dot(X[i],X[j]), X[j]);

    // normalize vector, checking for linear dependence
    Xnorm = TwoNorm(X[i]);
    if (std::abs(Xnorm) < TOL) {
      std::cerr << "GramSchmidt error: vectors are linearly-dependent!" << std::endl;
      return 1;
    } else {
      X[i].Scale(1.0/Xnorm);
    }
  }

  // return success
  return 0;
}
