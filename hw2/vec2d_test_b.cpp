/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include "vec2d_b.hpp"

// prototype for Gram-Schmidt routine
int GramSchmidt2d(vec2d X[], int numvectors);

// Example routine to test the vec1d class
int main(int argc, char* argv[]) {

  // create some vecs of size 5x5, and set some entries
  vec2d a(5, 5), b(5, 5), c(5, 5);
  for (int i=0; i<5; i++) for (int j=0; j<5; j++) b(i, j) = (i*5+j+1)*0.1;
  for (int i=0; i<5; i++) for (int j=0; j<5; j++) c(i, j) = (i*5+j+1);

  // output to screen
  std::cout << "writing 2d array of zeros" << std::endl;
  a.Write();
  std::cout << "writing 2d array of 0.1, 0.2, ..., 2.4, 2.5:" << std::endl;
  b.Write();
  std::cout << "writing 2d array of 1,2,...,24,25:" << std::endl;
  c.Write();

  // verify that b has a size of 25
  if (b.Size() != 25)
    std::cerr << "error: incorrect vector length" << std::endl;

  // access a's data array, update entries, and write to file
  double *dat = a.GetData();
  dat[0] = 10.0;
  dat[6] = 15.0;
  dat[12] = 20.0;
  dat[18] = 25.0;
  dat[24] = 30.0;
  a.Write("a_data");
  std::cout << "the new file 'a_data' on disk should have entries 10, 15, 20, 25, 30 on the diagonal and zeroes elsewhere" << std::endl << std::endl;

  // access each entry of a and write to screen
  std::cout << "entries of a, one at a time (via []): should give 10, 15, 20, 25, 30" << std::endl;
  for (int i=0; i<5; i++)
    std::cout << "  " << a(i, i) << std::endl;
  std::cout << std::endl;

  // update one entry of a
  std::cout << "updating the last entry of a to be 31 (via a[4]):" << std::endl;
  a(4,4) = 31.0;
  a.Write();
  a(4,4) = 30.0;  // reset to original

  // Test arithmetic operators
  std::cout << "Testing vector constant, should all be -1" << std::endl;
  b.Constant(-1.0);
  b.Write();

  std::cout << "Testing vector copy, should be 1, 2, ..., 24, 25" << std::endl;
  a.Copy(c);
  a.Write();

  std::cout << "Testing scalar multiply, should be 5, 10, ..., 120, 125" << std::endl;
  c.Scale(5.0);
  c.Write();

  // create a few vecs of size 10x10
  vec2d X[] = {vec2d(10, 10), vec2d(10, 10), vec2d(10, 10), vec2d(10, 10), vec2d(10, 10)};

  // fill in the vectors
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      X[0](i,j) = 1.0*(i + j);
      X[1](i,j) = -5.0 + 1.0*(i + j);
      X[2](i,j) = 2.0 + 2.0*(i + j);
      X[3](i,j) = 20.0 - 1.0*(i + j);
      X[4](i,j) = -20.0 + 1.0*(i + j);
    }
  }

  // check the LinearSum routine
  X[0].LinearSum(-2.0, X[1], 1.0, X[2]);
  std::cout << "Testing LinearSum, should be all 12:" << std::endl;
  X[0].Write();

  // check the various scalar output routines
  std::cout << "Testing TwoNorm:  " << TwoNorm(b) << std::endl;

  std::cout << "Testing RmsNorm:  " << RmsNorm(c) << std::endl;

  std::cout << "Testing MaxNorm (should be 1):  " << MaxNorm(b) << std::endl;

  std::cout << "Testing Min (should be 1):  " << a.Min() << std::endl;

  std::cout << "Testing Max (should be 125):  " << c.Max() << std::endl;

  std::cout << "Testing Dot:  " << Dot(a, c) << std::endl;
  std::cout << "Testing Linspace, should be 0 1 ... 23 24" << std::endl;
  vec2d d = Linspace(0.0, 24.0, 5, 5);
  d.Write();

  std::cout << "Testing Random" << std::endl;
  vec2d f = Random(5, 5);
  f.Write();

  /// performance/validity tests (Gram-Schmidt process)
  std::pair<long int, long int> tests[] = {{10000, 1000}, {1000, 10000}, {100, 100000}, {10, 1000000}, {100000, 100}, {1000000, 10}};
  std::cout << "Running GramSchmidt2d processes" << std::endl;
  for(size_t q = 0; q < 6; q++) {
    long int m = tests[q].first;
    long int n = tests[q].second;
    vec2d Y[] = {Random(m, n), Random(m, n), Random(m, n), Random(m, n), Random(m, n)};
    std::chrono::time_point<std::chrono::system_clock> stime = std::chrono::system_clock::now();
    if (GramSchmidt2d(Y,5))
      std::cerr << "GramSchmidt2d returned error" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> ftime = std::chrono::system_clock::now();
    std::chrono::duration<double> rtime = ftime-stime;
    std::cout << "  GramSchmidt2d time (m =" << m << ",n=" << n << "): " << rtime.count() << std::endl << std::endl;

    std::cout << "Resulting vectors should be orthonormal:" << std::endl;
    bool pass = true;
    double tolerance = 1e-12;
    for (int i=0; i<5; i++) {
      if (std::abs(Dot(Y[i], Y[i]) - 1.0) > tolerance) {
        pass = false;
        std::cout << "  <Y[" << i << "],Y[" << i << "]> = " << Dot(Y[i], Y[i]) << std::endl;
      }
      for (int j=i+1; j<5; j++)
        if (std::abs(Dot(Y[i], Y[j])) > tolerance) {
          pass = false;
          std::cout << "  <Y[" << i << "],Y[" << j << "]> = " << Dot(Y[i], Y[j]) << std::endl;
        }
    }
    if (pass)
      std::cout << "  passed orthonormality check" << std::endl << std::endl;
    else
      std::cout << "  failed orthonormality check" << std::endl << std::endl;

  }
  return 0;
} // end main
