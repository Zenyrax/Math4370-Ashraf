/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include "vec1d.hpp"

// prototype for Gram-Schmidt routine
int GramSchmidt1d(vec1d X[], int numvectors);


// Example routine to test the vec1d class
int main(int argc, char* argv[]) {

  // create some vecs of length 5, and set some entries
  vec1d a(5), b(5), c(5);
  for (int i=0; i<5; i++)  b[i] = (i+1)*0.1;
  for (int i=0; i<5; i++)  c[i] = (i+1);

  // output to screen
  std::cout << "writing array of zeros" << std::endl;
  a.Write();
  std::cout << "writing array of 0.1, 0.2, 0.3, 0.4, 0.5:" << std::endl;
  b.Write();
  std::cout << "writing array of 1,2,3,4,5:" << std::endl;
  c.Write();

  // verify that b has length 5
  if (b.Length() != 5)
    std::cerr << "error: incorrect vector length" << std::endl;

  // access a's data array, update entries, and write to file
  double *dat = a.GetData();
  dat[0] = 10.0;
  dat[1] = 15.0;
  dat[2] = 20.0;
  dat[3] = 25.0;
  dat[4] = 30.0;
  a.Write("a_data");
  std::cout << "the new file 'a_data' on disk should have entries 10, 15, 20, 25, 30" << std::endl << std::endl;

  // access each entry of a and write to screen
  std::cout << "entries of a, one at a time (via []): should give 10, 15, 20, 25, 30" << std::endl;
  for (int i=0; i<5; i++)
    std::cout << "  " << a[i] << std::endl;
  std::cout << std::endl;

  // update one entry of a
  std::cout << "updating the last entry of a to be 31 (via a[4]):" << std::endl;
  a[4] = 31.0;
  a.Write();
  a[4]= 30.0;  // reset to original

  // Test arithmetic operators
  std::cout << "Testing vector constant, should all be -1" << std::endl;
  b.Constant(-1.0);
  b.Write();

  std::cout << "Testing vector copy, should be 1, 2, 3, 4, 5" << std::endl;
  a.Copy(c);
  a.Write();

  std::cout << "Testing scalar multiply, should be 5, 10, 15, 20, 25" << std::endl;
  c.Scale(5.0);
  c.Write();

  // create a few vecs of length 10
  vec1d X[] = {vec1d(10), vec1d(10), vec1d(10), vec1d(10), vec1d(10)};

  // fill in the vectors
  for (int i=0; i<10; i++) {
    X[0][i] = 1.0*i;
    X[1][i] = -5.0 + 1.0*i;
    X[2][i] = 2.0 + 2.0*i;
    X[3][i] = 20.0 - 1.0*i;
    X[4][i] = -20.0 + 1.0*i;
  }

  // check the LinearSum routine
  X[0].LinearSum(-2.0, X[1], 1.0, X[2]);
  std::cout << "Testing LinearSum, should be all 12:" << std::endl;
  X[0].Write();

  // check the various scalar output routines
  std::cout << "Testing TwoNorm (should be 2.2360679774997898):  " << TwoNorm(b) << std::endl;

  std::cout << "Testing RmsNorm (should be 16.583123951777001):  " << RmsNorm(c) << std::endl;

  std::cout << "Testing MaxNorm (should be 1):  " << MaxNorm(b) << std::endl;

  std::cout << "Testing Min (should be 1):  " << a.Min() << std::endl;

  std::cout << "Testing Max (should be 25):  " << c.Max() << std::endl;

  std::cout << "Testing Dot (should be 275):  " << Dot(a, c) << std::endl;

  std::cout << "Testing Linspace, should be 0 1 2 3 4" << std::endl;
  vec1d d = Linspace(0.0, 4.0, 5);
  d.Write();

  std::cout << "Testing Random" << std::endl;
  vec1d f = Random(5);
  f.Write();


  /// performance/validity tests (Gram-Schmidt process)
  int n=1000000;
  std::cout << "Running GramSchmidt1d process" << std::endl;
  vec1d Y[] = {Random(n), Random(n), Random(n), Random(n), Random(n)};
  std::chrono::time_point<std::chrono::system_clock> stime = std::chrono::system_clock::now();
  if (GramSchmidt1d(Y,5))
    std::cerr << "GramSchmidt1d returned error" << std::endl;
  std::chrono::time_point<std::chrono::system_clock> ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> rtime = ftime-stime;
  std::cout << "  GramSchmidt1d time: " << rtime.count() << std::endl << std::endl;

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

  return 0;
} // end main
