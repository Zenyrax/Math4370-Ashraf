/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

#ifndef VEC1D_DEFINED__
#define VEC1D_DEFINED__

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>


// This defines a simple arithmetic vector class
class vec1d {

 private:

  ///// Contents /////
  long int length;
  double* data;

 public:

  ///// General class routines /////
  // constructor (initializes values to 0.0)
  vec1d(long int len);

  // destructor
  ~vec1d();

  // write myself to stdout
  void Write() const;

  // write myself to a file
  void Write(const char* outfile) const;

  // returns my overall length
  long int Length() const {return length;};

  // access my data array
  double* GetData() const {return data;};

  // access my data -- access value by location in array,
  // returns handle to value so that it can be changed as well,
  // does not check for legal 'i'
  double& operator[](long int i) {return data[i];};
  double operator[](long int i) const {return data[i];};


  ///// Arithmetic operations defined on a vec1d /////

  // in-place operations (x is the calling vec1d) -- 0=success, 1=fail
  void LinearSum(double a, const vec1d& y,      // x = a*y + b*z
                 double b, const vec1d& z);
  void Scale(double a);                         // x = x*a
  void Copy(const vec1d& y);                    // x = y
  void Constant(double a);                      // x = a

  // scalar quantites derived from vectors
  double Min();                                 // min x_i
  double Max();                                 // max x_i

};  // end vec1d


// independent constructor routines
vec1d Linspace(double a, double b, long int n);
vec1d Random(long int n);

// independent arithmetic routines
double Dot(const vec1d& x, const vec1d& y);     // sum_i (x_i * y_i)
double TwoNorm(const vec1d& x);                 // sqrt(sum_i x_i^2)
double RmsNorm(const vec1d& x);                 // sqrt(sum_i x_i^2 / n)
double MaxNorm(const vec1d& x);                 // max_i |x_i|

#endif
