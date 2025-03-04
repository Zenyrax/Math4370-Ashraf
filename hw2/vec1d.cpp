/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include "vec1d.hpp"


// This file implements the operations defined in the vec1d class


///// General class routines /////

// constructor (initializes values to 0.0)
vec1d::vec1d(long int len) {
  // if len is illegal, create an empty vector
  if (len < 1) {
    this->length = 0;
    this->data = NULL;
  } else {
    this->length = len;
    this->data = new double[length];
    for (long int i=0; i<this->length; i++)
      this->data[i] = 0.0;
  }
}


// destructor (frees space associated with a vec1d)
vec1d::~vec1d() {
  if (this->data!=NULL)
    delete[] this->data;
  this->length = 0;
}


// write myself to stdout
void vec1d::Write() const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Write error, empty data array" );

  // print data to screen
  for (long int i=0; i<this->length; i++)
    std::cout << std::setprecision(17) << "  " << this->data[i] << std::endl;
  std::cout << std::endl;
}


// write myself to a file
void vec1d::Write(const char *outfile) const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Write error, empty data array" );

  // throw exception if 'outfile' is empty
  if (strlen(outfile) < 1)
    throw std::invalid_argument( "vec1d::Write error, empty outfile" );

  // open output file
  std::ofstream fstr;
  fstr.open(outfile);
  if (!fstr.is_open())
    throw std::invalid_argument( "vec1d::Write error, unable to open file for writing" );

  // print data to file
  for (long int i=0; i<this->length; i++)
    fstr << std::setprecision(17) << "  " << this->data[i] << std::endl;
  fstr << std::endl;

  // close output file
  fstr.close();
}



///// Arithmetic operations defined on a given vec1d /////

// x = a*y + b*z
void vec1d::LinearSum(double a, const vec1d& y, double b, const vec1d& z) {
  // check that array sizes match
  if (y.length != this->length  ||  z.length != this->length)
    throw std::invalid_argument( "vec1d::LinearSum error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL || z.data == NULL)
    throw std::invalid_argument( "vec1d::LinearSum error: empty data array" );

  // perform operation
  for (long int i=0; i<this->length; i++)
    this->data[i] = a*y.data[i] + b*z.data[i];
}


//   x = x*a  (scales my data by scalar a)
void vec1d::Scale(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Scale error: empty data array" );

  // perform operation
  for (long int i=0; i<this->length; i++)
    this->data[i] *= a;
}


//   x = y  (copies y into x)
void vec1d::Copy(const vec1d& y) {
  // check that array sizes match
  if (y.length != this->length)
    throw std::invalid_argument( "vec1d::Copy error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL)
    throw std::invalid_argument( "vec1d::Copy error: empty data array" );

  // perform operation
  for (long int i=0; i<this->length; i++)
    this->data[i] = y.data[i];
}


//   x = a  (sets all entries of x to the scalar a)
void vec1d::Constant(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Copy error: empty data array" );

  // perform operation and return
  for (long int i=0; i<this->length; i++)
    this->data[i] = a;
}


///// scalar quantities derived from vectors /////

// min x_i
double vec1d::Min() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Min error: empty data array" );

  // perform operation and return
  double mn = this->data[0];
  for (long int i=0; i<this->length; i++)
    mn = std::min(mn, this->data[i]);
  return mn;
}


// max x_i
double vec1d::Max() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Max error: empty data array" );

  // perform operation and return
  double mx = this->data[0];
  for (long int i=0; i<this->length; i++)
    mx = std::max(mx, this->data[i]);
  return mx;
}



///// independent constructor routines /////

// create a vector of linearly spaced data
vec1d Linspace(double a, double b, long int n) {
  vec1d *x = new vec1d(n);
  double *xd = x->GetData();
  double h = (b-a)/(n-1);
  for (long int i=0; i<n; i++)
    xd[i] = a + h*i;
  return *x;
}


// create a vector of uniformly-distributed random data
vec1d Random(long int n) {
  vec1d *x = new vec1d(n);
  double *xd = x->GetData();
  for (long int i=0; i<n; i++)
    xd[i] = random() / (pow(2.0,31.0) - 1.0);
  return *x;
}


///// independent arithmetic routines /////

// dot-product of x and y
double Dot(const vec1d& x, const vec1d& y) {
  // check that array sizes match
  if (y.Length() != x.Length())
    throw std::invalid_argument( "vec1d::Dot error, vector sizes do not match" );

  // perform operation and return
  double sum = 0.0;
  for (long int i=0; i<x.Length(); i++)
    sum += x[i]*y[i];
  return sum;
}



// ||x||_2
double TwoNorm(const vec1d& x) {
  double sum = 0.0;
  for (long int i=0; i<x.Length(); i++)
    sum += x[i]*x[i];
  return sqrt(sum);
}


// ||x||_RMS
double RmsNorm(const vec1d& x) {
  double sum = 0.0;
  for (long int i=0; i<x.Length(); i++)
    sum += x[i]*x[i];
  return sqrt(sum/x.Length());
}


// ||x||_infty
double MaxNorm(const vec1d& x) {
  double mx = 0.0;
  for (long int i=0; i<x.Length(); i++)
    mx = std::max(mx, std::abs(x[i]));
  return mx;
}
