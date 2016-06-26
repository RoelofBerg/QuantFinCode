/** \file  NewtonRaphsonExample.cpp
    \brief C++ source file demonstrating the use of the NewtonRaphson template.
           Copyright 2005 by Erik Schloegl 
     */

#include <iostream>
#include "NewtonRaphson.hpp"

using namespace quantfin;

// The function
class my_func {
private:
  int dim;
public:
  inline my_func(int xdim) : dim(xdim) { };
  Array<double,1> operator()(Array<double,1>& x);
  Array<double,2> Jacobian(Array<double,1> x);    
  inline int argdim() { return dim; }; 
  inline int retdim() { return dim; }; 
};

Array<double,1> my_func::operator()(Array<double,1>& x)
{
  int i,j;
  Array<double,1> result(retdim());
  result = 0.0;
  for (i=0;i<retdim();i++) {                
    for (j=0;j<argdim();j++) result(i) += pow(x(j),i+1); }
  return result;                
}

Array<double,2> my_func::Jacobian(Array<double,1> x)
{
  int i,j;
  Array<double,2> result(retdim(),argdim());
  for (i=0;i<retdim();i++) {                
    for (j=0;j<argdim();j++) result(i,j) = (i+1) * pow(x(j),i); }
  return result;                
}

int main(int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  
  try {
  // Instantiate functor
  my_func f(3);
  // Instantiate Rootsearch object
  opt::NewtonRaphson<my_func> nr(f,1E-12);
  // Target value
  Array<double,1> target(f.argdim());
  target = 1.0;
  // Starting point
  Array<double,1> x(f.argdim());
  firstIndex i;
  x = 0.525 * i;
  cout << "The solution with starting point\n" << x << "\n is:\n" << nr.solve(x,target) << endl;
  cout << "The initial Jacobian was: " << f.Jacobian(x) << endl;
  x = 0.525 * (3-i);
  cout << "The solution with starting point\n" << x << "\n is:\n" << nr.solve(x,target) << endl;
  cout << "The initial Jacobian was: " << f.Jacobian(x) << endl;
  // Test with numerical Jacobian
  opt::NumericalJacobian<my_func> nj(f);
  // Instantiate Rootsearch object
  opt::NewtonRaphson<opt::NumericalJacobian<my_func> > nrnj(nj,1E-12);
  cout << "The solution with starting point\n" << x << "\n is:\n" << nrnj.solve(x,target) << endl;
  cout << "The initial Jacobian was: " << nj.Jacobian(x) << endl;
  }

  catch (std::logic_error xcpt) {
    cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    cerr << xcpt.what() << endl; }
  return EXIT_SUCCESS;
}
