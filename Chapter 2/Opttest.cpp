/** \file  Opttest.cpp
    \brief C++ source file implementing test suite for optimisation routines.
           Copyright 2005 by Erik Schlögl

		   Redistribution and use in source and binary forms, with or without modification, are permitted provided
		   that the following conditions are met:
		
           -# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

           -# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

           -# The name of the copyright holder may not be used to endorse or promote products derived from this software without specific prior written permission.
		
           THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
		   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
		   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
		   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
		   OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
		   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
		   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
		   POSSIBILITY OF SUCH DAMAGE.
     */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <blitz/array.h>
#include "Powell.hpp"
#include "Linesearch.hpp"

using namespace quantfin;

class testfunc {
public:
  double operator()(const blitz::TinyVector<double,1> x) const;    
};

double testfunc::operator()(const blitz::TinyVector<double,1> x) const
{
  // return -std::exp(-sum(sqr(x)));
  return std::log(1.0+sum(sqr(x+10.0)));
}
      
class testfunc2 {
public:
  double operator()(const blitz::TinyVector<double,2> x) const;    
};

double testfunc2::operator()(const blitz::TinyVector<double,2> x) const
{
  return -std::exp(-sum(sqr(x)));
}
      
class testfunc2A {
public:
  double operator()(const blitz::Array<double,1>& x) const;    
};

double testfunc2A::operator()(const blitz::Array<double,1>& x) const
{
  return -std::exp(-sum(sqr(x)));
}
      
class rosenbrock {
public:
  double operator()(const blitz::TinyVector<double,2> x) const;    
};

double rosenbrock::operator()(const blitz::TinyVector<double,2> x) const
{
  double f0 = 10.0*(x(1)-x(0)*x(0));
  double f1 = 1 - x(0);
  return 0.5 * (f0*f0 + f1*f1);
}
      
class rosenbrockA {
public:
  double operator()(const blitz::Array<double,1>& x) const;    
};

double rosenbrockA::operator()(const blitz::Array<double,1>& x) const
{
  double f0 = 10.0*(x(1)-x(0)*x(0));
  double f1 = 1 - x(0);
  return 0.5 * (f0*f0 + f1*f1);
}
      
int main(int argc, char *argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::flush;
  using blitz::Array;
  using blitz::TinyVector;
  using namespace quantfin;

  size_t i,j,k;
  try {
	/*
    testfunc f;
    TinyVector<double,1> currpos(0.2);
    TinyVector<double,1> direction(1.0);
    opt::BrentLinesearch<testfunc,double,1> ls;
    double result = ls(f,currpos,direction,1.0E-10);
    cout << "One-dimensional argument: ";
    cout << currpos << ' ' << result << endl;
    testfunc2 f2;
    TinyVector<double,2> start2(0.2);
    TinyVector<double,2> currpos2(0.2);
    TinyVector<double,2> direction2(1.0);
    opt::BrentLinesearch<testfunc2,double,2> ls2;
    result = ls2(f2,currpos2,direction2,1.0E-10);
    cout << "Two-dimensional argument: " << endl;
    cout << direction2 << currpos2 << ' ' << result << endl;
    direction2(1) = 0.0;
    currpos2 = start2;
    result = ls2(f2,currpos2,direction2,1.0E-10);
    cout << direction2 << currpos2 << ' ' << result << endl;
    direction2(0) = 0.0;
    direction2(1) = 1.0;
    currpos2 = start2;
    result = ls2(f2,currpos2,direction2,1.0E-10);
    cout << direction2 << currpos2 << ' ' << result << endl;
    opt::Powell<testfunc2,double,opt::BrentLinesearch<testfunc2,double,2>,2> powell(f2,1.0E-15,ls2,100);
    currpos2 = start2;
    cout << "Powell: " << endl;
    result = powell.solve(currpos2);
    cout << currpos2 << ' ' << result << endl;
    currpos2 = start2;
    cout << "Rosenbrock using Powell: " << endl;
    opt::BrentLinesearch<rosenbrock,double,2> rls;
    rosenbrock ros;
    opt::Powell<rosenbrock,double,opt::BrentLinesearch<rosenbrock,double,2>,2> rpowell(ros,1.0E-8,rls,1000);
    result = rpowell.solve(currpos2);
    cout << currpos2 << ' ' << result << endl;
	*/
    cout << "Testing general Array-based version..." << endl;
    testfunc2A f2A;
    Array<double,1> start2A(2);
    start2A = 0.2, 0.2;
    Array<double,1> currpos2A  = start2A.copy();
    Array<double,1> direction2A(2);
    direction2A = 1.0, 1.0;
    opt::GeneralBrentLinesearch<testfunc2A,double> ls2A;
    double result = ls2A(f2A,currpos2A,direction2A,1.0E-10);
    cout << "Two-dimensional argument: " << endl;
    cout << direction2A << currpos2A << ' ' << result << endl;
    direction2A(1) = 0.0;
    currpos2A = start2A;
    result = ls2A(f2A,currpos2A,direction2A,1.0E-10);
    cout << direction2A << currpos2A << ' ' << result << endl;
    direction2A(0) = 0.0;
    direction2A(1) = 1.0;
    currpos2A = start2A;
    result = ls2A(f2A,currpos2A,direction2A,1.0E-10);
    cout << direction2A << currpos2A << ' ' << result << endl;
    opt::GeneralPowell<testfunc2A,double,opt::GeneralBrentLinesearch<testfunc2A,double> > gpowell(f2A,1.0E-15,ls2A,100);
    currpos2A = start2A;
    cout << "General Powell: " << endl;
    result = gpowell.solve(currpos2A);
    cout << currpos2A << ' ' << result << endl;
    currpos2A = start2A;
    cout << "Rosenbrock using general Powell: " << endl;
    opt::GeneralBrentLinesearch<rosenbrockA,double> rlsA;
    rosenbrockA rosA;
    opt::GeneralPowell<rosenbrockA,double,opt::GeneralBrentLinesearch<rosenbrockA,double> > grpowell(rosA,1.0E-8,rlsA,1000);
    result = grpowell.solve(currpos2A);
    cout << currpos2A << ' ' << result << endl;
    
  }  // end try block

  catch (std::logic_error xcpt) {
    cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    cerr << xcpt.what() << endl; }
  return EXIT_SUCCESS;
}
