/** \file  BlitzScalarFunction.cpp
    \brief C++ source file demonstrating the use of the Blitz++ template library.
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

#include <iostream>
#include <blitz/array.h>

using namespace blitz;

double add_one(double x)
{
  return x+1.0;      
}

double half(int x)
{
  double result = x;
  return x/2.0;      
}

double hypotenuse(double a,double b)
{
  return sqrt(a*a+b*b);      
}

int fit_into(double x,double y)
{
  return y/x;   
}

BZ_DECLARE_FUNCTION(add_one)        
BZ_DECLARE_FUNCTION_RET(half,double)
BZ_DECLARE_FUNCTION2(hypotenuse)
BZ_DECLARE_FUNCTION2_RET(fit_into,int)

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  //using blitz::Array;
  //using blitz::Range;
  
  Array<double,2> A(4,4);
  A = 3.0;
  Array<double,2> B(add_one(A));
  Array<double,2> C(hypotenuse(A,B));
  Array<double,2> D(4,4);
  blitz::firstIndex i;
  blitz::secondIndex j;
  D = 1.0 + j + i*4.0;
  Array<int,2>    E(fit_into(A,D));
  Array<double,2> F(half(E));
  cout << A << endl << B << endl << C << endl << D << endl << E << endl << F << endl;

  return EXIT_SUCCESS;
}
