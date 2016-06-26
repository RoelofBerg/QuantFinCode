/** \file  BlitzLinearAlgebra.cpp
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
#include "InterfaceCLAPACK.hpp"

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  using blitz::Array;
  using blitz::Range;
  
  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  Array<double,1> x(4),y(4);
  x = i;
  y = 4.0 - i;
  cout << x << endl << y << endl;
  double dot_product = sum(x*y);
  cout << "Dot product: " << dot_product << endl;
  Array<double,2> A(4,4),B(2,4),D(4,4);
  A = x(i)*y(j);
  B = 1.0;
  cout << A << endl << B << endl;
  Array<double,2> C(B.extent(blitz::firstDim),A.extent(blitz::secondDim));
  C = sum(B(i,k)*A(k,j),k);
  cout << "Matrix multiplication: " << C << endl;
  
  A = 1.0,0.0,2.0,-1.0,
      1.0,3.0,0.0,1.0,
      0.0,-1.0,2.0,2.0,
      0.5,8.0,1.0,0.0;
  Array<double,2> X(4,1),R(4,1);
  X = 1.0;
  R = sum(A(i,k)*X(k,j),k);    
  cout << A << endl << X << endl << R << endl;
  Array<double,2> calcX(4,1);
  quantfin::interfaceCLAPACK::SolveLinear(A,calcX,R);
  cout << calcX << endl;

  cout << "Test tridiagonal solver:" << endl;
  D = 1.0, 2.0, 0.0, 0.0,
	  3.0, 4.0, 5.0, 0.0,
	  0.0, 6.0, 7.0, 8.0,
	  0.0, 0.0, 9.0, 10.0;
  R = sum(D(i,k)*X(k,j),k);    
  quantfin::interfaceCLAPACK::SolveTridiagonal(D,calcX,R);
  cout << calcX << endl;

  return EXIT_SUCCESS;
}
