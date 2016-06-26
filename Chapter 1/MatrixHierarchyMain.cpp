/** \file  MatrixHierarchyMain.cpp
    \brief Example program demonstrating usage of overloaded
           operator () for subscripting.
           Copyright 2013 by Erik Schlögl

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
#include "MatrixHierarchy.h"

using namespace std;

int main(int argc, char *argv[])
{
    int i,j;
    DenseMatrix A(3,5),B(3,5),C(3,5),D(3,5);
    // assign values to matrix elements
    for (i=0;i<3;i++) {
      for (j=0;j<5;j++) A(i,j) = B(i,j) = C(i,j) = 0.1*i*j; }
	D = A + B + C;
	cout << "Matrix D is: " << D << endl;
	DenseSquareMatrix E(4);
	DiagonalMatrix F(4);
    for (i=0;i<4;i++) {
	  F(i,i) = -1.5;
      for (j=0;j<4;j++) E(i,j) = 1.3*i+j; }
    cout << "The sum of E and F is " <<  E+F << endl;

    return EXIT_SUCCESS;
}
