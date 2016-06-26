/** \file  Matrix.cpp
    \brief Implementing a C++ class for matrices.
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

#include "Matrix.h"

Matrix::Matrix(int nrows,int ncols,double ini)
{
  int i;
  r = nrows;
  c = ncols;
  d = new double [nrows*ncols];
  double* p = d;
  for (i=0;i<nrows*ncols;i++) *p++ = ini;
}

Matrix::Matrix(const Matrix& mat)
{
  int i;
  r = mat.r;
  c = mat.c;
  d = new double [r*c];
  double* p  = d;
  double* pm = mat.d;
  for (i=0;i<r*c;i++) *p++ = *pm++;
}

Matrix::~Matrix()
{
  delete[] d;
}

Matrix operator+(const Matrix& A,const Matrix& B)
{
  int i;
  Matrix result(A);
  double* p  = result.d;
  double* pB = B.d;
  for (i=0;i<A.r*A.c;i++) *p++ += *pB++;
  return result;
}

std::ostream& operator<<(std::ostream& os,const Matrix& A)
{
  int i,j;
  for (i=0;i<A.rows();i++) {
	for (j=0;j<A.columns()-1;j++) os << A(i,j) << ',';
	os << A(i,j) << std::endl; }
  return os;
}
