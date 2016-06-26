/** \file  Matrix.h
    \brief Declaring a C++ class for matrices.
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

class Matrix {
private:
  int     r;    ///< number of rows
  int     c;    ///< number of columns
  double* d;    ///< array of doubles for matrix contents
public:
  Matrix(int nrows,int ncols,double ini = 0.0);
  Matrix(const Matrix& mat);
  ~Matrix();
  inline int rows() const { return r; };
  inline int columns() const { return c; };
  inline double operator()(int i,int j) const;
  inline double& operator()(int i,int j);
  friend Matrix operator+(const Matrix& A,const Matrix& B);
  friend Matrix operator+(const Matrix& A,double x);
  friend Matrix operator*(const Matrix& A,const Matrix& B);
  friend Matrix operator*(const Matrix& A,double x);
  /* ... */
};

std::ostream& operator<<(std::ostream& os,const Matrix& A);

inline double Matrix::operator()(int i,int j) const
{
  return d[i*c + j];
}

inline double& Matrix::operator()(int i,int j) 
{
  return d[i*c + j];
}