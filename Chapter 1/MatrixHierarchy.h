/** \file  MatrixHierarchy.h
    \brief Declaring a C++ class hierarchy for matrices.
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

/// Abstract base class
class Matrix {
public:
  virtual int rows() const = 0;
  virtual int columns() const = 0;
  virtual double operator()(int i,int j) const = 0;
  virtual double& operator()(int i,int j) = 0;
  /* ... */
};

// Base class must be declared "virtual in order to allow multiple inheritance from classes derived from the class Matrix. 
class SquareMatrix : public virtual Matrix {
public:
  double determinant() const;
  /* ... */
};

// Base class must be declared "virtual in order to allow multiple inheritance from classes derived from the class Matrix. 
class DenseMatrix : public virtual Matrix {
private:
  int     r;    ///< number of rows
  int     c;    ///< number of columns
  double* d;    ///< array of doubles for matrix contents
public:
  DenseMatrix(int nrows,int ncols,double ini = 0.0);
  DenseMatrix(const DenseMatrix& mat);
  DenseMatrix(const Matrix& mat);
  DenseMatrix& operator=(const Matrix& mat);
  inline DenseMatrix& operator=(const DenseMatrix& mat) { return operator=((const Matrix&)mat); };
  virtual ~DenseMatrix();
  virtual int rows() const;
  virtual int columns() const;
  virtual double operator()(int i,int j) const;
  virtual double& operator()(int i,int j);
  /* ... */
};

class DenseSquareMatrix : public DenseMatrix, public SquareMatrix {
public:
  inline DenseSquareMatrix(int nrows,double ini = 0.0) : DenseMatrix(nrows,nrows,ini) { };
  /* ... */
};

class DiagonalMatrix : public SquareMatrix {
private:
  int     r;    ///< number of rows and columns
  double* d;    ///< array of doubles for matrix contents
public:
  DiagonalMatrix(int nrows,double ini = 0.0);
  virtual ~DiagonalMatrix();
  virtual int rows() const;
  virtual int columns() const;
  virtual double operator()(int i,int j) const;
  virtual double& operator()(int i,int j);
  /* ... */
};

DenseMatrix operator+(const Matrix& A,const Matrix& B);
std::ostream& operator<<(std::ostream& os,const Matrix& A);

