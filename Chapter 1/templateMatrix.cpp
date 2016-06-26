/** \file  templateMatrix.cpp
    \brief Implementing a hierarchy of C++ templates for matrices.
           Copyright 2005, 2013 by Erik Schlögl

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
#include <complex>

using namespace std;

template<class T> class Matrix {
public:
  virtual ~Matrix();
  virtual T operator()(int i,int j) const = 0;
  virtual T& operator()(int i,int j) = 0;
  /* ... */
};

template<class T> Matrix<T>::~Matrix() { }

template<class T> class DenseMatrix : public Matrix<T> {
private:
  int     r;    ///< number of rows
  int     c;    ///< number of columns
  T*      d;    ///< array of type T for matrix contents
public:
  DenseMatrix(int nrows,int ncols,T ini = 0);
  virtual ~DenseMatrix();
  virtual T operator()(int i,int j) const;
  virtual T& operator()(int i,int j);
  /* ... */
};

template<class T> DenseMatrix<T>::~DenseMatrix()
{
  delete[] d;               
}

template<class T> DenseMatrix<T>::DenseMatrix(int nrows,int ncols,T ini)
{
  r = nrows;               
  c = ncols;
  d = new T[r*c];
  int i;
  for (i=0;i<r*c;i++) d[i] = ini;
}

template<class T> T DenseMatrix<T>::operator()(int i,int j) const
{
  return d[i*c+j];               
}

template<class T> T& DenseMatrix<T>::operator()(int i,int j)
{
  return d[i*c+j];               
}

int main(int argc, char *argv[])
{
    int n = 5;
    int m = 3;
    int i,j;
    DenseMatrix<double> D(n,m);
	// instead of creating a class for complex numbers from scratch, use the complex template from the Standard Template Library
    DenseMatrix<complex<double> > C(n,m);
    // assign values to matrix elements
    for (i=0;i<n;i++) {
      for (j=0;j<m;j++) {
		 D(i,j) = 0.1*i*j; 
		 C(i,j) = complex<double>(i,j); }}
    // access matrix elements
    double sum = 0.0;
    complex<double> csum = 0.0;
    for (i=0;i<n;i++) {
      for (j=0;j<m;j++) {
		sum  += D(i,j); 
		csum += C(i,j); }}
    cout << "The sum of the matrix elements in D is " <<  sum << endl;
    cout << "The sum of the matrix elements in C is " <<  csum << endl;

    return EXIT_SUCCESS;
}
