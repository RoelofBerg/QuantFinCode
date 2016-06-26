/** \file  Regression.cpp
    \brief C++ source file implementing classes to perform statistical regression.
           Copyright 2006, 2010 by Erik Schlögl

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

#include "Regression.hpp"
#include "InterfaceCLAPACK.hpp"

using namespace quantfin;

PolynomialLeastSquaresFit::PolynomialLeastSquaresFit(const Array<double,1>& y,const Array<double,1>& x,int degree)
  : coeff(degree+1,1)
{
  int i,j,k;
  // set up linear equations
  int n = degree+1;
  Array<double,2> A(n,n);
  Array<double,2> B(n,1);
  A = 0.0;
  B = 0.0;
  for (i=0;i<y.extent(firstDim);i++) {
    double xr = 1.0;
    for (j=0;j<n;j++) {
      B(j,0) += xr*y(i);
      int row = 0;
      for (k=j;k>=0;k--) {
        A(row,k) += xr;
        row++; }
      xr *= x(i); }
    for (j=n;j<=2*degree;j++) {
      int m = degree;
      for (k=j-degree;k<n;k++) {
        A(k,m) += xr;
        m--; }
      xr *= x(i); }}
  // solve
  interfaceCLAPACK::SolveLinear(A,coeff,B);
}

double PolynomialLeastSquaresFit::predict(double x) const
{
  int i;
  double result = 0.0;
  for (i=coeff.extent(firstDim)-1;i>=0;i--) result = x*result + coeff(i,0);
  return result;
}

LeastSquaresFit::LeastSquaresFit(const Array<double,1>& y,
								 const Array<double,2>& x,
								 const std::vector<boost::function<double (const Array<double,1>&)> >& xbasis_functions)
  : coeff(xbasis_functions.size(),1),basis_functions(xbasis_functions)
{
  int i,j,k;
  // set up linear equations
  int n = basis_functions.size();
  Array<double,2> A(n,n);
  Array<double,2> B(n,1);
  Array<double,1> xr(n);
  A = 0.0;
  B = 0.0;
  for (i=0;i<y.extent(firstDim);i++) {
	Array<double,1> current_x = x(i,Range::all());
    for (j=0;j<n;j++) {
	  xr(j) = (basis_functions[j])(current_x);
      B(j,0) += xr(j)*y(i);
	  for (k=0;k<=j;k++) A(j,k) += xr(j)*xr(k); }}
  for (j=0;j<n;j++) {
	for (k=0;k<j;k++) A(k,j) = A(j,k); }
  // solve
  interfaceCLAPACK::SolveLinear(A,coeff,B);
}

double LeastSquaresFit::predict(const Array<double,1>& x) const
{
  int i;
  double result = 0.0;
  for (i=coeff.extent(firstDim)-1;i>=0;i--) result += coeff(i,0) * (basis_functions[i])(x);
  return result;
}
