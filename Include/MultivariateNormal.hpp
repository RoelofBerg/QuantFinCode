/** \file  MultivariateNormal.hpp
    \brief Header file declaring a random number generator class for the multivariate normal distribution.
           Copyright 2003, 2006, 2008, 2010, 2011 by Erik Schlögl

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

#ifndef MULTIVARIATENORMAL_HPP
#define MULTIVARIATENORMAL_HPP
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <random/normal.h>
#include "QFQuasiRandom.hpp"

namespace quantfin {

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;
using ranlib::NormalUnit;

/// Random number generator for multivariate normal distribution with zero mean and given covariance matrix.
class MultivariateNormal {
private:
  int                  dim;     ///< Number of variates.
  int               eigdim;     ///< Number of significant eigenvalues.
  NormalUnit<double>&  rnd;     ///< Univariate normal random number generator.
  Array<double,1>   lambda;     ///< Square roots of non-zero eigenvalues of covariance matrix.
  Array<double,1>   eigval;     ///< Non-zero eigenvalues of covariance matrix.
  Array<double,2>   eigvec;     ///< Eigenvectors of the covariance matrix.
  Array<double,2>        v;     ///< Workspace.
  Array<double,2> covariance_matrix;     
  Array<double,2>  inverse;     ///< Moore-Penrose inverse of the covariance matrix.
  double              detR;     ///< Rank R determinant of the covariance matrix.
  void initialise(const Array<double,2>& covar);
  void PDF4CDF(const Array<double,1>& x,Array<double,1>& result,const Array<double,1>& c);
public:
  /// Constructor.
  MultivariateNormal(const Array<double,2>& covar,NormalUnit<double>& xrnd);
  MultivariateNormal(const Array<double,2>& covar);
  inline void transform_random_variables(Array<double,1>& x);
  /// Fill array x with the next set of variates.
  void operator()(Array<double,1>& x);
  /// Evaluate multivariate normal CDF by Monte Carlo integration
  double CDF(const Array<double,1>& d,unsigned long n = 1000000,bool quasi = false);
  void test(const Array<double,2>& covar,int n);
  // Accessors
  inline const Array<double,1>& get_lambda() const { return lambda; }; 
  inline const Array<double,1>& eigenvalues() const { return eigval; };
  inline int rank() const { return eigdim; };
};

inline void MultivariateNormal::transform_random_variables(Array<double,1>& x)
{
  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;
  int c;
  for (c=0;c<eigdim;c++) {
    double tmp = x(c);
    v(c,0) = tmp * lambda(c); }
  Array<double,2> tmp(dim,1);
  tmp = sum(eigvec(i,k)*v(k,j),k);
  x = tmp(blitz::Range::all(),0);
}

}

#endif
