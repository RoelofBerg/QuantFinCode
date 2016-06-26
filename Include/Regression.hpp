/** \file  Regression.hpp
    \brief Header file declaring classes to perform statistical regression.
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

#ifndef REGRESSION_HPP
#define REGRESSION_HPP
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <boost/function.hpp>

namespace quantfin {

using blitz::Array;
using blitz::Range;
using blitz::firstDim;
using blitz::firstIndex;
using blitz::secondIndex;

class PolynomialLeastSquaresFit {
private:
  Array<double,2> coeff;   ///< Regression coefficients as column vector.
public:
  /// Polynomial regression constructor.
  PolynomialLeastSquaresFit(const Array<double,1>& y,const Array<double,1>& x,int degree = 1);
  double predict(double x) const;
  inline const Array<double,2>& coefficients() const { return coeff; };
};

class LeastSquaresFit {
private:
  Array<double,2> coeff;   ///< Regression coefficients as column vector.
  std::vector<boost::function<double (const Array<double,1>&)> > basis_functions;
public:
  /// Regression constructor.
  LeastSquaresFit(const Array<double,1>& y,
	              const Array<double,2>& x, ///< Independent variables: (number of points) X (number of independent variables)
				  const std::vector<boost::function<double (const Array<double,1>&)> >& xbasis_functions);
  double predict(const Array<double,1>& x) const;
  inline const Array<double,2>& coefficients() const { return coeff; };
};

}

#endif
