/** \file  CIRprocess.hpp
    \brief Header file defining class representing a multidimensional process of the type considered by Cox/Ingersoll/Ross (1985).  
           Copyright 2008 by Erik Schlögl

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

#ifndef CIRPROCESS_HPP
#define CIRPROCESS_HPP
#include <stdexcept>
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>

namespace quantfin { 

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;
using blitz::Range;
using blitz::toEnd;
using blitz::firstIndex;
using blitz::secondIndex;

/** Multidimensional process of the type considered by Cox/Ingersoll/Ross (1985):
    \f[ dV_t = \kappa(\theta-V_t)dt + \sigma\sqrt{V_t}dW_t \f]
	*/
class CIRprocess { 
private:
  double                 kappa;  ///< Speed of mean reversion
  double                 theta;  ///< Level of mean reversion
  double                 Vzero;  ///< Initial value of the process
  const Array<double,1>& sigma;  ///< Volatility vector
public:
  inline CIRprocess(double xkappa,  ///< Speed of mean reversion
                    double xtheta,  ///< Level of mean reversion
                    double xVzero,  ///< Initial value of the process
                    const Array<double,1>& xsigma  ///< Volatility vector
					) : kappa(xkappa),theta(xtheta),Vzero(xVzero),sigma(xsigma) { };
  /// Query the number of factors driving the process.
  inline int factors() const { return sigma.extent(firstDim); };
  inline double get_kappa() const { return kappa; };
  inline double get_theta() const { return theta; };
  inline double get_sigma_level() const { return std::sqrt(blitz::sum(sigma*sigma)); };
  inline double get_initial() const { return Vzero; };
};

}
 
#endif

