/** \file  GramCharlier.hpp
    \brief Header file declaring class for a distribution given by a Gram/Charlier expansion.
           Copyright 2004, 2005, 2007, 2011 by Erik Schlögl

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

#ifndef GRAMCHARLIER_HPP
#define GRAMCHARLIER_HPP
#include <stdexcept>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include "Polynomial.hpp"
#include "Constants.hpp"
#include <boost/function.hpp>
#include "LogisticMap.hpp"
#include "QFArrayUtil.hpp"
#include "Powell.hpp"

namespace quantfin { 

using blitz::Array;
using blitz::firstDim;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

/// Class for a distribution given by a Gram/Charlier expansion.
class GramCharlier {
private:
  Array<double,1> coeff;     ///< Coefficients of the Gram/Charlier expansion - coefficients beyond length of array are deemed to be zero.
  Polynomial       poly;     ///< Polynomial part of the density.
  LogisticMap       map;     ///< Logistic mapping to enforce bounds for sigma;
  class kz {
  public:
    //inline kz() : double_f_double(std::sqrt(3.0)) { };
    inline kz() { };
    virtual double operator()(double x) const; };
  void initialise_poly();
  /// Set bounds for lambda in line search.
  void set_bounds(Array<double,1>& currpos,const Array<double,1>& direction,double& lambda_min,double& lambda_max);
  class GCcalibration_class {
  private:
    GramCharlier&                                parent;
	boost::function<double (double)> objective_function;
  public:
    inline GCcalibration_class(GramCharlier& xparent,boost::function<double (double)> xobjective_function) 
      : parent(xparent),objective_function(xobjective_function) { };
    double operator()(const Array<double,1>& xpar);
  };
public:
  /// Default constructor (standard normal distribution).
  inline GramCharlier() : coeff(1),map(1e-6,10.0) { coeff = 1.0; };
  /// Constructor.
  GramCharlier(const Array<double,1>& xcoeff     ///< Coefficients of the Gram/Charlier expansion - coefficients beyond length of array are deemed to be zero.
               );
  inline const Array<double,1>& coefficients() const { return coeff; };
  void set_coefficient(int i,double ci);
  inline double get_coefficient(int i) const;
  /// Density.
  inline double operator()(double x) const { return RSQRT2PI * std::exp(-x*x/2.0) * real(poly(x)); };
  /// Query skewness.
  double skewness() const;
  /// Query excess kurtosis.
  double excess_kurtosis() const;
  /// Verify whether Gram/Charlier expansion is a valid density using Jondeau/Rockinger (2001).
  bool JRverify() const;
  /// Verify whether Gram/Charlier expansion is a valid density by checking polynomial roots.
  bool verify();
  /// For a given excess kurtosis, the skewness on the frontier defining the set of valid Gram/Charlier densities as described in Jondeau/Rockinger (2001).
  double JRfrontier_skewness(double k) const;
  /// Transform unconstrained variables into constrained skewness and excess kurtosis, yielding a valid Gram/Charlier density as described in Jondeau/Rockinger (2001).
  void JRmap(double& s,double& k) const;
  void logistic_map(double& x,double a,double b) const;
  inline int highest_moment() const { return coeff.extent(firstDim)-1; };
  /// Find GC coefficients to minimise a given objective function.
  double fit_coefficients(boost::function<double (double)> objective_function,double sgm,int dim,double eps = 1e-12);
};

inline double GramCharlier::get_coefficient(int i) const
{
  if (i<coeff.extent(firstDim)) return coeff(i);
  else return 0.0;
}


}

#endif
