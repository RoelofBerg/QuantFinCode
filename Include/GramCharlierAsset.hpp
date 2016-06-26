/** \file  GramCharlierAsset.hpp
    \brief Header file declaring class for an asset whose price risk neutral distribution is given a Gram/Charlier expansion.
           Copyright 2007, 2009, 2011 by Erik Schlögl

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

#ifndef GRAMCHARLIERASSET_HPP
#define GRAMCHARLIERASSET_HPP
#include <boost/math/distributions/normal.hpp>
#include "GramCharlier.hpp"

namespace quantfin {

/// Class for an asset whose price risk neutral distribution is given a Gram/Charlier expansion.
class GramCharlierAsset {
private:
  GramCharlier&          gc;        ///< Gram/Charlier expanded density for the standardised risk-neutral distribution.
  double                 sgm;       ///< Volatility level.
  double                 xzero;     ///< Initial ("time zero") value.
  double                 T;         ///< Maturity for which the risk-neutral distribution is valid.
  boost::math::normal       N;
  double genericBlackScholes(double S1,double S2,double vol,int sign) const;
  /// Objective function for best fit calibration to a given set of Black/Scholes implied volatilities.
  double calibration_objective_function(double xsgm);
  // Variables used during calibration.
  const Array<double,1>* strikes_;
  Array<double,1>        prices_;
  double                 domestic_;
  double                 foreign_;
public:
  /// Constructor.
  inline GramCharlierAsset(GramCharlier& xgc,          ///< Gram/Charlier expanded density for the standardised risk-neutral distribution.
	                       double xsgm,                ///< Volatility level.
                           double ini,                 ///< Initial ("time zero") value.
                           double xT                   ///< Maturity for which the risk-neutral distribution is valid.
                           ) : gc(xgc),sgm(xsgm),xzero(ini),T(xT),prices_(3) { };
  inline double initial_value() const { return xzero; };  ///< Query the initial ("time zero") value.
  inline void initial_value(double ini) { xzero = ini; }; ///< Set the initial ("time zero") value.
  inline double maturity() const { return T; };
  /// Price a European call option.
  double call(double K,double domestic_discount,double foreign_discount = 1.0) const;
  /// Best fit calibration to a given set of Black/Scholes implied volatilities.
  double calibrate(const Array<double,1>* xstrikes,const Array<double,1>* xvols,double domestic_discount,double foreign_discount,int highest_moment);
  inline double standard_deviation() const { return sgm; };
  inline double skewness() const { return gc.skewness(); };
  inline double excess_kurtosis() const { return gc.excess_kurtosis(); };
};

}
#endif
