/** \file  HestonAsset.hpp
    \brief Header file declaring class for an asset whose price process follows a stochastic volatility process of the Heston (1993) type.
           Copyright 2008, 2015 by Erik Schlögl

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

#ifndef HESTONASSET_HPP
#define HESTONASSET_HPP
#include "CIRprocess.hpp"
#include "BlackScholesAsset.hpp"
#include "GaussianQuadrature.hpp"

namespace quantfin {

using blitz::TinyVector;

/// Class for an asset whose price process follows a stochastic volatility process of the Heston (1993) type.
class HestonAsset {
private:
  CIRprocess&      vol_process;     ///< Volatility process.
  double                 xzero;     ///< Initial ("time zero") value.
  double                   rho;     ///< Correlation between Brownian motion driving the volatility process and the Brownian motion driving the asset price process. 
  double                lambda;     ///< Market price of volatility risk parameter.
  double kappa;
  double theta;
  double sigma;
  GaussLaguerre gausslaguerre;
  int n;
  double P_Gatheral(double phi,int j,double mat,double x,double r) const;
public:
  /// Constructor.
  HestonAsset(CIRprocess& xvol_process,     ///< Volatility process.
              double ini,                   ///< Initial ("time zero") value.
              double xrho,                  ///< Correlation between Brownian motion driving the volatility process and the Brownian motion driving the asset price process. 
              double xlambda,               ///< Market price of volatility risk parameter.
			  int xn = 25                   ///< Gauss/Laguerre quadrature refinement
              );
  inline double initial_value() const { return xzero; };  ///< Query the initial ("time zero") value.
  inline void initial_value(double ini) { xzero = ini; }; ///< Set the initial ("time zero") value.
  inline double get_rho() const { return rho; };          ///< Query the correlation between Brownian motion driving the volatility process and the Brownian motion driving the asset price process. 
  inline double get_initial_volatility() const { return std::sqrt(vol_process.get_initial()); };
  /// Price a European call or put option (call is default).
  double option(double mat,double K,double r,int sign = 1) const;
};

/// Class to convert HestonAsset into BlackScholesAsset. Makes the conversion explicit, but avoids the need for BlackScholesAsset to "know" about HestonAsset.
class HestonAsset_as_BlackScholesAsset : public BlackScholesAsset {
private:
  Array<double,1> vol_lvl;
  ConstVol*     const_vol;
public:
  HestonAsset_as_BlackScholesAsset(const HestonAsset& heston_asset);
  virtual ~HestonAsset_as_BlackScholesAsset();
};

}
#endif
