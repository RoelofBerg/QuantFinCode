/** \file  BlackScholesAsset.hpp
    \brief Header file declaring class for an asset whose price process follows a geometric Brownian motion.
           Copyright 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013 by Erik Schlögl

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

#ifndef BLACKSCHOLESASSET_HPP
#define BLACKSCHOLESASSET_HPP
#include <boost/math/distributions/normal.hpp>
#include "DeterministicVol.hpp"
#include "TermStructure.hpp"

namespace quantfin {

class GenericBlackScholes {
private:
  boost::math::normal N;
public:
  double operator()(double S1,double S2,double vol,int sign) const;
};

/// Class for an asset whose price process follows a geometric Brownian motion.
class BlackScholesAsset {
private:
  DeterministicAssetVol* v;                    ///< Volatility function.
  boost::shared_ptr<DeterministicAssetVol> sv; ///< For the case that the volatility function object needs to be created by the BlackScholesAsset constructor.
  double                 xzero;                ///< Initial ("time zero") value.
  boost::shared_ptr<TermStructure>  dividend;  ///< Term structure of dividend yields.
  boost::math::normal N;
  GenericBlackScholes genericBlackScholes;
public:
  /// Constructor.
  inline BlackScholesAsset(DeterministicAssetVol* xv,  ///< Volatility function.
                           double ini,                 ///< Initial ("time zero") value.
                           double div = 0.0            ///< Dividend yield.
                           ) : v(xv),xzero(ini),dividend(new FlatTermStructure(div,0.0,100.0)) { };
  BlackScholesAsset(const char* path);
  /// Copy constructor.
  inline BlackScholesAsset(const BlackScholesAsset& xasset) : v(xasset.v),sv(xasset.sv),xzero(xasset.xzero),dividend(xasset.dividend->pointer_to_copy()) { };
  inline double initial_value() const { return xzero; };  ///< Query the initial ("time zero") value.
  inline void initial_value(double ini) { xzero = ini; }; ///< Set the initial ("time zero") value.
  /// (Forward) dividend yield between times u and v.
  inline double dividend_yield(double u,double v) const { return dividend->forward_yield(u,v); };  
  /// (Forward) dividend discount factor between times u and v.
  inline double dividend_discount(double u,double v) const { return (*dividend)(v)/(*dividend)(u); };  
  inline void dividend_yield(double ini) { dividend.reset(new FlatTermStructure(ini,0.0,100.0)); }; 
  /// Price a European call or put option (call is default).
  double option(double mat,double K,double r,int sign = 1,double today = 0.0) const;
  /// Delta of a European call or put option (call is default).
  double delta(double mat,double K,double r,int sign = 1) const;
  /// Gamma of a European call or put option (call is default).
  double gamma(double mat,double K,double r,int sign = 1) const;
  /// Option to exchange one asset for another.
  double Margrabe(const BlackScholesAsset& S,double mat,double K, int sign = 1) const;
  double DoleansExp(double t,double T,const Array<double,1>& dW) const;
  inline double volatility(double t,double ttm) const { return std::sqrt(v->volproduct(t,ttm,*v)/ttm); };
  /// Change the volatility function.
  inline void set_volatility(DeterministicAssetVol* xv) { v = xv; };
  /// Covariance of logarithmic asset prices.
  inline double covariance(double t,double ttm,const BlackScholesAsset& S) const { return v->volproduct(t,ttm,*(S.v)); };
  /// Access the volatility function.
  inline const DeterministicAssetVol& volatility_function() const { return *v; };
  /// Access the volatility function.
  inline boost::shared_ptr<DeterministicAssetVol> volatility_function_shared() const { return sv; };
  /// Calculate the implied volatility for a given price.
  double implied_volatility(double price,double mat,double K,double r,int sign = 1) const;
};

}
#endif
