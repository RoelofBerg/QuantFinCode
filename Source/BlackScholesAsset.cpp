/** \file  BlackScholesAsset.cpp
    \brief C++ source file implementing class for an asset whose price process follows a geometric Brownian motion.
           Copyright 2003, 2005, 2006, 2007, 2009, 2011, 2012, 2013 by Erik Schlögl

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

#include <stdexcept>
#include <functional>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "BlackScholesAsset.hpp"
#include "Rootsearch.hpp"

using namespace quantfin;

double GenericBlackScholes::operator()(double S1,double S2,double vol,int sign) const
{
  double sd  = std::sqrt(vol);
  double h1  = (std::log(S1/S2) + 0.5*vol) / sd;
  return sign * (S1*boost::math::cdf(N,sign*h1)-S2*boost::math::cdf(N,sign*(h1-sd)));
}   

/// Price a European call or put option (call is default).
double BlackScholesAsset::option(double mat,                ///< Time to maturity of the option.
                                 double K,                  ///< Exercise price.
                                 double r,                  ///< Deterministic interest rate.
                                 int sign,                  ///< Call = 1 (default), put = -1.
								 double today
                                 ) const
{
  double vol = v->volproduct(today,mat,*v);
  return genericBlackScholes(xzero*(*dividend)(today+mat)/(*dividend)(today),K*exp(-mat*r),vol,sign);
}

/// Delta of a European call or put option (call is default).
double BlackScholesAsset::delta(double mat,                ///< Maturity of the option.
                                double K,                  ///< Exercise price.
                                double r,                  ///< Deterministic interest rate.
                                int sign                   ///< Call = 1 (default), put = -1.
                                ) const
{
  double vol = v->volproduct(0.0,mat,*v);
  double sd  = std::sqrt(vol);
  double h1  = (std::log((xzero*(*dividend)(mat)/(*dividend)(0.0))/(K*exp(-mat*r))) + 0.5*vol) / sd;
  return sign*(*dividend)(mat)*boost::math::cdf(N,sign*h1);
}

/// Gamma of a European call or put option (call is default).
double BlackScholesAsset::gamma(double mat,                ///< Maturity of the option.
                                double K,                  ///< Exercise price.
                                double r,                  ///< Deterministic interest rate.
                                int sign                   ///< Call = 1 (default), put = -1.
                                ) const
{
  double vol = v->volproduct(0.0,mat,*v);
  double sd  = std::sqrt(vol);
  double h1  = (std::log((xzero*(*dividend)(mat)/(*dividend)(0.0))/(K*exp(-mat*r))) + 0.5*vol) / sd;
  return (*dividend)(mat)*boost::math::pdf(N,h1)/(xzero*sd);
}

double BlackScholesAsset::Margrabe(const BlackScholesAsset& S,double mat,double K, int sign) const
{
  DeterministicAssetVolDiff voldiff(*v,*(S.v));
  double vol = voldiff.volproduct(0.0,mat,voldiff);      
  return genericBlackScholes(xzero*(*dividend)(mat)/(*dividend)(0.0),K*S.xzero*exp(-mat*S.dividend_yield(0.0,mat)),vol,sign);
}

double BlackScholesAsset::DoleansExp(double t,double T,const Array<double,1>& dW) const
{
  Array<double,1> vol_lvl(dW.extent(firstDim));
  if (!v->get_volatility_level(t,T,vol_lvl)) throw std::logic_error("Volatility not constant in BlackScholesAsset::DoleansExp");
  return std::exp(blitz::sum(dW*vol_lvl)-0.5*v->volproduct(t,T-t,*v));
}

/// Calculate the implied volatility for a given price.
double BlackScholesAsset::implied_volatility(double price,double mat,double K,double r,int sign) const
{
  boost::function<double (double)> objective_function = boost::bind(&GenericBlackScholes::operator(),&genericBlackScholes,(double)(xzero*(*dividend)(mat)/(*dividend)(0.0)),(double)(K*exp(-mat*r)),_1,sign);
  Rootsearch<boost::function<double (double)>,double,double> rs(objective_function,price,0.3,0.2,1e-12);
  double vol = rs.solve();
  return sqrt(vol/mat);
}
