/** \file  HestonAsset.cpp
    \brief C++ source file implementing class for an asset whose price process follows a stochastic volatility process of the Heston (1993) type.
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

#include <complex>
#include <boost/bind.hpp>
#include "HestonAsset.hpp"
#include "Constants.hpp"

using namespace quantfin; 
using std::complex;

HestonAsset::HestonAsset(CIRprocess& xvol_process,double ini,double xrho,double xlambda,int xn) 
  : vol_process(xvol_process),xzero(ini),rho(xrho),lambda(xlambda),n(xn),gausslaguerre(xn) 
{ 
  kappa = vol_process.get_kappa();
  theta = vol_process.get_theta();
  sigma = vol_process.get_sigma_level();
}

HestonAsset_as_BlackScholesAsset::HestonAsset_as_BlackScholesAsset(const HestonAsset& heston_asset)
  : BlackScholesAsset(NULL,heston_asset.initial_value()),vol_lvl(2)
{
  double rho = heston_asset.get_rho();
  double sgm = heston_asset.get_initial_volatility();
  vol_lvl(0) = rho*sgm;
  vol_lvl(1) = std::sqrt(1.0-rho*rho)*sgm;
  const_vol  = new ConstVol(vol_lvl);
  set_volatility(const_vol);
}

HestonAsset_as_BlackScholesAsset::~HestonAsset_as_BlackScholesAsset()
{
  if (const_vol) delete const_vol;
}

double HestonAsset::option(double mat,double K,double r,int sign) const
{
  double discount = std::exp(-mat*r);
  double fwd      = xzero/discount;
  double x        = std::log(fwd/K);
  boost::function<double (double t)> f0 = boost::bind(&HestonAsset::P_Gatheral,this,_1,0,mat,x,r);
  boost::function<double (double t)> f1 = boost::bind(&HestonAsset::P_Gatheral,this,_1,1,mat,x,r);
  double P1 = gausslaguerre.integrate(f1)/PI + 0.5;
  double P0 = gausslaguerre.integrate(f0)/PI + 0.5;
  double result   = discount*(fwd*P1 - K*P0);
  if (sign==-1) result -= xzero - K*discount;
  return result;
}

double HestonAsset::P_Gatheral(double phi,int j,double mat,double x,double r) const
{
  complex<double> iphi(0.0,phi);
  complex<double> alpha(0.0,j*phi-0.5*phi);
  alpha -= 0.5*phi*phi;
  complex<double> beta   = kappa - j*rho*sigma - iphi*rho*sigma;
  double gamma = 0.5*sigma*sigma;
  complex<double> d      = std::sqrt(beta*beta-4.0*alpha*gamma);
  complex<double> rplus  = (beta+d)/(sigma*sigma);
  complex<double> rminus = (beta-d)/(sigma*sigma);
  complex<double> g      = rminus/rplus;
  complex<double> D      = rminus * (1.0-std::exp(-d*mat))/(1.0-g*std::exp(-d*mat));
  complex<double> C      = kappa * (rminus*mat - std::log((1.0-g*std::exp(-d*mat))/(1.0-g))/gamma);
  complex<double> f      = C*theta + D*vol_process.get_initial();
  complex<double> result = std::exp(iphi*x+f)/iphi;
  return result.real();
}

