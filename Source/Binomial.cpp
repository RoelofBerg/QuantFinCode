/** \file  Binomial.cpp
    \brief C++ source file implementing classes for option pricing using binomial trees.
           Copyright 2006 by Erik Schlögl

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

#include "Binomial.hpp"
#include "QFUtil.hpp"

using namespace quantfin;

BinomialLattice::BinomialLattice(const BlackScholesAsset& xS,double xr,double xT,int xN)
  : S(xS),r(xr),T(xT),N(xN),gridslice(xN),tmpslice(xN)
{
  // precompute constants
  dt = T / (N-1);
  discount = std::exp(-r*dt);
  set_CoxRossRubinstein();
}

void BinomialLattice::set_JarrowRudd()
{
  double sgm = S.volatility(0.0,T);
  double mu  = (r - S.dividend_yield(0.0,T) - 0.5*sgm*sgm) * dt;
  double sd  = sgm*std::sqrt(dt);
  u  = std::exp(mu+sd);
  d  = std::exp(mu-sd);
  p  = 0.5;
}

void BinomialLattice::set_CoxRossRubinstein()
{
  double sgm = S.volatility(0.0,T);
  u  = std::exp(sgm*std::sqrt(dt));
  d  = 1.0/u;
  p  = (1.0/discount*std::exp(-S.dividend_yield(0.0,T)*dt)-d)/(u-d);
}

void BinomialLattice::set_Tian()
{
  double sgm = S.volatility(0.0,T);
  double v   = std::exp(sgm*sgm*dt);
  double tmp = std::sqrt(v*v + 2.0*v - 3.0); 
  double R   = std::exp(-S.dividend_yield(0.0,T)*dt)/discount;
  u  = 0.5*v*R * (v + 1.0 + tmp);
  d  = 0.5*v*R * (v + 1.0 - tmp);
  p  = (R-d)/(u-d);
}

/// Binomial lattice method of Leisen/Reimer (1996). K is the target strike.
void BinomialLattice::set_LeisenReimer(double K)
{
  // N-1 is the number of steps. This must be odd.
  int half = N/2;
  N = 2*half;
  gridslice.resize(N);
  tmpslice.resize(N);
  dt = T / (N-1);
  discount = std::exp(-r*dt);
  double sgm    = S.volatility(0.0,T);
  double R      = std::exp(-S.dividend_yield(0.0,T)*dt)/discount;
  double d2     = (std::log(S.initial_value()/K)+(r-S.dividend_yield(0.0,T)-0.5*sgm*sgm)*T)/(sgm*std::sqrt(T));
  double pprime = PeizerPratt(d2+sgm*std::sqrt(T));
  p  = PeizerPratt(d2);
  u  = pprime*R/p;
  d  = (R-p*u) / (1.0-p);
}

double BinomialLattice::PeizerPratt(double z) const
{
  double tmp = z/(double(N-1)+1.0/3.0+0.1/N);
  return 0.5 + sign(z)*std::sqrt(0.25 - 0.25*std::exp(-tmp*tmp*(double(N-1)+1.0/6.0)));
}

std::ostream& quantfin::operator<<(std::ostream& os,const BinLattice& bl)
{
  int j,k;
  int m = bl.dimension();
  for (k=0;k<m;k++) {
    for (j=0;j<k;j++) os << bl(k,j) << ',';
	os << bl(k,j) << std::endl; }
  return os;
}
