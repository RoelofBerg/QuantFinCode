/** \file  PiecewiseVol.cpp
    \brief C++ source file implementing deterministic volatility function classes with piecewise constant parameters.
           Copyright 2005, 2006, 2010 by Erik Schl�gl

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

#include "PiecewiseVol.hpp"
#include "QFArrayUtil.hpp"
#include "DeterministicVolMediator.hpp"

using namespace quantfin;

/// Constructor.
PiecewiseConstVol::PiecewiseConstVol(const Array<double,1>& xT, ///< Time line defining calendar time segments on which parameters are constant.
                                     const Array<double,2>& xv  ///< Matrix of volatility scale parameters. Dimensions: time segment (row) by driving factor (column).
                                     ) : timeline(xT),v(xv),vol_sq(xv.extent(firstDim)) 
{ 
  int i,j;
  for (i=0;i<v.extent(firstDim);i++) {
    vol_sq(i) = 0.0;
	for (j=0;j<v.extent(secondDim);j++) vol_sq(i) += v(i,j)*v(i,j); }
}

boost::shared_ptr<DeterministicAssetVol> PiecewiseConstVol::component_vol(int i) const
{
  Array<double,2> cv(v.extent(firstDim),v.extent(secondDim)); // Matrix of volatility scale parameters. Dimensions: time segment (row) by driving factor (column).
  cv = 0.0;
  cv(blitz::Range::all(),i) = v(blitz::Range::all(),i);
  boost::shared_ptr<DeterministicAssetVol> result(new PiecewiseConstVol(timeline,cv));
  return result;
}

bool PiecewiseConstVol::get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const
{
  if (vol_lvl.extent(firstDim)!=v.extent(secondDim)) throw std::logic_error("Volatility dimension mismatch");
  int i = find_segment(t,timeline);
  if (t==timeline(i+1)) i++;
  int j = find_segment(T,timeline);
  if ((i!=j)&&(t!=T)) return false;
  vol_lvl = v(i,blitz::Range::all());
  return true;   
}

/// The integral over the scalar product between two piecewise constant volatility vectors.
double PiecewiseConstVol::volproduct(double t,double dt,const DeterministicAssetVol& xv) const
{
  Array<double,1> lvl(v.extent(secondDim));
  double result  = 0.0;
  int    i       = find_segment(t,timeline);
  double covered = 0.0;
  while ((covered<dt-1e-12)&&(i<v.extent(firstDim))) {
	double tstep = std::min(timeline(i+1)-std::max(t,timeline(i)),dt-covered);
    lvl      = v(i,blitz::Range::all());
    result  += DeterministicVolMediator::volproduct_ConstVol(t+covered,tstep,lvl,xv);                                  
    covered += tstep;
    i++; }
  if (covered<dt-1e-12) throw std::logic_error("Cannot extrapolate");
  return result;
}

Array<double,1> PiecewiseConstVol::integral(double t,double dt) const
{
  Array<double,1> result(v.extent(secondDim));
  result = 0.0; 
  int    i       = find_segment(t,timeline);
  double covered = 0.0;
  while ((covered<dt-1e-12)&&(i<v.extent(firstDim))) {
    double tstep = std::min(timeline(i+1)-std::max(t,timeline(i)),dt-covered);
    result  += tstep * v(i,blitz::Range::all()); 
    covered += tstep;
    i++; }
  if (covered<dt-1e-12) throw std::logic_error("Cannot extrapolate");
  return result;           
}

int PiecewiseConstVol::factors() const
{
  return v.extent(secondDim);
}

/// Integral over the square of the forward zero coupon bond volatility.
double PiecewiseConstVol::FwdBondVol(double t,double T1,double T2) const
{
  throw std::logic_error("Can't do PiecewiseConstVol::FwdBondVol() yet");
}

Array<double,1> PiecewiseConstVol::segments(double t,double dt) const
{
  Array<double,1> tmp(2);
  tmp(0) = t;
  tmp(1) = t + dt;
  int i = find_segment(t,timeline);
  int j = find_segment(t+dt,timeline);
  Array<double,1> tmp2(timeline(blitz::Range(i+1,j)));
  return unique_merge(tmp,tmp2);              
}

/// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
double PiecewiseConstVol::bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const
{
  throw std::logic_error("Can't do PiecewiseConstVol::bondvolproduct() yet");
}

/** Integral over the scalar product between the bond volatility given by (*this) and a 
    bond volatility given by another DeterministicVol object. */
double PiecewiseConstVol::bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const
{
  throw std::logic_error("Can't do PiecewiseConstVol::bondbondvolproduct() yet");
}

/// Function needed for the exponential-affine representation of zero coupon bond prices.
double PiecewiseConstVol::A(double t,double T) const
{
  throw std::logic_error("PiecewiseConstVol::A() not yet implemented");       
}

/// Function needed for the exponential-affine representation of zero coupon bond prices.
Array<double,1> PiecewiseConstVol::B(double t,double T) const
{
  throw std::logic_error("PiecewiseConstVol::B() not yet implemented");       
}

/// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
double PiecewiseConstVol::var(int i,double t,double dt) const
{
  throw std::logic_error("PiecewiseConstVol::var() not implemented");
  return 0.0;       
}

/// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
double PiecewiseConstVol::covar_dW(int i,double t,double dt) const
{
  throw std::logic_error("PiecewiseConstVol::covar_dW() not implemented");
  return 0.0;       
}

Array<double,1> PiecewiseConstVol::bondvolintegral(double t,double dt,double T) const
{
  Array<double,1> result(1);
  throw std::logic_error("PiecewiseConstVol::bondvolintegral() not implemented");
  return result;                       
}

Array<double,1> PiecewiseConstVol::StateVariableMean(double t,double dt,double T) const
{
  Array<double,1> result(1);
  throw std::logic_error("PiecewiseConstVol::StateVariableMean() not implemented");
  return result;                       
}


