/** \file  DeterministicVol.cpp
    \brief C++ source file implementing deterministic volatility function classes.
           Copyright 2003, 2006, 2010 by Erik Schlögl

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

#include "DeterministicVol.hpp"
#include "DeterministicVolMediator.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

int DeterministicAssetVol::type() const
{
  return DeterministicVolMediator::FLAT;
}

Array<double,1> DeterministicAssetVol::segments(double t,double dt) const
{
  Array<double,1> result(2);
  result(0) = t;
  result(1) = t + dt;
  return result;              
}

bool DeterministicAssetVol::get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const
{
  return false;   
}

double DeterministicAssetVol::z_integral(int i,double t,double dt) const
{
  throw std::logic_error("z_integral not implemented");
  return 0.0;
}

Array<double,1> DeterministicAssetVol::z_bondintegral(double t,double dt,double T,const DeterministicAssetVol& xmodel) const
{
  throw std::logic_error("z_bondintegral not implemented");
  Array<double,1> result(1);
  return result;
}

Array<double,1> DeterministicAssetVol::z_volintegral(double t,double dt,double T,const DeterministicAssetVol& fxvol) const
{
  throw std::logic_error("z_volintegral not implemented");
  Array<double,1> result(1);
  return result;
}

double DeterministicAssetVol::bondvolexponential(double t,double dt,double T,const Array<double,1>& dZ,const Array<double,1>& dW) const
{
  throw std::logic_error("bondvolexponential not implemented");
  return 0.0;
}

double DeterministicAssetVol::FwdFXexponential(double t,double dt,const Array<double,1>& dWj,const Array<double,1>& dZj0,const Array<double,1>& dZjk,DeterministicAssetVol* xv,DeterministicAssetVol* fxvol) const
{
  throw std::logic_error("FwdFXexponential not implemented");
  return 0.0;       
}

double DeterministicAssetVol::covar(int j,double t,double dt,const DeterministicAssetVol& xv) const
{
  throw std::logic_error("covar() not implemented");
  return 0.0;       
}

boost::shared_ptr<DeterministicAssetVol> DeterministicAssetVolDiff::component_vol(int i) const
{
  boost::shared_ptr<DeterministicAssetVol> result(new DeterministicAssetVolDiff(v1.component_vol(i),v2.component_vol(i)));
  return result;
}

Array<double,1> DeterministicAssetVolDiff::segments(double t,double dt) const
{
  Array<double,1> seg_v1 = v1.segments(t,dt);
  Array<double,1> seg_v2 = v2.segments(t,dt);
  Array<double,1> result = unique_merge(seg_v1,seg_v2);
  return result;              
}

bool DeterministicAssetVolDiff::get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const
{
  Array<double,1> tmp(vol_lvl.extent(firstDim));
  bool result = v1.get_volatility_level(t,T,vol_lvl);
  result = result && v2.get_volatility_level(t,T,tmp);
  vol_lvl -= tmp;
  return result;   
}

double DeterministicAssetVolDiff::volproduct(double t,double dt,const DeterministicAssetVol& xv) const
{
  return DeterministicVolMediator::volproduct_DeterministicAssetVolDiff(t,dt,xv,v1,v2);                                  
}

double ConstVol::volproduct(double t,double dt,const DeterministicAssetVol& xv) const
{
  return DeterministicVolMediator::volproduct_ConstVol(t,dt,lvl,xv);                                  
}

Array<double,1> DeterministicAssetVolDiff::integral(double t,double dt) const
{
  Array<double,1> result(v1.integral(t,dt)-v2.integral(t,dt));   
  return result;
}

Array<double,1> ConstVol::integral(double t,double dt) const
{
  Array<double,1> result(dt*lvl); 
  return result;           
}

int DeterministicAssetVolDiff::factors() const
{
  return v1.factors();    
}

boost::shared_ptr<DeterministicAssetVol> ConstVol::component_vol(int i) const
{
  Array<double,1> clvl(lvl.extent(firstDim));
  clvl = 0.0;
  clvl(i) = lvl(i);
  boost::shared_ptr<DeterministicAssetVol> result(new ConstVol(clvl));
  return result;
}

int ConstVol::factors() const
{
  return lvl.extent(firstDim);    
}

double DeterministicAssetVolDiff::FwdBondVol(double t,double T1,double T2) const
{
  throw std::logic_error("DeterministicAssetVolDiff::FwdBondVol not implemented");
}

/** Integral over the square of the forward zero coupon bond volatility.
  *
  * In the case of a constant volatility vector v, the value of the 
  * integral \f[ \sqrt{\int_t^{T_1}(\sigma^*(u,T_2)-\sigma^*(u,T_1))^2du} \f] is given
  * by \f[ \sqrt{v^2(T_2-T_1)^2(T_1-t)} \f]
  */
double ConstVol::FwdBondVol(double t,double T1,double T2) const
{
  return sqrt(blitz::sum(lvl*lvl) * (T1-T2) * (T1-T2) * (T1-t));
}

/// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
double DeterministicAssetVolDiff::bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const
{
  return v1.bondvolproduct(t,dt,bondmat,xv)-v2.bondvolproduct(t,dt,bondmat,xv);
}

/** Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.

    Here this is 
    \f[ \int_t^{t+\Delta t}-\sigma_B\cdot(T-s)\cdot\sigma_Xds=-\sigma_B\sigma_X\left(T\Delta t-\frac12((t+\Delta t)^2-t^2)\right) \f]
  */
double ConstVol::bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const
{
  return DeterministicVolMediator::bondvolproduct_ConstVol(t,dt,bondmat,lvl,xv);                                  
}

/** Integral over the scalar product between the bond volatility given by (*this) and a 
    bond volatility given by another DeterministicAssetVol object. */
double DeterministicAssetVolDiff::bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const
{
  return v1.bondbondvolproduct(t,dt,T1,T2,xv)-v2.bondbondvolproduct(t,dt,T1,T2,xv);
}

/** Integral over the scalar product between the bond volatility given by (*this) and a 
    bond volatility given by another DeterministicAssetVol object. */
double ConstVol::bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const
{
  return DeterministicVolMediator::bondbondvolproduct_ConstVol(t,dt,T1,T2,lvl,*this,xv);                                  
}

/// Function needed for the exponential-affine representation of zero coupon bond prices.
double DeterministicAssetVolDiff::A(double t,double T) const
{
  throw std::logic_error("DeterministicAssetVolDiff::A() not implemented");       
}

/// Function needed for the exponential-affine representation of zero coupon bond prices.
Array<double,1> DeterministicAssetVolDiff::B(double t,double T) const
{
  throw std::logic_error("DeterministicAssetVolDiff::B() not implemented");       
}

/** Function needed for the exponential-affine representation of zero coupon bond prices.
  * 
  * Zero coupon bond prices can be expressed as an exponential-affine function of the
  * state variables, i.e. as \f[ B(t,T) = \frac{B(0,T)}{B(0,t)}{\mathcal{A}}(t,T)\exp\{-{\mathcal{B}}(t,T)Z(t)\} \f]
  * This is the function \f$ {\mathcal{A}}(t,T) \f$, which for the case of a constant
  * volatility vector v is given by \f[ \exp\left\{-\frac12v^2Tt(T-t)\right\} \f]
  */
double ConstVol::A(double t,double T) const
{
  return exp(-0.5*blitz::sum(lvl*lvl)*T*t*(T-t));
}

/** Function needed for the exponential-affine representation of zero coupon bond prices.
  * 
  * Zero coupon bond prices can be expressed as an exponential-affine function of the
  * state variables, i.e. as \f[ B(t,T) = \frac{B(0,T)}{B(0,t)}{\mathcal{A}}(t,T)\exp\{-{\mathcal{B}}(t,T)Z(t)\} \f]
  * This is the function \f$ {\mathcal{B}}(t,T) \f$, which for the case of a constant
  * volatility vector v is given by \f$ v(T-t). \f$
  */
Array<double,1> ConstVol::B(double t,double T) const
{
  Array<double,1> result(lvl.extent(firstDim));
  result = T-t;
  return result;
}

/// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
double ConstVol::var(int i,double t,double dt) const
{
  return lvl(i)*lvl(i)*dt;       
}

/// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
double ConstVol::covar_dW(int i,double t,double dt) const
{
  return lvl(i)*dt;       
}

/// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
double DeterministicAssetVolDiff::var(int i,double t,double dt) const
{
  throw std::logic_error("DeterministicAssetVolDiff::var() not implemented");
  return 0.0;       
}

/// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
double DeterministicAssetVolDiff::covar_dW(int i,double t,double dt) const
{
  throw std::logic_error("DeterministicAssetVolDiff::covar_dW() not implemented");
  return 0.0;       
}

Array<double,1> DeterministicAssetVolDiff::bondvolintegral(double t,double dt,double T) const
{
  Array<double,1> result(1);
  throw std::logic_error("DeterministicAssetVolDiff::bondvolintegral() not implemented");
  return result;                       
}

Array<double,1> DeterministicAssetVolDiff::StateVariableMean(double t,double dt,double T) const
{
  Array<double,1> result(1);
  throw std::logic_error("DeterministicAssetVolDiff::StateVariableMean() not implemented");
  return result;                       
}

Array<double,1> ConstVol::bondvolintegral(double t,double dt,double T) const
{
  Array<double,1> result(lvl.copy());
  double tdt = t+dt;
  result *= 0.5*(tdt*tdt-t*t) - T*dt;
  return result; 
}

Array<double,1> ConstVol::StateVariableMean(double t,double dt,double T) const
{
  Array<double,1> result(lvl*bondvolintegral(t,dt,T));
  return result; 
}

bool ConstVol::get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const
{
  if (lvl.extent(firstDim)!=vol_lvl.extent(firstDim)) throw std::logic_error("Volatility dimension mismatch");
  vol_lvl = lvl;
  return true;   
}

double ConstVol::z_integral(int i,double t,double dt) const
{
  return lvl(i)*dt;
}

Array<double,1> ConstVol::z_bondintegral(double t,double dt,double T,const DeterministicAssetVol& xmodel) const
{
  Array<double,1> result(lvl*xmodel.bondvolintegral(t,dt,T));
  return result;
}

Array<double,1> ConstVol::z_volintegral(double t,double dt,double T,const DeterministicAssetVol& fxvol) const
{
  Array<double,1> result(lvl*fxvol.integral(t,dt));
  return result;
}

