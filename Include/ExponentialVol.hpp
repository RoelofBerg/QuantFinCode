/** \file  ExponentialVol.hpp
    \brief Header file declaring exponentially decaying deterministic volatility function classes.
           Copyright 2006, 2010, 2011 by Erik Schlögl

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

#ifndef EXPONENTIALVOL_HPP
#define EXPONENTIALVOL_HPP

#include "DeterministicVol.hpp"

namespace quantfin {

/* This class represents the exponentially decaying volatility \f$ \sigma e^{-a t} \f$. */
class ExponentialVol : public DeterministicAssetVol {
private:
  Array<double,1> lvl;
  Array<double,1> decay;
public:
  inline ExponentialVol(const Array<double,1>& xlvl,const Array<double,1>& xdecay) : lvl(xlvl.copy()),decay(xdecay.copy()) { };
  inline ExponentialVol(double v,double a) : lvl(1),decay(1) { lvl = v; decay = a; };
  /// Returns the corresponding volatility of a state variable in a Gaussian HJM model where the forward rate volatility is given by *this - note that the sign of the mean reversion coefficient changes in this case.
  boost::shared_ptr<DeterministicAssetVol> component_vol(int i) const;
  virtual double volproduct(double t,double dt,const DeterministicAssetVol& xv) const;
  virtual Array<double,1> integral(double t,double dt) const;
  /// Dimension of the volatility vector.
  virtual int factors() const;
  virtual int type() const;
  inline double mean_reversion(int i) const { return decay(i); };
  virtual double volatility_level(int i,double T1,double T2) const;
  virtual double FwdBondVol(double t,double T1,double T2) const;
  double volproduct_ExponentialVol(double t,double dt,const Array<double,1>& xlvl,const Array<double,1>& xdecay) const;
  /// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
  virtual double bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another DeterministicVol object. */
  virtual double bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another ExponentialVol object. */
  virtual double bondbondvolproduct_ExponentialVol(double t,double dt,double T1,double T2,const ExponentialVol& xv) const;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      volatility (not bond volatility) given by another ExponentialVol object. */
  virtual double bondvolproduct_ExponentialVol(double t,double dt,double T,const ExponentialVol& xv) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual double A(double t,double T) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual Array<double,1> B(double t,double T) const;
  /// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
  virtual double var(int i,double t,double dt) const;
  /// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
  virtual double covar_dW(int i,double t,double dt) const;
  /** For Monte Carlo simulation: Covariance between the j-th component of the state variable increment and the j-th 
      component of the state variable increment for the DeterministicAssetVol xv */
  virtual double covar(int j,double t,double dt,const DeterministicAssetVol& xv) const;
  /** Integral over zero coupon bond volatility. 
  
      \f[ \int_{t}^{t+\Delta t}\sigma^*(s,T)ds \f] */
  virtual Array<double,1> bondvolintegral(double t,double dt,double T) const;
  /// Expected value of the state variable increments between t and t+dt under the time T forward measure.
  virtual Array<double,1> StateVariableMean(double t,double dt,double T) const;
  virtual double bondvolexponential(double t,double dt,double T,const Array<double,1>& dZ,const Array<double,1>& dW) const;
  virtual double z_integral(int i,double t,double dt) const;
  virtual Array<double,1> z_bondintegral(double t,double dt,double T,const DeterministicAssetVol& xmodel) const;
  virtual Array<double,1> z_volintegral(double t,double dt,double T,const DeterministicAssetVol& fxvol) const;
  Array<double,1> z_bondintegral_ExponentialVol(double t,double dt,double T,const ExponentialVol& bondvol) const;
  virtual double FwdFXexponential(double t,double dt,const Array<double,1>& dWj,const Array<double,1>& dZj0,const Array<double,1>& dZjk,DeterministicAssetVol* xv,DeterministicAssetVol* fxvol) const;
};

}

#endif
