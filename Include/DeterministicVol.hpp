/** \file  DeterministicVol.hpp
    \brief Header file declaring deterministic volatility function classes.
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

#ifndef DETERMINISTICVOL_HPP
#define DETERMINISTICVOL_HPP
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

namespace quantfin {

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;

/** Abstract base class which defines additional functionality for deterministic volatility functions.
  *
  * It is used in particular for asset price processes modelled as geometric Brownian motion.
  */
class DeterministicAssetVol {
public:
  /// Returns the corresponding volatility of a state variable in a Gaussian HJM model where the forward rate volatility is given by *this.
  virtual boost::shared_ptr<DeterministicAssetVol> component_vol(int i) const = 0;
  /// The integral over the scalar product between two volatility vectors.
  virtual double volproduct(double t,double dt,const DeterministicAssetVol& xv) const = 0;
  virtual Array<double,1> integral(double t,double dt) const = 0;
  /// Dimension of the volatility vector.
  virtual int factors() const = 0;
  virtual int type() const;
  /// Division of time period into segments over which volatility is constant.
  virtual Array<double,1> segments(double t,double dt) const;
  /// Get constant volatility levels over time interval [t,T]. Returns false if volatility is not constant over this interval.
  virtual bool get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const;
  /** Integral over the square of the forward zero coupon bond volatility.
    * The value of the integral \f[ \sqrt{\int_t^{T_1}(\sigma^*(u,T_2)-\sigma^*(u,T_1))^2du} \f] is
    * calculated in classes derived from this abstract base class.
    */
  virtual double FwdBondVol(double t,double T1,double T2) const = 0;
  /// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
  virtual double bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const = 0;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another DeterministicVol object. */
  virtual double bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const = 0;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual double A(double t,double T) const = 0;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual Array<double,1> B(double t,double T) const = 0;
  /// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
  virtual double var(int i,double t,double dt) const = 0;
  /// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
  virtual double covar_dW(int i,double t,double dt) const = 0;
  /** For Monte Carlo simulation: Covariance between the j-th component of the state variable increment and the j-th 
      component of the state variable increment for the DeterministicAssetVol xv */
  virtual double covar(int j,double t,double dt,const DeterministicAssetVol& xv) const;
  /** Integral over zero coupon bond volatility. 
  
      \f[ \int_{t}^{t+\Delta t}\sigma^*(s,T)ds \f] */
  virtual Array<double,1> bondvolintegral(double t,double dt,double T) const = 0;
  /// Expected value of the state variable increments between t and t+dt under the time T forward measure.
  virtual Array<double,1> StateVariableMean(double t,double dt,double T) const = 0;
  virtual double bondvolexponential(double t,double dt,double T,const Array<double,1>& dZ,const Array<double,1>& dW) const;
  virtual double z_integral(int i,double t,double dt) const;
  virtual Array<double,1> z_bondintegral(double t,double dt,double T,const DeterministicAssetVol& xmodel) const;
  virtual Array<double,1> z_volintegral(double t,double dt,double T,const DeterministicAssetVol& fxvol) const;
  virtual double FwdFXexponential(double t,double dt,const Array<double,1>& dWj,const Array<double,1>& dZj0,const Array<double,1>& dZjk,DeterministicAssetVol* xv,DeterministicAssetVol* fxvol) const;
};

/** The deterministic asset volatility that results as the difference between two deterministic
    asset volatilities. This is the volatility of the quotient of two assets, which have deterministic
    volatility.
  */
class DeterministicAssetVolDiff : public DeterministicAssetVol {
private:
  const DeterministicAssetVol& v1;
  const DeterministicAssetVol& v2;
  boost::shared_ptr<DeterministicAssetVol> pv1,pv2;
public:
  /// Constructor.
  inline DeterministicAssetVolDiff(const DeterministicAssetVol& xv1,const DeterministicAssetVol& xv2);
  inline DeterministicAssetVolDiff(boost::shared_ptr<DeterministicAssetVol> xv1,boost::shared_ptr<DeterministicAssetVol> xv2);
  virtual boost::shared_ptr<DeterministicAssetVol> component_vol(int i) const;
  /// Get constant volatility levels over time interval [t,T]. Returns false if volatility is not constant over this interval.
  virtual bool get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const;
  /// The integral over the scalar product between two volatility vectors.
  virtual double volproduct(double t,double dt,const DeterministicAssetVol& xv) const;
  virtual Array<double,1> integral(double t,double dt) const;
  /// Dimension of the volatility vector.
  virtual int factors() const;
  /// Division of time period into segments over which volatility is constant.
  virtual Array<double,1> segments(double t,double dt) const;
  virtual double FwdBondVol(double t,double T1,double T2) const;
  /// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
  virtual double bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another DeterministicVol object. */
  virtual double bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual double A(double t,double T) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual Array<double,1> B(double t,double T) const;
  /// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
  virtual double var(int i,double t,double dt) const;
  /// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
  virtual double covar_dW(int i,double t,double dt) const;
  /** Integral over zero coupon bond volatility. 
  
      \f[ \int_{t}^{t+\Delta t}\sigma^*(s,T)ds \f] */
  virtual Array<double,1> bondvolintegral(double t,double dt,double T) const;
  /// Expected value of the state variable increments between t and t+dt under the time T forward measure.
  virtual Array<double,1> StateVariableMean(double t,double dt,double T) const;
};

inline DeterministicAssetVolDiff::DeterministicAssetVolDiff(const DeterministicAssetVol& xv1,const DeterministicAssetVol& xv2) 
  : v1(xv1),v2(xv2) 
{ 
  if (xv1.factors()!=xv2.factors()) throw std::logic_error("Dimension mismatch in DeterministicAssetVolDiff");  
}

inline DeterministicAssetVolDiff::DeterministicAssetVolDiff(boost::shared_ptr<DeterministicAssetVol> xv1,boost::shared_ptr<DeterministicAssetVol> xv2)
  : v1(*xv1),v2(*xv2),pv1(xv1),pv2(xv2) 
{ 
  if (xv1->factors()!=xv2->factors()) throw std::logic_error("Dimension mismatch in DeterministicAssetVolDiff");  
}

class ConstVol : public DeterministicAssetVol {
private:
  Array<double,1> lvl;
public:
  inline ConstVol(const Array<double,1>& xlvl) : lvl(xlvl.copy()) { };
  inline ConstVol(double v) : lvl(1) { lvl = v; };
  /// Returns the corresponding volatility of a state variable in a Gaussian HJM model where the forward rate volatility is given by *this.
  virtual boost::shared_ptr<DeterministicAssetVol> component_vol(int i) const;
  /// Get constant volatility levels over time interval [t,T]. Returns false if volatility is not constant over this interval.
  virtual bool get_volatility_level(double t,double T,Array<double,1>& vol_lvl) const;
  virtual double volproduct(double t,double dt,const DeterministicAssetVol& xv) const;
  virtual Array<double,1> integral(double t,double dt) const;
  /// Dimension of the volatility vector.
  virtual int factors() const;
  virtual double FwdBondVol(double t,double T1,double T2) const;
  /// Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
  virtual double bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const;
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another DeterministicVol object. */
  virtual double bondbondvolproduct(double t,double dt,double T1,double T2,const DeterministicAssetVol& xv) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual double A(double t,double T) const;
  /// Function needed for the exponential-affine representation of zero coupon bond prices.
  virtual Array<double,1> B(double t,double T) const;
  /// For Monte Carlo simulation: variance of increment from t to dt of state variable i.
  virtual double var(int i,double t,double dt) const;
  /// For Monte Carlo simulation: covariance of increment from t to dt of state variable i with increment of Brownian motion \f$ \Delta W^{(i)} \f$.
  virtual double covar_dW(int i,double t,double dt) const;
  /** Integral over zero coupon bond volatility. 
  
      \f[ \int_{t}^{t+\Delta t}\sigma^*(s,T)ds \f] */
  virtual Array<double,1> bondvolintegral(double t,double dt,double T) const;
  /// Expected value of the state variable increments between t and t+dt under the time T forward measure.
  virtual Array<double,1> StateVariableMean(double t,double dt,double T) const;
  virtual double z_integral(int i,double t,double dt) const;
  virtual Array<double,1> z_bondintegral(double t,double dt,double T,const DeterministicAssetVol& xmodel) const;
  virtual Array<double,1> z_volintegral(double t,double dt,double T,const DeterministicAssetVol& fxvol) const;
};

}

#endif


