/** \file  PiecewiseVol.hpp
    \brief Header file declaring deterministic volatility function classes with piecewise constant parameters.
           Copyright 2005, 2006, 2010 by Erik Schlögl

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

#ifndef PIECEWISEVOL_H
#define PIECEWISEVOL_H
#include "DeterministicVol.hpp"

namespace quantfin {

class PiecewiseConstVol : public DeterministicAssetVol {
private:
  Array<double,1> timeline; ///< Time line defining calendar time segments on which parameters are constant.
  Array<double,2> v;        ///< Matrix of volatility scale parameters. Dimensions: time segment (row) by driving factor (column).
  Array<double,1> vol_sq;   ///< Scalar products of the volatility scale parameter vector for each time segment with itself.
public:
  /// Constructor.
  PiecewiseConstVol(const Array<double,1>& xT, ///< Time line defining calendar time segments on which parameters are constant.
                    const Array<double,2>& xv  ///< Matrix of volatility scale parameters. Dimensions: time segment (row) by driving factor (column).
                    );
  boost::shared_ptr<DeterministicAssetVol> component_vol(int i) const;
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

}

#endif


