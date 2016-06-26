/** \file  GaussianHJM.hpp
    \brief Header file declaring class for a Gauss/Markov Heath/Jarrow/Morton (HJM) model.
           Copyright 2003, 2005, 2006, 2009, 2011 by Erik Schlögl

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

#ifndef GAUSSIANHJM_H
#define GAUSSIANHJM_H
#include <random/normal.h>
#include <boost/math/distributions/normal.hpp>
#include "TermStructure.hpp"
#include "DeterministicVol.hpp"
#include "BlackScholesAsset.hpp"

namespace quantfin {

using blitz::Array;
using blitz::firstDim;
using ranlib::NormalUnit;

/** \brief Class for a Gauss/Markov Heath/Jarrow/Morton (HJM) model.

    A particular instance of such a model is characterised by a volatility function for the 
    instantaneous forward rates and by an initial term structure of interest rates.
  */
class GaussianHJM {
private:
  DeterministicAssetVol*    v;  ///< Pointer to deterministic volatility function object.
  TermStructure*    initialTS;  ///< Pointer to term structure object containing the initial term structure of interest rates.
  boost::math::normal       N;
  int                    GH_n;  ///< Refinement for Gaussian quadrature (where numerical integration is needed). The default is 40.
  /// Class for rootsearch required to implement coupon bond pricing in the single factor model.
  class r_star { 
  private:
    const GaussianHJM& model;
    const Array<double,1>&  tenor;
	const Array<double,1>&  coupon;
  public:
	inline r_star(const GaussianHJM& xmodel,const Array<double,1>& xtenor,const Array<double,1>& xcoupon) 
	  : model(xmodel),tenor(xtenor),coupon(xcoupon) { };
	virtual double operator()(double r) const;
  };
  /// Class for numerical integration required to implement coupon bond pricing in the two-factor model.
  class cb_option {
  private:
    const GaussianHJM& model;
    const Array<double,1>&  tenor;
	Array<double,1> coupon;
	Array<double,1> tmpB;
	double e0,e1,v0,v1;
	double fix_point_problem(double z2);
	double current_z1;
	double K;
	int sign;
  public:
    cb_option(const GaussianHJM& xmodel,const Array<double,1>& xtenor,const Array<double,1>& xcoupon,double strike,int xsign);
    double two_factor_CBoption_integrand(double z);
  };
public:
  /// Constructor.
  inline GaussianHJM(DeterministicAssetVol* xv,  ///< Pointer to deterministic volatility function object.
                     TermStructure* ini          ///< Pointer to term structure object containing the initial term structure of interest rates.
                     ) : v(xv),initialTS(ini),GH_n(40) { };
  inline void set_quadrature_refinement(int n) { GH_n = n; };
  inline double time_horizon() const { return (initialTS->timeline())(initialTS->timeline().extent(firstDim)-1); };
  /** Number of factors, i.e. dimension of the driving Brownian motion.
      Passed through to the deterministic volatility function object DeterministicVol. */
  inline int factors() const { return v->factors(); }; 
  /// Access the volatility function.
  inline const DeterministicAssetVol& volatility_function() const { return *v; };
  /** Initial continuously compounded forward yields for accrual periods defined by the vector of dates T.
      Passed through to the initial TermStructure object. */
  inline Array<double,1> initial_forward_rates(const Array<double,1>& T) const { return initialTS->ccfwd(T); };
  /// Calculate a zero coupon bond price given the state variables.
  double bond(const Array<double,1> &W,double t,double ttm) const;
  /** Initial zero coupon bond price for a given maturity.
      Passed through to the initial TermStructure object. */
  inline double bond(double mat) const { return initialTS->operator()(mat); };
  /// Calculate the ("time zero") value of a zero coupon bond option.
  double ZCBoption(double T1,double T2,double K,int sign = 1) const;
  /// Calculate the ("time zero") value of a coupon bond option.
  double CBoption(const Array<double,1>& tenor,const Array<double,1>& coupon,double K,int sign = 1) const;
  /// Calculate the ("time zero") value of swaption.
  double swaption(const Array<double,1>& tenor,  ///< Swap tenor.
	              double K,                 ///< Strike (fixed side) level.
				  int sign = 1              ///< 1 = receiver swaption, -1 = payer swaption.
				  ) const;
  /// Calculate the ("time zero") value of a caplet.
  double caplet(double mat,double delta,double K) const;
  /// Calculate the ("time zero") value of an equity option.
  double option(const BlackScholesAsset& S,double mat,double K,int sign = 1) const;
  /// Calculate the ("time zero") value of an FX option.
  double FXoption(double X0,double T1,double K,const GaussianHJM& fmodel,const DeterministicAssetVol& fxvol,int sign = 1) const;
  /// Swap netting the difference between the foreign and domestic floating rates (reverse if sign = -1), paid on a domestic notional.
  double DiffSwap(const Array<double,1>& T,const GaussianHJM& fmodel,const DeterministicAssetVol& fxvol,int sign = 1) const;
  /// Foreign caplet with payoff converted to domestic currency at a guaranteed exchange rate \e X0.
  double QuantoCaplet(double X0,double T,double delta,double lvl,const GaussianHJM& fmodel,const DeterministicAssetVol& fxvol) const;
  double QuantoCaplet_old(double X0,double T,double delta,double lvl,const GaussianHJM& fmodel,const DeterministicAssetVol& fxvol) const;
  /** Integral over zero coupon bond volatility. 
  
      \f[ \int_{t}^{t+\Delta t}\sigma^*(s,T)ds \f]
      Passed through to the deterministic volatility function object DeterministicVol. */
  inline Array<double,1> bondvolintegral(double t,double dt,double T) const { return v->bondvolintegral(t,dt,T); };
  /// Expected value of the state variable increments between t and t+dt under the time T forward measure.
  inline Array<double,1> StateVariableMean(double t,double dt,double T) const { return v->StateVariableMean(t,dt,T); };
  /** Expected value of the state variable increments between t and t+dt under the time T forward measure
      associated with *xmodel, where the volatility of the spot exchange rate (of units of the *xmodel currency
      per unit of *this currency) is fxvol. */
  Array<double,1> StateVariableMean(double t,double dt,double T,const GaussianHJM* xmodel,const DeterministicAssetVol& fxvol) const;
  /// Expected value of the state variable increments between t and t+dt under the physical measure, where the market price of risk is given by mpr.
  inline Array<double,1> StateVariableMean(double t,double dt,const DeterministicAssetVol& mpr) const;
  /** Generate the state variable increments in a way consistent with a given set of Brownian motion increments.
      Passed through to the deterministic volatility function object DeterministicVol. */
  inline Array<double,1> StateVariableTransform(const Array<double,1>& dWi,double t,double dt,NormalUnit<double>& rnd) const;
  /** \brief 
      Covariance between the j-th component of the Brownian motion increment and the j-th 
      component of the state variable increment.

      Passed through to the deterministic volatility function object DeterministicVol. */
  inline double covar_dW(int j,double t,double dt) const { return v->covar_dW(j,t,dt); };
  /** \brief 
      Variance of the j-th component of the state variable increment.

      Passed through to the deterministic volatility function object DeterministicVol. */
  inline double var(int j,double t,double dt) const { return v->var(j,t,dt); };
  /** \brief 
      Covariance between the j-th component of the state variable increment and the j-th 
      component of the state variable increment of the GaussianHJM object \e m.

      Passed through to the deterministic volatility function object DeterministicVol. */
  inline double covar(int j,double t,double dt,const GaussianHJM& m) const { return v->covar(j,t,dt,*(m.v)); };
  /** Integral over the scalar product between the bond volatility given by (*this) and a deterministic asset volatility.
  
      Passed through to the deterministic volatility function object DeterministicVol. */
  inline double bondvolproduct(double t,double dt,double bondmat,const DeterministicAssetVol& xv) const
    { return v->bondvolproduct(t,dt,bondmat,xv); };
  /** Integral over the scalar product between the bond volatility given by (*this) and a 
      bond volatility given by another DeterministicVol object. 
  
      Passed through to the deterministic volatility function object DeterministicVol. */
  inline double bondbondvolproduct(double t,double dt,double T1,double T2,const GaussianHJM& xm) const
    { return v->bondbondvolproduct(t,dt,T1,T2,*(xm.v)); };
  inline double FwdFXexponential(double t,double dt,const Array<double,1>& dWj,const Array<double,1>& dZj0,const Array<double,1>& dZjk,const GaussianHJM* xm,DeterministicAssetVol* fxvol) const { return v->FwdFXexponential(t,dt,dWj,dZj0,dZjk,xm->v,fxvol); };
};

class GaussMarkovTermStructure : public TermStructure {
private:
  double                   today;
  const GaussianHJM&         hjm;
  Array<double,1> statevariables;
public:
  /// Constructor.
  GaussMarkovTermStructure(double xtoday,                         ///< Current time for which this term structure applies.
                           const GaussianHJM& xhjm,               ///< Underlying Gauss/Markov Heath/Jarrow/Morton model.
                           const Array<double,1>& xstatevariables ///< State variables describing the state of the model at the current time.
 						   );
  /// Virtual copy constructor
  virtual GaussMarkovTermStructure* pointer_to_copy() const;
  /// Reinitialise (if necessary) when zero coupon bond prices have changed
  virtual void reinitialise();
  /** Function returning the interpolated value, which for term structures of
      interest rates is the maturity t, time T(0) forward zero coupon bond price. */
  virtual double operator()(double t) const;    
  /// For nested Visitor pattern.
  virtual void accept(QFNestedVisitor& visitor) const;
  virtual const std::string& name() const;
};

inline Array<double,1> GaussianHJM::StateVariableTransform(const Array<double,1>& dWi,double t,double dt,NormalUnit<double>& rnd) const 
{ 
  int i;
  Array<double,1> dZi(dWi.extent(firstDim));
  double sqrtdt = std::sqrt(dt);
  for (i=0;i<dWi.extent(firstDim);i++) {
    double sd_dZi = std::sqrt(v->var(i,t,dt));
    double rho    = v->covar_dW(i,t,dt) / (sd_dZi*sqrtdt);
    dZi(i) = sd_dZi * (rho*dWi(i)/sqrtdt + std::sqrt(1.0-rho*rho)*rnd.random()); }
  return dZi;
}

inline Array<double,1> GaussianHJM::StateVariableMean(double t,double dt,const DeterministicAssetVol& mpr) const
{
  int i;
  Array<double,1> mu(factors());
  if (!mpr.get_volatility_level(t,t+dt,mu)) throw std::logic_error("Volatility not constant over required interval");
  for (i=0;i<factors();i++) mu(i) *= v->z_integral(i,t,dt);
  return mu;
}

}

#endif
