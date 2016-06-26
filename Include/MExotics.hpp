/** \file  MExotics.hpp
    \brief Header file declaring classes to price exotic options using the "Quintessential Option Pricing Formula" 
           of Skipper and Buchen (2003).
           Copyright 2006, 2007, 2010, 2013 by Erik Schlögl

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

#ifndef MEXOTICS_HPP
#define MEXOTICS_HPP
#include "MBinaryPayoffs.hpp"
#include "MCEngine.hpp"

namespace quantfin { namespace exotics {

/** Option to exchange one asset for another.

    Payoff: \f[ [K_1S_1(T)-K_2S_2(T)]^+ \f]
    */
class Margrabe {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff     part1;
  MBinaryPayoff     part2;
  MBinary*             M1;
  MBinary*             M2;
  double K1,K2;
public:
  Margrabe(const BlackScholesAsset& S1,const BlackScholesAsset& S2,double t,double T,double xK1,double xK2,const TermStructure& xts);
  ~Margrabe();
  inline double price() { return K1*M1->price(0)-K2*M2->price(0); };
  inline double price(const Array<double,1>& T,const Array<double,2>& history) { 
	return K1*M1->price(T,history,0)-K2*M2->price(T,history,0); };
  /// For American option pricing by Monte Carlo
  double early_exercise_payoff(const Array<double,1>& T,const Array<double,2>& history) const;
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

/** Standard Option.

    Payoff: \f[ [s\cdot(S(T)-K)]^+ \f]
    */
class StandardOption {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff     part1;
  MBinaryPayoff     part2;
  MBinary*             M1;
  MBinary*             M2;
  int                sign;
  double K;
public:
  StandardOption(const BlackScholesAsset& S,double t,double T,double xK,const TermStructure& xts,int xsign = 1);
  ~StandardOption();
  inline double price() { return sign*(M1->price(0)-K*M2->price(0)); };
  inline double price(double t,double S) { 
	return sign*(M1->price(t,S,0)-K*M2->price(t,S,0)); };
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

/** Power Option.

    Payoff: \f[ [S(T)^{\alpha}-K]^+ \f]
    */
class PowerOption {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff     part1;
  MBinaryPayoff     part2;
  MBinary*             M1;
  MBinary*             M2;
  double K;
public:
  PowerOption(const BlackScholesAsset& S,double alpha,double t,double T,double xK,const TermStructure& xts);
  ~PowerOption();
  inline double price() { return M1->price(0)-K*M2->price(0); };
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

/** Option on the discretely sampled geometric mean of an asset price process.

    Payoff: \f[ [\sqrt[n]{\prod_{i=1}^n S(T_i)}-K]^+ \f]
    */
class DiscreteGeometricMeanFixedStrike {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff       geo;
  MBinaryPayoff        BB;        
  MBinary*           Mgeo;
  MBinary*            MBB;
  double                K;
  double           factor;  ///< Factor for existing average.
public:
  DiscreteGeometricMeanFixedStrike(const BlackScholesAsset& S,const Array<double,1>& T,double xK,const TermStructure& xts,int number_of_observations_in_existing_average,double existing_average);
  DiscreteGeometricMeanFixedStrike(const GaussMarkovWorld& world,const Array<double,1>& T,double xK,int reportable_asset_index,int number_of_observations_in_existing_average,double existing_average);
  ~DiscreteGeometricMeanFixedStrike();
  inline double price() { return factor*Mgeo->price(0)-K*MBB->price(0); };
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

/** Option on the discretely sampled geometric mean of an asset price process.

    Payoff: \f[ [S(T_n)-K\sqrt[n]{\prod_{i=1}^n S(T_i)}]^+ \f]
    */
class DiscreteGeometricMeanFloatingStrike {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff       geo;
  MBinaryPayoff        AB;        
  MBinary*           Mgeo;
  MBinary*            MAB;
  double                K;
  double           factor;  ///< Factor for existing average.
public:
  DiscreteGeometricMeanFloatingStrike(const BlackScholesAsset& S,const Array<double,1>& T,double xK,const TermStructure& xts,int number_of_observations_in_existing_average,double existing_average);
  ~DiscreteGeometricMeanFloatingStrike();
  inline double price() { return MAB->price(0)-K*factor*Mgeo->price(0); };
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

class DiscreteBarrierOut {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff        AB;
  MBinaryPayoff        BB;        
  MBinary*            MAB;
  MBinary*            MBB;
  double                K;
  int             callput;
public:
  DiscreteBarrierOut(const BlackScholesAsset& S,  ///< Underlying asset.
                     const Array<double,1>& T,    ///< Barrier monitoring time points.
                     double xK,                   ///< Strike.
                     double barrier,              ///< Barrier.
                     const TermStructure& xts,    ///< (Deterministic) term structure of interest rates.
                     int xcallput = 1,            ///< Call (1) or put (-1). Call is default. 
                     int updown   = -1            ///< Up (1) or down (-1) option. Down is default.
                     );
  ~DiscreteBarrierOut();
  inline double price(unsigned long n = 100) { return callput*(MAB->price(n)-K*MBB->price(n)); };
  boost::shared_ptr<MCPayoffList> get_payoff() const;
  // Accessors
  inline const Array<double,2>& covariance_matrix() const { return MAB->covariance_matrix(); };
  inline const Array<double,1>& eigenvalues() const { return MAB->eigenvalues(); };
};

class ProductOption {
private:      
  const TermStructure& ts;  ///< (Deterministic) term structure of interest rates.
  MBinaryPayoff     part1;
  MBinaryPayoff     part2;
  MBinary*             M1;
  MBinary*             M2;
  double K1;
  int sign;
public:
  ProductOption(const BlackScholesAsset& S1,const BlackScholesAsset& S2,double t,double T,double xK1,const TermStructure& xts,int xsign = 1);
  ~ProductOption();
  inline double price() { return sign*(M1->price(0)-K1*M2->price(0)); };
  inline double price(const Array<double,1>& T,const Array<double,2>& history) { 
	return sign*(M1->price(T,history,0)-K1*M2->price(T,history,0)); };
  /// For American option pricing by Monte Carlo
  double early_exercise_payoff(const Array<double,1>& T,const Array<double,2>& history) const;
  boost::shared_ptr<MCPayoffList> get_payoff() const;
};

}}

#endif
