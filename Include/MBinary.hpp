/** \file  MBinary.hpp
    \brief Header file declaring a class for the "Quintessential Option Pricing Formula" of Skipper and Buchen (2003).
           Copyright 2006, 2008, 2010, 2013 by Erik Schlögl

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

#ifndef MBINARY_HPP
#define MBINARY_HPP
#include <vector>
#include "MultivariateNormal.hpp"
#include "BlackScholesAsset.hpp"
#include "TermStructure.hpp"
#include "MCEngine.hpp"
#include "GaussianEconomy.hpp"

namespace quantfin {

using blitz::Array;
using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;

class MBinaryPayoff : public MCPayoff {
public:
  std::vector<const BlackScholesAsset*>  underlying;  ///< Vector of pointers to underlying assets.
  const TermStructure&                           ts;  ///< (Deterministic) term structure of interest rates.
  Array<double,1>                             alpha;  ///< Payoff powers.
  Array<int,2>                                    S;  ///< Exercise indicators.
  Array<double,2>                                 A;  ///< Exercise condition matrix.  
  Array<double,1>                                 a;  ///< Strike vector.
  double                                   notional;
  inline MBinaryPayoff(const TermStructure&     xts,  ///< (Deterministic) term structure of interest rates.
                       int                ntimeline,  ///< Number of points on the time line collecting all event dates.
                       int         payoff_dimension,  ///< Number of (asset,time) combinations.
                       int       exercise_dimension,  ///< Number of indicator functions.
					   double xnotional = 1.0
                       ) : MCPayoff(ntimeline-1,payoff_dimension),alpha(payoff_dimension),
                           ts(xts),S(exercise_dimension,exercise_dimension),A(exercise_dimension,payoff_dimension),
                           a(exercise_dimension),notional(xnotional) { };
  /// Calculate discounted payoff. 
  virtual double operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values);
};
  
class MBinary {
private:
  // Inputs
  int                              payoff_dimension;  ///< Number of (asset,time) combinations.
  int                            exercise_dimension;
  int                              number_of_assets;
  bool                                    worthless;
  const std::vector<const BlackScholesAsset*>& underlying;  ///< Vector of pointers to underlying assets.
  const TermStructure&                           ts;  ///< (Deterministic) term structure of interest rates.
  Array<double,1>                          timeline;  ///< Time line collecting all event dates.
  Array<int,2>                                index;  ///< A 2 x N matrix of indices, where each column represents the 
                                                      ///< indices of an (asset,time) combination affecting the payoff.
  Array<double,1>                             alpha;  ///< Payoff powers.
  Array<int,2>                                    S;  ///< Exercise indicators.
  Array<double,2>                                 A;  ///< Exercise condition matrix.  
  Array<double,1>                                 a;  ///< Strike vector.
  double                                   notional;
  // Values pre-calculated in the constructor
  Array<double,1>  mu;  ///< Vector of drifts.
  Array<double,2> sgm;  ///< N x N covariance matrix of logarithmic asset prices at event dates.
  double         beta;
  Array<double,1>  Sd;
  Array<double,2> SCS;
  Array<double,1>   x;  ///< Initial prices of the underlying assets to match payoff vector - thus repeat entries are possible.
  MultivariateNormal* mvrnd;
  // Saved values to reduce calculations
  double P_t,I_t;
  unsigned long current_n;
  void initialise();
  void initialise(const Array<double,1>& T,const Array<double,2>& history);
  double dPdX(int i);
  double dIdX(int i);
public:    
  MBinary(const std::vector<const BlackScholesAsset*>& xunderlying,  ///< Vector of pointers to underlying assets.
          const TermStructure&                           xts,  ///< (Deterministic) term structure of interest rates.
          const Array<double,1>&                   xtimeline,  ///< Time line collecting all event dates.
          const Array<int,2>&                         xindex,  ///< A 2 x N matrix of indices, where each column represents the 
                                                               ///< indices of an (asset,time) combination affecting the payoff.
          const Array<double,1>&                      xalpha,  ///< Payoff powers.
          const Array<int,2>&                             xS,  ///< Exercise indicators.
          const Array<double,2>&                          xA,  ///< Exercise condition matrix.  
          const Array<double,1>&                          xa,  ///< Strike vector.
		  double xnotional = 1.0
          );
  MBinary(MBinaryPayoff& xpayoff);
  MBinary(const GaussMarkovWorld& world,MBinaryPayoff& xpayoff);
  MBinary(const MBinary& original,const Array<double,1>& T,const Array<double,2>& history);
  ~MBinary();
  boost::shared_ptr<MBinaryPayoff> get_payoff() const;
  double price(unsigned long n = 1000000);
  double price(double t,double single_underlying,unsigned long n = 1000000);
  double price(const Array<double,1>& T,const Array<double,2>& history,unsigned long n = 1000000);
  // "Greeks"
  /// Delta with respect to the i-th asset.
  double delta(int i,unsigned long n = 1000000);
  // Accessors
  inline const Array<double,2>& covariance_matrix() const { return SCS; };
  inline const Array<double,1>& eigenvalues() const { return mvrnd->eigenvalues(); };
};

}

#endif

