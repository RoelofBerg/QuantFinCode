/** \file  GaussianEconomy.hpp
    \brief Header file defining classes for representing a set of "economies" where asset prices are driving by multidimensional 
           geometric Brownian motion.
           Copyright 2010, 2011, 2012, 2013 by Erik Schlögl

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

#ifndef GAUSSIANECONOMY_HPP
#define GAUSSIANECONOMY_HPP
#include <stdexcept>
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include "BlackScholesAsset.hpp"
#include "GaussianHJM.hpp"
#include "MultivariateNormal.hpp"


namespace quantfin { 

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;
using blitz::Range;
using blitz::toEnd;
using blitz::firstIndex;
using blitz::secondIndex;

/** Class representing assets quoted in a particular "currency." These assets include a term structure of zero coupon bond prices
    and assets with deterministic continuous dividend yield (which may be zero). */
class GaussianEconomy {
public:
  GaussianEconomy(std::vector<boost::shared_ptr<BlackScholesAsset> >& xunderlying,
	              boost::shared_ptr<DeterministicAssetVol> xv,boost::shared_ptr<TermStructure> xinitialTS);
  GaussianEconomy(const char* path);
  std::vector<boost::shared_ptr<BlackScholesAsset> >        underlying; ///< Vector of pointers to underlying assets.
  boost::shared_ptr<DeterministicAssetVol>                           v; ///< Pointer to deterministic volatility function object (for Gaussian term structure dynamics).
  std::vector<boost::shared_ptr<DeterministicAssetVol> > component_vol; ///< Vector of pointers to DeterministicAssetVol objects where all components except one are zero - corresponds to volatilies of the state variables z.
  boost::shared_ptr<TermStructure>                           initialTS; ///< Pointer to term structure object containing the initial term structure of interest rates.
  boost::shared_ptr<GaussianHJM>                                   hjm; ///< for closed-form solutions
};

/** Class representing a multicurrency Gauss/Markov HJM term structure model. Structure parallels that of the class GeometricBrownianMotion
    (which is for geometric Brownian motion of multiple assets with deterministic volatility under deterministic interest rates).
*/
class GaussMarkovWorld {
private:
  std::vector<boost::shared_ptr<GaussianEconomy> >  economies; ///< "Economies" in different currencies. First entry is taken to be the reference (domestic) currency.
  std::vector<boost::shared_ptr<DeterministicAssetVol> > vols; ///< Volatilities of exchange rates (in terms of units of "domestic" currency per unit of each of the "foreign" currencies.
  Array<double,1>          initial_exchange_rates;
  Array<double,1>      terminal_fwd_exchange_rate;
  Array<double,2>              terminal_fwd_asset;
  Array<double,1>*                              T; ///< Process timeline.
  Array<double,2>*                state_variables; ///< Simulation path is represented internally as a (number of state variables) x (number of process time points) Array
  Array<double,2>*          state_variable_drifts; ///< (number of state variables) x (number of process time steps) Array
  int                   number_of_state_variables;
  int                         diffusion_dimension; ///< Dimension of the driving Brownian motion
  int                                    max_rank; ///< maximum rank of covariance matrices of the multivariate normal distributions of state variable increments
  int                        currency_start_index; ///< Index of first exchange rate state variable
  Array<double,1>    termstructure_statevariables;
  Array<int,1>                economy_start_index; ///< index of first state variable relevant to each economy
  std::vector<boost::shared_ptr<MultivariateNormal> > mvn; ///< multivariate normal distributions of state variable increments
  int                           current_numeraire;
  int                  numeraire_reportable_index; ///< index of the numeraire asset in the list of reportable assets (if applicable)
  int                    numeraire_reportable_ZCB; ///< index of the domestic terminal zero coupon bond in the list of reportable assets (if applicable), to convert terminal forward to spot value of the numeraire
  bool is_asset_index(int i) const;
  inline int economy_index(int j) const;  
  void initialise();
  bool check_inputs();
  const DeterministicAssetVol* get_vol(int index);
  double zsum,z2sum; // for debugging
  int zn; // for debugging
  void propagate_state_variables(const Array<double,2>& x);
public:
  class reportable {
  public:
    int currency_index;
	int    asset_index;
	double    maturity; };
  std::vector<reportable> reportable_list;
  /// Constructor
  GaussMarkovWorld(std::vector<boost::shared_ptr<GaussianEconomy> >& xeconomies,
	               std::vector<boost::shared_ptr<DeterministicAssetVol> >& xvols,
				   Array<double,1> xinitial_exchange_rates);
  /// Construct from a CSV file specification
  GaussMarkovWorld(const char* path);
  /// Destructor
  ~GaussMarkovWorld();
  /// Determine which asset values should be reported back by operator()
  int set_reporting(int currency_index,int asset_index,double maturity = -1.0);
  /// To do: Set drifts, etc. to simulate under a particular choice of numeraire.
  bool set_numeraire(int num);
  /// Query the dimension of the process.
  int dimension() const;
  /// To do: Query the number of random variables (per period) driving the process.
  inline int factors() const { return max_rank; };
  /// To do: Set process timeline - need var/covar matrix for state variable increments.
  bool set_timeline(const Array<double,1>& timeline);
  /// Get timeline of asset price realisations reported by operator().
  inline const Array<double,1>& get_timeline() const { return *T; };
  /// Get number of steps in process time discretisation.
  inline int number_of_steps() const { return T->extent(firstDim)-1; };
  /// Generating a realisation of the process under the martingale measure associated with deterministic bond prices is not applicable in a stochastic term structure model.
  void operator()(Array<double,2>& underlying_values,const Array<double,2>& x,const TermStructure& ts);
  /// Generate a realisation of the process under the martingale measure associated with a given numeraire asset. underlying_values is an asset x (time points) Array.
  void operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const Array<double,2>& x,const TermStructure& ts,int numeraire_index);
  /** Generate a realisation of the process under the martingale measure associated with a given numeraire asset, 
      returning the state variables. underlying_values is an asset x (time points) Array. */
  const Array<double,2>& generate_state_variables(Array<double,1>& numeraire_values,const Array<double,2>& x,int numeraire_index);
  /// Get the current (simulated term structure at the i-th time in the time line, for the j-th currency
  GaussMarkovTermStructure get_TermStructure(int i,int j);
  inline const Array<double,1>& get_initial_exchange_rates() const { return initial_exchange_rates; };
  inline const std::vector<boost::shared_ptr<GaussianEconomy> >&  get_economies() const { return economies; };
  inline const std::vector<boost::shared_ptr<DeterministicAssetVol> >& get_FXvolatilities() const { return vols; };
  inline double get_forward_exchange_rate(int i,double mat) const 
  { return initial_exchange_rates(i-1)*(*(economies[i]->initialTS))(mat)/(*(economies[0]->initialTS))(mat); };
  inline double get_terminal_forward_asset(int currency,int asset) const { return terminal_fwd_asset(currency,asset); };
  inline double time_horizon() const { return (*T)(T->extent(firstDim)-1); };
};

}
 
#endif

