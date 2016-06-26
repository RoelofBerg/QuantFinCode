/** \file  MCDownandout.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2007, 2009, 2010 by Erik Schlögl

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

#include <iostream>
#include <cstdlib>
#include <random/normal.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/normal.hpp>
#include "BlackScholesAsset.hpp"
#include "PiecewiseVol.hpp"
#include "GeometricBrownianMotion.hpp"
#include "Payoff.hpp"
#include "StringForm.hpp"
#include "MBinary.hpp"
#include "MExotics.hpp"
#include "MCEngine.hpp"
#include "MCGeneric.hpp"
#include "QFRandom.hpp"

using namespace quantfin;

/// Direct implementation of class representing discounted payoff of a down-and-out barrier call option.
class down_and_out_call : public MCPayoff {
public:
  double strike,barrier;
  down_and_out_call(const Array<double,1>& T,int underlying_index,double xstrike,double xbarrier);
  /// Calculate discounted payoff. 
  virtual double operator()(const Array<double,1>& underlying_values, ///< Underlying values for the (asset,time) combinations in index Array.
	                        const Array<double,1>& numeraire_values   ///< Numeraire values for the dates in timeline Array.
							);
};

/// Constructor.
down_and_out_call::down_and_out_call(const Array<double,1>& T,int underlying_index,double xstrike,double xbarrier)
: MCPayoff(T,T.extent(firstDim)-1),strike(xstrike),barrier(xbarrier)
{
  firstIndex  idx;
  secondIndex jdx;
  index = idx * (jdx+1);
}

/// Calculated discounted payoff for a given path and numeraire.
double down_and_out_call::operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values)
{
  int i;
  bool indicator = true;
  double result  = 0.0;
  for (i=0;i<underlying_values.extent(firstDim);i++) if (underlying_values(i)<=barrier) indicator = false;
  if (indicator) result = numeraire_values(0)/numeraire_values(numeraire_values.extent(firstDim)-1) * std::max(0.0,underlying_values(underlying_values.extent(firstDim)-1)-strike);
  return result;
}
 
/** Test and demonstration for Monte Carlo simulation.

    Command-line arguments:
      -# initial stock price.
         The default is 100.
      -# interest rate.
         The default is 5%.
      -# volatility.
         The default is 30%.
      -# maturity.
         The default is 1.5.
      -# moneyness.
         The default is 1 (at the money).
      -# N, time line refinement.
         The default is 10.
      -# Numeraire asset index.
         The default is -1.
      -# minimum MC paths.
         The default is 100.
      -# maximum MC paths.
         The default is 100.
*/
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  try {
    // read command line parameters
    double S = 100.0;
    if (argc>1) S = atof(argv[1]);
    double r = 0.05;
    if (argc>2) r = atof(argv[2]);
    double sgm = 0.3;
    if (argc>3) sgm = atof(argv[3]);
    double maturity = 1.5;
    if (argc>4) maturity = atof(argv[4]);
    double K = 1.0;
    if (argc>5) K = atof(argv[5]);
    K *= S;
	double call_strike = K;
    int N = 10;
    if (argc>6) N = atoi(argv[6]);
    int numeraire_index = -1;
    if (argc>7) numeraire_index = atoi(argv[7]);
    size_t minpaths = 100;
    if (argc>8) minpaths = atoi(argv[8]);
    size_t maxpaths = 100;
    if (argc>9) maxpaths = atoi(argv[9]);
	// set up timeline - N is the number of time steps, mat is the maturity of the option
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = maturity/N;
    T = idx*dt;
	// set up asset - volatility level is sgm on first factor only
    Array<double,1> sgm1(2);
    sgm1 = sgm, 0.0;
    ConstVol vol(sgm1);
    BlackScholesAsset stock(&vol,S);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << maturity << "\nsgm: " << sgm << endl;
	// closed-form vanilla option prices
    double CFcall = stock.option(maturity,K,r);
    double CFput  = stock.option(maturity,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
 	// initialise vector of underlying assets
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
	// flat term structure for discounting
    FlatTermStructure ts(r,0.0,maturity+10.0);

	// Down and out call option price by Monte Carlo and in "closed form"
	MCGatherer<Array<double,1> > mcgatherer(2);
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	// univariate standard normal pseudo-random number generator
    ranlib::NormalUnit<double> normalRNG;
	// collection of independent random variates required to generate a single path
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,2,N);  // 2 factors, N time steps
	// asset price process
    GeometricBrownianMotion gbm(underlying);
	double barrier = 0.95*K;
	exotics::DiscreteBarrierOut down_and_out_call_option(stock,T,K,barrier,ts,1,-1); 
	// Convergence of "closed form" price
	cout << "DiscreteBarrierOut closed-form formula\nPrice: " << endl;
	i = minpaths;
	while (i<=maxpaths) {
	  // i is the number of simulations in the evaluation of the multivariate normal CDF by quasi-random Monte Carlo
	  cout << i << ',' << down_and_out_call_option.price(i) << endl;
	  i *= 2; }
	cout << "Covariance matrix dimension: " << down_and_out_call_option.covariance_matrix().extent(firstDim);
	cout << "\nEigenvalues: " << down_and_out_call_option.eigenvalues() << endl;
	// create Monte Carlo payoff objects - one using MBinaries and the other using the down_and_out_call class defined above
	boost::shared_ptr<MCPayoffList> dopayoff = down_and_out_call_option.get_payoff();
	boost::shared_ptr<MCPayoff> down_and_out_call_instance(new down_and_out_call(T,0,K,barrier));
	MCPayoffList both_options;
    both_options.push_back(down_and_out_call_instance);
    both_options.push_back(dopayoff);
	// MCMapping to map random numbers to asset price realisations to discounted payoffs
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_do(both_options,gbm,ts,numeraire_index);
	// mapping functor
	boost::function<Array<double,1> (Array<double,2>)> func_do = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mappingArray,&mc_mapping_do,_1);
	// collection of independent random variates required to generate a single path
    RandomArray<ranlib::NormalUnit<double>,double> random_container_do(normalRNG,gbm.factors(),gbm.number_of_steps()); 
	// generic Monte Carlo algorithm object
	MCGeneric<Array<double,2>,Array<double,1>,RandomArray<ranlib::NormalUnit<double>,double> > mcg(func_do,random_container_do);
	mcgatherer.reset();
	size_t n = minpaths;
	cout << "Paths,MC value 1,MC value 2\n";
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mcg.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << (mcgatherer.mean())(0) << ',' << (mcgatherer.mean())(1) << ',';
	  cout << (mcgatherer.mean())(0)-d*(mcgatherer.stddev())(0) << ',' << (mcgatherer.mean())(0)+d*(mcgatherer.stddev())(0) << ',' << endl;
	  n = mcgatherer.number_of_simulations(); }
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
