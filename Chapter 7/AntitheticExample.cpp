/** \file  AntitheticExample.cpp
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
      -# Numeraire asset index.
         The default is -1.
      -# minimum MC paths.
         The default is 40.
      -# maximum MC paths.
         The default is 40.
*/
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i,j,k;
  try {
    // --------------- read command line parameters ---------------
    double S = 100.0;
    if (argc>1) S = atof(argv[1]);
    double r = 0.05;
    if (argc>2) r = atof(argv[2]);
    double sgm = 0.3;
    if (argc>3) sgm = atof(argv[3]);
    double mat = 1.5;
    if (argc>4) mat = atof(argv[4]);
    double K = 1.0;
    if (argc>5) K = atof(argv[5]);
    K *= S;
	double call_strike = K;
    int numeraire_index = -1;
    if (argc>6) numeraire_index = atoi(argv[6]);
    size_t minpaths = 40;
    if (argc>7) minpaths = atoi(argv[7]);
    size_t maxpaths = 40;
    if (argc>8) maxpaths = atoi(argv[8]);
    // --------------- create underlying asset ---------------
    Array<double,1> sgm1(1);
    sgm1 = sgm;
    ConstVol vol(sgm1);
    BlackScholesAsset stock(&vol,S,0.03);
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
    // --------------- closed form option price ---------------
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nvolatility: " << sgm << endl;
    double CFcall = stock.option(mat,K,r);
    cout << "Closed form call: " << CFcall << endl;
    double CFprice = CFcall;
    FlatTermStructure ts(r,0.0,mat+10.0);  // flat term structure for discounting
	// --------------- European call option price by Monte Carlo ---------------
	unsigned long n = minpaths;
	// instantiate random number generator
    ranlib::NormalUnit<double> normalRNG;
    RandomArray<ranlib::NormalUnit<double>,double> random_container2(normalRNG,1,1); // 1 factor, 1 time step
	// instantiate stochastic process
    GeometricBrownianMotion gbm(underlying);
	// 95% quantile for confidence interval
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	// boost functor to convert random variates to their antithetics (instantiated from template)
	boost::function<Array<double,2> (Array<double,2>)> antithetic = normal_antithetic<Array<double,2> >;
	// instantiate MCGatherer objects to collect simulation results
	MCGatherer<double> mcgatherer;
	MCGatherer<double> mcgatherer_antithetic;
	// instantiate MCPayoff object
	MCEuropeanCall mc_call(0,mat,0,call_strike);
	// instantiate MCMapping and bind to functor
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping2(mc_call,gbm,ts,numeraire_index);
	boost::function<double (Array<double,2>)> func2 = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mapping,&mc_mapping2,_1);
	// instantiate generic Monte Carlo algorithm object
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc2(func2,random_container2);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc2_antithetic(func2,random_container2,antithetic);
	cout << "European call option\nPaths,Closed form value,MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors,";
	cout << "Antithetic MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors,CI width,Antithetic CI width\n";
	// run Monte Carlo for different numbers of simulations
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc2.simulate(mcgatherer,n);
	  // half as many paths for antithetic
	  mc2_antithetic.simulate(mcgatherer_antithetic,n/2); 
	  cout << mcgatherer.number_of_simulations() << ',' << CFprice << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev() << ',' << (mcgatherer.mean()-CFprice)/mcgatherer.stddev() << ',';
	  cout << mcgatherer_antithetic.mean() << ',' << mcgatherer_antithetic.mean()-d*mcgatherer_antithetic.stddev() << ',' << mcgatherer_antithetic.mean()+d*mcgatherer_antithetic.stddev() << ',' << (mcgatherer_antithetic.mean()-CFprice)/mcgatherer_antithetic.stddev() << ',';
      cout << 2.0*d*mcgatherer.stddev() << ',' << 2.0*d*mcgatherer_antithetic.stddev() << ',' << endl;
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
