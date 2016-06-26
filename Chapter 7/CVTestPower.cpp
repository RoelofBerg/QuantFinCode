/** \file  CVTestPower.cpp
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

/** Test and demonstration for control variates.

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
         The default is 40.
      -# maximum MC paths.
         The default is 40.
      -# Switch for antithetic.
         The default is 0 (off).
*/
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  try {
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
	K *= 0.5;
	double call_strike = K;
    int N = 10; 
    if (argc>6) N = atoi(argv[6]);
    int numeraire_index = -1;
    if (argc>7) numeraire_index = atoi(argv[7]);
    size_t minpaths = 40;
    if (argc>8) minpaths = atoi(argv[8]);
    size_t maxpaths = 40;
    if (argc>9) maxpaths = atoi(argv[9]);
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
    Array<double,1> sgm1(2);
    sgm1 = sgm;
    ConstVol vol(sgm1);
    FlatTermStructure ts(r,0.0,mat+10.0);

    ranlib::NormalUnit<double> normalRNG;
	MCGatherer<double> mcgatherer;
	// boost functor to convert random variates to their antithetics (instantiated from template)
	boost::function<Array<double,2> (Array<double,2>)> antithetic = normal_antithetic<Array<double,2> >;
	MCGatherer<double> mcgatherer_antithetic;
	unsigned long n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);

	// Power option price by Monte Carlo 
	Array<double,1> gsgm(1);
	gsgm = sgm;
	ConstVol gvol(gsgm);
    BlackScholesAsset gasset(&gvol,S);
	// Fixed strike
	double alpha = 0.1;
	call_strike = std::pow(call_strike,alpha);
	exotics::PowerOption powopt(gasset,alpha,T(0),T(T.extent(firstDim)-1),call_strike,ts);
	double CFpow = powopt.price();
	cout << "Closed form value of power option: " << CFpow << "\nStrike: " << call_strike << "\nPower: " << alpha << endl;
	boost::shared_ptr<MCPayoffList> powpayoff = powopt.get_payoff();
    std::vector<const BlackScholesAsset*> g_underlying;
    g_underlying.push_back(&gasset);
    GeometricBrownianMotion ggbm(g_underlying);
	ggbm.set_timeline(T);
	cout << ggbm.number_of_steps() << endl;
    RandomArray<ranlib::NormalUnit<double>,double> random_container_g(normalRNG,ggbm.factors(),ggbm.number_of_steps()); 

	// Power option price by Monte Carlo - including using standard options as control variates
	// Build objects for control variate simulation
	Array<double,1> cv_values(5),cv_strikes(5);
	cv_strikes = 0.75, 1.0, 1.25, 1.5, 2.5;
	cv_strikes *= K;
	for (i=0;i<5;i++) cv_values(i) = gasset.option(T(T.extent(firstDim)-1)-T(0),cv_strikes(i),r);
	cout << "Strikes: " << cv_strikes << endl;

	// using MCGatherer<Array<double,1> > specialisation
	MCGatherer<Array<double,1> > cvgatherer2(6);
	MCPayoffList                        cvlist;
	cvgatherer2.set_control_variate(true);
	cvlist.push_back(powpayoff);
	for (i=0;i<5;i++) {
	  boost::shared_ptr<MCEuropeanCall> callpayoff(new MCEuropeanCall(T(0),T(T.extent(firstDim)-1),0,cv_strikes(i)));
	  cvlist.push_back(callpayoff); }
	// Test list
	Array<double,1> test_underlying(2),test_numeraire(2);
    test_underlying = 1000.0,1000.0;
	test_numeraire  = 1.0,1.0;
	Array<double,1> test_strikes(cvlist.payoffArray(test_underlying,test_numeraire));
	std::cout << "Testing strikes: " << test_strikes << std::endl;
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_cv2(cvlist,ggbm,ts,numeraire_index);
	boost::function<Array<double,1> (Array<double,2>)> func_cv = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mappingArray,&mc_mapping_cv2,_1);
	MCGeneric<Array<double,2>,Array<double,1>,RandomArray<ranlib::NormalUnit<double>,double> > mccv2(func_cv,random_container_g);
	Array<int,1> CVidx(5);
	CVidx = 1, 2, 3, 4, 5;
	Array<double,2> CV_weights(5,1);
	n = minpaths;
	StringForm mystringform(12);
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound,";
	cout << "5 CV MC value,Optimal CV weight 1,Optimal CV weight 2\n";
	while (cvgatherer2.number_of_simulations()<maxpaths) {
	  mccv2.simulate(cvgatherer2,n);
	  cout << cvgatherer2.number_of_simulations() << ',' << mystringform(cvgatherer2.mean(0)) << ',';
	  cout << mystringform(cvgatherer2.mean(0)-d*cvgatherer2.stddev(0)) << ',';
	  cout << mystringform(cvgatherer2.mean(0)+d*cvgatherer2.stddev(0)) << ',';
	  cout << mystringform(cvgatherer2.CVestimate(0,CVidx,cv_values)) << ',';
      CV_weights = cvgatherer2.CVweight(0,CVidx);
	  for (i=0;i<5;i++) cout << mystringform(CV_weights(i,0)) << ','; 
	  cout << endl;
	  n = cvgatherer2.number_of_simulations(); }
	cvgatherer2.fix_weights(0,CVidx,cv_values);
	n = minpaths;
	cout << "Paths,Fixed weight CV MC value,95% CI lower bound,95% CI upper bound\n";
	while (cvgatherer2.number_of_simulations()<maxpaths) {
	  mccv2.simulate(cvgatherer2,n);
	  cout << cvgatherer2.number_of_simulations() << ',' << mystringform(cvgatherer2.mean(0)) << ',';
	  cout << mystringform(cvgatherer2.mean(0)-d*cvgatherer2.stddev(0)) << ',';
	  cout << mystringform(cvgatherer2.mean(0)+d*cvgatherer2.stddev(0)) << endl;
	  n = cvgatherer2.number_of_simulations(); }
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
