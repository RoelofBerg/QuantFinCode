/** \file  CVTest.cpp
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
    BlackScholesAsset stock(&vol,S,0.03);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl;
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
    FlatTermStructure ts(r,0.0,mat+10.0);

	// Margrabe option price by Monte Carlo
    ranlib::NormalUnit<double> normalRNG;
	MCGatherer<double> mcgatherer;
	// boost functor to convert random variates to their antithetics (instantiated from template)
	boost::function<Array<double,2> (Array<double,2>)> antithetic = normal_antithetic<Array<double,2> >;
	MCGatherer<double> mcgatherer_antithetic;
	unsigned long n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);

	// Geometric average option price by Monte Carlo 
	Array<double,1> gsgm(1);
	gsgm = sgm;
	ConstVol gvol(gsgm);
    BlackScholesAsset gasset(&gvol,S);
	// Fixed strike
	exotics::DiscreteGeometricMeanFixedStrike geo(gasset,T,call_strike,ts,1,S);
	double CFgeo = geo.price();
	double CFgeoFixedStrike = CFgeo; // used as control variate for arithmetic average option below
	cout << "Closed form value of geometric average option: " << CFgeo << endl;
	boost::shared_ptr<MCPayoffList> geopayoff = geo.get_payoff();
    std::vector<const BlackScholesAsset*> g_underlying;
    g_underlying.push_back(&gasset);
    GeometricBrownianMotion ggbm(g_underlying);
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_g(*geopayoff,ggbm,ts,numeraire_index);
    RandomArray<ranlib::NormalUnit<double>,double> random_container_g(normalRNG,ggbm.factors(),ggbm.number_of_steps()); 

	// Arithmetic average option price by Monte Carlo - including using geometric mean as control variate
	MCDiscreteArithmeticMeanFixedStrike avgpayoff(0,T,call_strike,1,S);
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_avg(avgpayoff,ggbm,ts,numeraire_index);
	boost::function<double (Array<double,2>)> func_g = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mapping,&mc_mapping_avg,_1);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mcavg(func_g,random_container_g);
	mcgatherer.reset();
	// Build objects for antithetic simulation
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mcavg_antithetic(func_g,random_container_g,antithetic);
	mcgatherer_antithetic.reset();
	// Build objects for control variate simulation
	Array<double,1> cv_values(1);
	cv_values = CFgeoFixedStrike;
	MCControlVariateMapping<GeometricBrownianMotion,GeometricBrownianMotion,Array<double,2> > mc_cvmapping(mc_mapping_avg,mc_mapping_g,cv_values);
	MCGatherer<double> mcgathererCV;
	MCGatherer<double> mcgathererCV_antithetic;
	func_g = boost::bind(&MCControlVariateMapping<GeometricBrownianMotion,GeometricBrownianMotion,Array<double,2> >::mapping,&mc_cvmapping,_1);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mcavg_cv(func_g,random_container_g);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mcavg_cv_antithetic(func_g,random_container_g,antithetic);
	n = minpaths;
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound,Antithetic MC value,95% CI lower bound,95% CI upper bound,CI width,Antithetic CI width,";
	cout << "CV MC value,95% CI lower bound,95% CI upper bound,Antithetic CV MC value,95% CI lower bound,95% CI upper bound,CI width,Antithetic CI width\n";
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mcavg.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev() << ',';
	  mcavg_antithetic.simulate(mcgatherer_antithetic,n/2);
	  cout << mcgatherer_antithetic.mean() << ',' << mcgatherer_antithetic.mean()-d*mcgatherer_antithetic.stddev() << ',' << mcgatherer_antithetic.mean()+d*mcgatherer_antithetic.stddev() << ',';
	  cout << 2.0*d*mcgatherer.stddev() << ',' << 2.0*d*mcgatherer_antithetic.stddev() << ',' << std::flush;
	  mcavg_cv.simulate(mcgathererCV,n);
	  cout << mcgathererCV.mean() << ',' << mcgathererCV.mean()-d*mcgathererCV.stddev() << ',' << mcgathererCV.mean()+d*mcgathererCV.stddev() << ',';
	  mcavg_cv_antithetic.simulate(mcgathererCV_antithetic,n/2);
	  cout << mcgathererCV_antithetic.mean() << ',' << mcgathererCV_antithetic.mean()-d*mcgathererCV_antithetic.stddev() << ',' << mcgathererCV_antithetic.mean()+d*mcgathererCV_antithetic.stddev() << ',';
      cout << 2.0*d*mcgathererCV.stddev() << ',' << 2.0*d*mcgathererCV_antithetic.stddev() << ',' << endl;
	  n = mcgatherer.number_of_simulations(); }

	// using MCGatherer<Array<double,1> > specialisation
	MCGatherer<Array<double,1> > cvgatherer(2);
	MCPayoffList                        cvlist;
	cvgatherer.set_control_variate(true);
	boost::shared_ptr<MCDiscreteArithmeticMeanFixedStrike> pavgpayoff(new MCDiscreteArithmeticMeanFixedStrike(0,T,call_strike,1,S));
	cvlist.push_back(pavgpayoff);
	cvlist.push_back(geopayoff);
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_cv(cvlist,ggbm,ts,numeraire_index);
	boost::function<Array<double,1> (Array<double,2>)> func_cv = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mappingArray,&mc_mapping_cv,_1);
	MCGeneric<Array<double,2>,Array<double,1>,RandomArray<ranlib::NormalUnit<double>,double> > mccv(func_cv,random_container_g);
	n = minpaths;
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound,";
	cout << "CV MC value,95% CI lower bound,95% CI upper bound,Optimal CV weight\n";
	while (cvgatherer.number_of_simulations()<maxpaths) {
	  mccv.simulate(cvgatherer,n);
	  cout << cvgatherer.number_of_simulations() << ',' << cvgatherer.mean(0) << ',' << cvgatherer.mean(0)-d*cvgatherer.stddev(0) << ',' << cvgatherer.mean(0)+d*cvgatherer.stddev(0) << ',';
	  cout << cvgatherer.CVestimate(0,1,CFgeoFixedStrike) << ',' << cvgatherer.CVestimate(0,1,CFgeoFixedStrike)-d*cvgatherer.CVestimate_stddev(0,1) << ',' << cvgatherer.CVestimate(0,1,CFgeoFixedStrike)+d*cvgatherer.CVestimate_stddev(0,1) << ',';
	  cout << cvgatherer.CVweight(0,1) << endl; 
	  n = cvgatherer.number_of_simulations(); }

	// using MCGatherer<Array<double,1> > specialisation with two control variates
	MCGatherer<Array<double,1> > cvgatherer2(3);
	cvgatherer2.set_control_variate(true);
	boost::shared_ptr<MCEuropeanCall> callpayoff(new MCEuropeanCall(T(0),T(T.extent(firstDim)-1),0,call_strike));
	std::cerr << "European call parameters: " << T(0) << ',' << T(T.extent(firstDim)-1) << ',' << call_strike << endl;
	cvlist.push_back(callpayoff);
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping_cv2(cvlist,ggbm,ts,numeraire_index);
	func_cv = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mappingArray,&mc_mapping_cv2,_1);
	MCGeneric<Array<double,2>,Array<double,1>,RandomArray<ranlib::NormalUnit<double>,double> > mccv2(func_cv,random_container_g);
	Array<int,1> CVidx(2);
	CVidx = 1, 2;
	Array<double,1> CV_expectation(2);
	Array<double,2> CV_weights(2,1);
	CV_expectation = CFgeoFixedStrike, gasset.option(T(T.extent(firstDim)-1)-T(0),call_strike,r);
	std::cerr << "European call parameters: " << T(T.extent(firstDim)-1)-T(0) << ',' << call_strike << endl;
	n = minpaths;
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound,";
	cout << "CV MC value,95% CI lower bound,95% CI upper bound,Optimal CV weight,";
	cout << "2 CV MC value,Optimal CV weight 1,Optimal CV weight 2\n";
	while (cvgatherer2.number_of_simulations()<maxpaths) {
	  mccv2.simulate(cvgatherer2,n);
	  cout << cvgatherer2.number_of_simulations() << ',' << cvgatherer2.mean(0) << ',' << cvgatherer2.mean(0)-d*cvgatherer2.stddev(0) << ',' << cvgatherer2.mean(0)+d*cvgatherer2.stddev(0) << ',';
	  cout << cvgatherer2.CVestimate(0,1,CFgeoFixedStrike) << ',' << cvgatherer2.CVestimate(0,1,CFgeoFixedStrike)-d*cvgatherer2.CVestimate_stddev(0,1) << ',' << cvgatherer2.CVestimate(0,1,CFgeoFixedStrike)+d*cvgatherer2.CVestimate_stddev(0,1) << ',';
	  cout << cvgatherer2.CVweight(0,1) << ','; 
	  cout << cvgatherer2.CVestimate(0,CVidx,CV_expectation) << ',';
      CV_weights = cvgatherer2.CVweight(0,CVidx);
	  cout << CV_weights(0,0) << ',' << CV_weights(1,0) << endl; 
	  n = cvgatherer2.number_of_simulations(); }
	cvgatherer2.fix_weights(0,CVidx,CV_expectation);
	n = minpaths;
	cout << "Paths,Fixed weight CV MC value,95% CI lower bound,95% CI upper bound\n";
	while (cvgatherer2.number_of_simulations()<maxpaths) {
	  mccv2.simulate(cvgatherer2,n);
	  cout << cvgatherer2.number_of_simulations() << ',' << cvgatherer2.mean(0) << ',' << cvgatherer2.mean(0)-d*cvgatherer2.stddev(0) << ',' << cvgatherer2.mean(0)+d*cvgatherer2.stddev(0) << endl;
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
