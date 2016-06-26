/** \file  MCBondOpt.cpp
    \brief Test and demonstration program.
           Copyright 2011, 2012 by Erik Schlögl

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
#include <boost/bind.hpp>
#include "MCEngine.hpp"
#include "QFRandom.hpp"
#include "MCGeneric.hpp"
#include "GaussianEconomy.hpp"
#include "ExponentialVol.hpp"

using namespace quantfin;
 
/** Test and demonstration for HJM zero coupon bond options by Monte Carlo.

    Command-line arguments:
      -# underlying bond maturity.
         The default is 10.
      -# interest rate.
         The default is 5%.
      -# volatility.
         The default is 3%.
      -# maturity.
         The default is 1.5.
      -# moneyness.
         The default is 1 (at the money).
      -# N, time line refinement.
         The default is 10.
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Maximum number of MC pricing paths.
         The default is 100.
      -# Mean reversion parameter.
         The default is 0.1.

	Usage e.g. mcbondopt 10 0.05 0.03 1.5 1 3 100 1000000000 0.1
*/
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  try {
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
    double bondmat = 10.0;
    if (argc>1) bondmat = atof(argv[1]);
    double r = 0.05;
    if (argc>2) r = atof(argv[2]);
    double sgm = 0.03;
    if (argc>3) sgm = atof(argv[3]);
    double mat = 1.5;
    if (argc>4) mat = atof(argv[4]);
    double K = 1.0;
    if (argc>5) K = atof(argv[5]);
    int N = 10;
    if (argc>6) N = atoi(argv[6]);
    size_t minpaths = 100;
    if (argc>7) minpaths = atoi(argv[7]);
    size_t maxpaths = 100;
    if (argc>8) maxpaths = atoi(argv[8]);
    double mean_reversion = 0.1;
    if (argc>9) mean_reversion = atof(argv[9]);
	cout << bondmat << ' ' << r << ' ' << sgm << ' ' << mat << ' ' << K << ' ' << N << ' ' << minpaths << ' ' << maxpaths;
	cout << ' ' << mean_reversion << endl;
	int numeraire_index = 0;
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
	boost::shared_ptr<DeterministicAssetVol> hjmvol(new ExponentialVol(sgm,mean_reversion));
    std::vector<boost::shared_ptr<BlackScholesAsset> > underlying;
    boost::shared_ptr<TermStructure> ts(new FlatTermStructure(r,0.0,bondmat+10.0));
	double strike = K * (*ts)(bondmat)/(*ts)(mat);
	boost::shared_ptr<GaussianEconomy> domestic_economy(new GaussianEconomy(underlying,hjmvol,ts));
	std::vector<boost::shared_ptr<GaussianEconomy> > economies;
	economies.push_back(domestic_economy);
    double CFcall = domestic_economy->hjm->ZCBoption(mat,bondmat,strike);
    cout << "Closed form call: " << CFcall << endl;
	std::vector<boost::shared_ptr<DeterministicAssetVol> > xvols;
	Array<double,1> initial_exchange_rates(0);
	GaussMarkovWorld world(economies,xvols,initial_exchange_rates);
	world.set_timeline(T);
	world.set_reporting(0,-1,bondmat);
    ranlib::NormalUnit<double> normalRNG;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,world.factors(),world.number_of_steps()); 
    MCEuropeanCall callpayoff(T(0),mat,0,strike);
	MCMapping<GaussMarkovWorld,Array<double,2> > mc_mapping(callpayoff,world,*ts,numeraire_index);
	boost::function<double (Array<double,2>)> func = boost::bind(&MCMapping<GaussMarkovWorld,Array<double,2> >::mapping,&mc_mapping,_1);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc(func,random_container);
	MCGatherer<double> mcgatherer;
	size_t n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	cout << "Number of time steps: " << world.number_of_steps() << endl;
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound,Std error" << endl;
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev();
	  cout << ',' << (mcgatherer.mean()-CFcall)/mcgatherer.stddev() << endl;
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
