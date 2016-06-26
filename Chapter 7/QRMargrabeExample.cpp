/** \file  QRMargrabeExample.cpp
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
#include <random/uniform.h>
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
#include "QMCEngine.hpp"
#include "MCGeneric.hpp"
#include "QFRandom.hpp"
#include "QFQuasiRandom.hpp"


using namespace quantfin;

/** Test and demonstration for quasirandom Monte Carlo.

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
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
    FlatTermStructure ts(r,0.0,mat+10.0);
    Array<int,2> xS(1,1);
    xS = 1;
    Array<double,1> sgm2(2);
    sgm2 = 0.1, 0.4;
    ConstVol vol2(sgm2);
    BlackScholesAsset stock2(&vol2,S,0.04);
    K = 1.0;
	double margrabe_price = stock.Margrabe(stock2,mat,K);
    cout << "Closed form Margrabe option: " << margrabe_price << endl;
    underlying.push_back(&stock2);
    Array<int,2> mindex(2,2);
    mindex = 0, 1, 
             N, N;
    Array<double,1> malpha(2);
    malpha = 1.0, 0.0;
    Array<double,2> mA(1,2);
    mA = 1.0, -1.0;
    Array<double,1> ma(1);
    ma = K;
    MBinary M1(underlying,ts,T,mindex,malpha,xS,mA,ma);
	double M1price = M1.price(1000000000);
    Array<double,1> malpha2(2);
    malpha2 = 0.0, 1.0;
    MBinary M2(underlying,ts,T,mindex,malpha2,xS,mA,ma);
    cout << "MBinary price of Margrabe option: " << M1price-K*M2.price(1000000000) << std::endl;
    exotics::Margrabe margrabe(stock,stock2,0.0,mat,1.0,K,ts);
    cout << "Price of Margrabe option via MExotics: " << margrabe.price() << std::endl;

	// Margrabe option price by Monte Carlo
    ranlib::NormalUnit<double> normalRNG;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,2,N);  // 2 factors, N time steps
	SobolArrayNormal sobol(2,N,maxpaths*2);
    GeometricBrownianMotion gbm(underlying);
	boost::shared_ptr<MBinaryPayoff> M1payoff = M1.get_payoff();
	boost::shared_ptr<MBinaryPayoff> M2payoff = M2.get_payoff();
	MCPayoffList mcpayofflist;
    mcpayofflist.push_back(M1payoff);
    mcpayofflist.push_back(M2payoff,-K);
	cout << "Strike coefficient: " << K << endl;
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping(mcpayofflist,gbm,ts,numeraire_index);
	boost::function<double (Array<double,2>)> func = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mapping,&mc_mapping,_1);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc(func,random_container);
	MCGatherer<double> mcgatherer;
	// boost functor to convert random variates to their antithetics (instantiated from template)
	boost::function<Array<double,2> (Array<double,2>)> antithetic = normal_antithetic<Array<double,2> >;
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc_antithetic(func,random_container,antithetic);
	MCGatherer<double> mcgatherer_antithetic;
	MCGeneric<Array<double,2>,double,SobolArrayNormal> mc_QR(func,sobol);
	MCGatherer<double> mcgathererQR;
	// For nested construction of randomised QMC using random shift
	MCGatherer<double> mcgathererQRran;
	ranlib::Uniform<double> ugen;
	RandomArray<ranlib::Uniform<double>,double> unigen(ugen,2,N);

	unsigned long n = minpaths;
	cout << "Minimum number of paths: " << minpaths << "\nMaximum number of paths: " << maxpaths << endl;
	cout << "Margrabe option\nPaths,Closed form value,MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors,";
	cout << "Antithetic MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors,CI width,Antithetic CI width,";
	cout << "QR MC value,randomised QR MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors,CI width\n";
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	double CFprice = margrabe.price();
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << CFprice << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev() << ',' << (mcgatherer.mean()-CFprice)/mcgatherer.stddev() << ',' << std::flush;
	  mc_antithetic.simulate(mcgatherer_antithetic,n/2);
	  cout << mcgatherer_antithetic.mean() << ',' << mcgatherer_antithetic.mean()-d*mcgatherer_antithetic.stddev() << ',' << mcgatherer_antithetic.mean()+d*mcgatherer_antithetic.stddev() << ',' << (mcgatherer_antithetic.mean()-CFprice)/mcgatherer_antithetic.stddev() << ',' << std::flush;
	  cout << 2.0*d*mcgatherer.stddev() << ',' << 2.0*d*mcgatherer_antithetic.stddev() << ',';
      mc_QR.simulate(mcgathererQR,n-1);
	  cout << mcgathererQR.mean() << ',';
	  // Nested construction of randomised QMC using random shift
	  SobolArrayNormal sobolr(2,N,n/32);
	  RandomShiftQMCMapping<GeometricBrownianMotion,Array<double,2>,SobolArrayNormal> QR_mapping(sobolr,mcpayofflist,gbm,ts,numeraire_index);
	  boost::function<double (Array<double,2>)> QRfunc = boost::bind(&RandomShiftQMCMapping<GeometricBrownianMotion,Array<double,2>,SobolArrayNormal>::mapping,&QR_mapping,_1);
	  MCGeneric<Array<double,2>,double,RandomArray<ranlib::Uniform<double>,double> > randomQMC(QRfunc,unigen);
	  mcgathererQRran.reset();
      randomQMC.simulate(mcgathererQRran,32);
	  cout << mcgathererQRran.mean() << ','<< mcgathererQRran.mean()-d*mcgathererQRran.stddev() << ',' << mcgathererQRran.mean()+d*mcgathererQRran.stddev() << ',' << (mcgathererQRran.mean()-CFprice)/mcgathererQRran.stddev() << ',' << 2.0*d*mcgathererQRran.stddev() << endl;
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
