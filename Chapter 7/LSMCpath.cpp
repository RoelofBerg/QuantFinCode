/** \file  LSMCpath.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2010 by Erik Schlögl

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
#include "GeometricBrownianMotion.hpp"
#include "QFRandom.hpp"
#include "BlackScholesAsset.hpp"
#include "Payoff.hpp"
#include "LongstaffSchwartz.hpp"
#include "MCAmerican.hpp"
#include "MCGeneric.hpp"
#include "MExotics.hpp"

using namespace quantfin;
 
/** Test and demonstration for American options by Monte Carlo.

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
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Longstaff/Schwartz training paths.
         The default is 100.
      -# Longstaff/Schwartz polynomial degree.
         The default is 2.
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Switch for including a European put as one of the basis functions.
         The default is 0 (false).
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
    int N = 10;
    if (argc>6) N = atoi(argv[6]);
    size_t minpaths = 100;
    if (argc>7) minpaths = atoi(argv[7]);
    size_t train = 100;
    if (argc>8) train = atoi(argv[8]);
    int degree = 2;
    if (argc>9) degree = atoi(argv[9]);
    size_t maxpaths = 100;
    if (argc>10) maxpaths = atoi(argv[10]);
    bool include_put = false;
    if (argc>11) include_put = atoi(argv[11]);
	int numeraire_index = -1;
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
    ConstVol vol(sgm);
    BlackScholesAsset stock(&vol,S);
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
    FlatTermStructure ts(r,0.0,mat+10.0);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl;
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
	exotics::StandardOption Mput(stock,T(0),mat,K,ts,-1);
    cout << "Closed form put via MBinary: " << Mput.price() << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
    BlackScholesAsset tmpstock(&vol,S*1.5);
	double tmp_t = (T(0)+mat)/2.0;
	exotics::StandardOption tmpMput(tmpstock,tmp_t,mat,K,ts,-1);
	cout << "Testing: " << tmpMput.price() << " == ";
	cout << Mput.price(tmp_t,S*1.5) << endl;
    GeometricBrownianMotion gbm(underlying);
	gbm.set_timeline(T);
    ranlib::NormalUnit<double> normalRNG;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,gbm.factors(),gbm.number_of_steps()); 
	MCTrainingPaths<GeometricBrownianMotion,RandomArray<ranlib::NormalUnit<double>,double> >
	  training_paths(gbm,T,train,random_container,ts,numeraire_index);
    Payoff put(K,-1);
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&put,_1);
	boost::function<double (const Array<double,1>&,const Array<double,2>&)> payoff = boost::bind(REBAdapter,_1,_2,f,0);
	std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> > basisfunctions;
	Array<int,1> p(1);
	for (i=0;i<=degree;i++) {
	  p(0) = i;
	  add_polynomial_basis_function(basisfunctions,p); }
    boost::function<double (double,double)> put_option;
    put_option = boost::bind(&exotics::StandardOption::price,&Mput,_1,_2);
	boost::function<double (const Array<double,1>&,const Array<double,2>&)> put_option_basis_function = boost::bind(REBAdapterT,_1,_2,put_option,0);
	if (include_put) basisfunctions.push_back(put_option_basis_function);
	RegressionExerciseBoundary boundary(T,training_paths.state_variables(),training_paths.numeraires(),payoff,basisfunctions);
	LSExerciseStrategy<RegressionExerciseBoundary> exercise_strategy(boundary);
	MCMapping<GeometricBrownianMotion,Array<double,2> > mc_mapping(exercise_strategy,gbm,ts,numeraire_index);
	boost::function<double (Array<double,2>)> func = boost::bind(&MCMapping<GeometricBrownianMotion,Array<double,2> >::mapping,&mc_mapping,_1);
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc(func,random_container);
	MCGatherer<double> mcgatherer;
	size_t n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	cout << "Paths,MC value,95% CI lower bound,95% CI upper bound" << endl;
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev() << endl;
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
