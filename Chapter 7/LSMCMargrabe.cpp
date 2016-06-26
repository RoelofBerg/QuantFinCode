/** \file  LSMCMargrabe.cpp
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
 
/** Test and demonstration for American Margrabe option by Monte Carlo.

    Command-line arguments:
      -# initial stock price 1.
         The default is 100.
      -# initial stock price 2.
         The default is 100.
      -# interest rate.
         The default is 5%.
      -# volatility 1.
         The default is 30%.
      -# volatility 2.
         The default is 30%.
	  -# correlation.
	     The default is 0.0.
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
      -# Switch for including a European option as one of the basis functions.
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
    double S1 = 100.0;
    if (argc>1) S1 = atof(argv[1]);
    double S2 = 100.0;
    if (argc>2) S2 = atof(argv[2]);
    double r = 0.05;
    if (argc>3) r = atof(argv[3]);
    double sgm1 = 0.3;
    if (argc>4) sgm1 = atof(argv[4]);
    double sgm2 = 0.3;
    if (argc>5) sgm2 = atof(argv[5]);
    double rho = 0.0;
    if (argc>6) rho = atof(argv[6]);
    double mat = 1.5;
    if (argc>7) mat = atof(argv[7]);
    double K = 1.0;
    if (argc>8) K = atof(argv[8]);
    K *= S1/S2;
    int N = 10;
    if (argc>9) N = atoi(argv[9]);
    size_t minpaths = 100;
    if (argc>10) minpaths = atoi(argv[10]);
    size_t train = 100;
    if (argc>11) train = atoi(argv[11]);
    int degree = 2;
    if (argc>12) degree = atoi(argv[12]);
    size_t maxpaths = 100;
    if (argc>13) maxpaths = atoi(argv[13]);
    bool include_put = false;
    if (argc>14) include_put = atoi(argv[14]);
	int numeraire_index = -1;
    double dividend1 = 0.0;
    if (argc>15) dividend1 = atof(argv[15]);
    double dividend2 = 0.0;
    if (argc>16) dividend2 = atof(argv[16]);
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
	Array<double,1> sgm_1(2),sgm_2(2);
	sgm_1 = sgm1, 0.0;
	sgm_2 = rho, std::sqrt(1-rho*rho);
	sgm_2 *= sgm2;
    ConstVol vol1(sgm_1);
    ConstVol vol2(sgm_2);
    BlackScholesAsset stock1(&vol1,S1,dividend1);
    BlackScholesAsset stock2(&vol2,S2,dividend2);
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock1);
    underlying.push_back(&stock2);
    FlatTermStructure ts(r,0.0,mat+10.0);
    double CFcall = stock1.Margrabe(stock2,mat,K);
    cout << "Closed form price: " << CFcall << endl;
	exotics::Margrabe Mopt(stock1,stock2,T(0),mat,1.0,K,ts);
    cout << "Closed form price via MBinary: " << Mopt.price() << endl;
    BlackScholesAsset tmpstock1(&vol1,S1*1.5,dividend1);
    BlackScholesAsset tmpstock2(&vol2,S2*1.6,dividend2);
	double tmp_t = T(N/2);
	// history as timepoints x assets
	Array<double,2> S_history(2,N/2+1);
	// in the case of the Margrabe option only the last timepoint in the history matters
	S_history(0,N/2) = S1*1.5;
	S_history(1,N/2) = S2*1.6;
	exotics::Margrabe tmpMopt(tmpstock1,tmpstock2,tmp_t,mat,1.0,K,ts);
	cout << "Testing: " << tmpMopt.price() << " == ";
	cout << Mopt.price(T(Range(fromStart,N/2)),S_history) << endl;
	cout << "Dividend 1: " << dividend1 << "   Dividend 2: " << dividend2 << endl;
	// American option by Monte Carlo
    GeometricBrownianMotion gbm(underlying);
	gbm.set_timeline(T);
    ranlib::NormalUnit<double> normalRNG;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,gbm.factors(),gbm.number_of_steps()); 
	MCTrainingPaths<GeometricBrownianMotion,RandomArray<ranlib::NormalUnit<double>,double> >
	  training_paths(gbm,T,train,random_container,ts,numeraire_index);
	cout << "Training paths created." << endl;
	// payoff requires (time points) x (state variables) Array as second argument
	boost::function<double (const Array<double,1>&,const Array<double,2>&)> payoff = boost::bind(&exotics::Margrabe::early_exercise_payoff,&Mopt,_1,_2);
	std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> > basisfunctions;
	Array<int,1> p(2);
	p(1) = 0.0;
	for (i=0;i<=degree;i++) {
	  p(0) = i;
	  add_polynomial_basis_function(basisfunctions,p); }
	p(0) = 0.0;
	for (i=1;i<=degree;i++) {
	  p(1) = i;
	  add_polynomial_basis_function(basisfunctions,p); }
    boost::function<double (const Array<double,1>&,const Array<double,2>&)> put_option;
	// basis functions require (state variables) X (observation time points) Array as second argument
    put_option = boost::bind(&exotics::Margrabe::price,&Mopt,_1,_2);
	if (include_put) basisfunctions.push_back(put_option);
	cout << "Fitting exercise boundary..." << endl;
	// training_paths is currently a paths x (time points) x (state variables) Array
	RegressionExerciseBoundary boundary(T,training_paths.state_variables(),training_paths.numeraires(),payoff,basisfunctions);
	cout << "Creating exercise strategy..." << endl;
	LSExerciseStrategy<RegressionExerciseBoundary> exercise_strategy(boundary);
	cout << "Setting up Monte Carlo simulation..." << endl;
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
