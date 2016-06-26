/** \file  LSMCHJM.cpp
    \brief Test and demonstration program.
           Copyright 2012, 2013 by Erik Schlögl

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
#include <fstream>
#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "TSBootstrap.hpp"
#include "QFRandom.hpp"
#include "CSV2Array.hpp"
#include "GaussianEconomy.hpp"
#include "MCGeneric.hpp"
#include "LongstaffSchwartz.hpp"
#include "MCAmerican.hpp"
#include "Payoff.hpp"
#include "MExotics.hpp"

using namespace quantfin;

/** Test and demonstration for Longstaff/Schwartz American option pricing in a Gauss/Markov HJM model.

    Command-line arguments:
	  -# option parameter CSV file name.
      -# GaussMarkovWorld CSV file name.
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Maximum number of MC pricing paths.
         The default is 100.

	Usage e.g. lsmchjm inputlsmchjm.csv worldbs.csv 100 4000000
  */
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i,j;
  try {
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
	if (argc<2) throw std::logic_error("Missing input parameter CSV file name");
	if (argc<3) throw std::logic_error("Missing GaussMarkovWorld CSV file name");
    size_t minpaths = 100;
    if (argc>3) minpaths = atoi(argv[3]);
    size_t maxpaths = 100;
    if (argc>4) maxpaths = atoi(argv[4]);
	// Read data from files and create multicurrency term structure model
    std::ifstream is_inputs(argv[1]);
    if (!is_inputs.is_open()) throw std::logic_error("Failed to open input parameter CSV file");
    blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
	GaussMarkovWorld world(argv[2]);
	double maturity = 10.0;
    size_t train = 100;
    int degree = 2;
	int IRdegree = 0;
	int crossdegree = 0;
	int steps = 60;
	double strike = 110.0;
	std::map<std::string,std::string> inputs_map;
	for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
    if (inputs_map.count("Strike"))                   strike = std::atof(inputs_map["Strike"].data());
    if (inputs_map.count("Expiry"))                 maturity = std::atof(inputs_map["Expiry"].data());
    if (inputs_map.count("Training paths"))            train = std::atoi(inputs_map["Training paths"].data());
    if (inputs_map.count("Polynomial degree"))        degree = std::atoi(inputs_map["Polynomial degree"].data());
    if (inputs_map.count("IR polynomial degree"))   IRdegree = std::atoi(inputs_map["IR polynomial degree"].data());
    if (inputs_map.count("Cross polynomial degree")) crossdegree = std::atoi(inputs_map["Cross polynomial degree"].data());
    if (inputs_map.count("Steps"))                     steps = std::atoi(inputs_map["Steps"].data());
	cout << "Training paths," << train << "\nPolynomial degree," << degree << "\nIR polynomial degree," << IRdegree << "\nCross polynomial degree," << crossdegree << "\nSteps," << steps << endl;
    Array<double,1> timeline(steps+1);
	double dt = maturity/steps;
    firstIndex idx;
    timeline = idx * dt;
	const std::vector<boost::shared_ptr<GaussianEconomy> >& ecv = world.get_economies();
	const GaussianEconomy& domestic_economy = *(ecv[0]);
	boost::shared_ptr<BlackScholesAsset> Sp(domestic_economy.underlying[0]);
	const BlackScholesAsset& S = *Sp;
	world.set_timeline(timeline);
	Array<double,1> numeraire_values(timeline.extent(firstDim));
	cout << "European option price: " << domestic_economy.hjm->option(S,maturity,strike,-1) << endl;

	boost::mt19937 boost_random_engine; // this should be the single instance of a random number generator in the whole system
	// Run the actual simulation
	int numeraire_index = 0;
	cout << "Numeraire index: " << numeraire_index << endl;
	Array<double,1> T(timeline);
	world.set_reporting(0,0);
	world.set_reporting(0,-1,world.time_horizon());
	boost::normal_distribution<double> normal;
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > normalRNG(boost_random_engine,normal);
	RandomWrapper<boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >,double> normalRNGwrap(normalRNG);
    RandomArray<RandomWrapper<boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >,double>,double> 
		random_container(normalRNGwrap,world.factors(),world.number_of_steps()); 
	MCTrainingPaths<GaussMarkovWorld,RandomArray<RandomWrapper<boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >,double>,double> >
	  training_paths(world,T,train,random_container,*(domestic_economy.initialTS),numeraire_index);
	std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> > basisfunctions;
	Array<int,1> p(2);
	p(1) = 0.0;
	for (i=0;i<=degree;i++) {
	  p(0) = i;
	  add_polynomial_basis_function(basisfunctions,p); }
	if (IRdegree) {
	  p(0) = 0.0;
	  for (i=1;i<=IRdegree;i++) {
	    p(1) = i;
		add_polynomial_basis_function(basisfunctions,p); }}
	if (crossdegree) {
	  for (i=1;i<=crossdegree;i++) {
	    p(1) = p(0) = i;
		add_polynomial_basis_function(basisfunctions,p); }}
	// BlackScholesAsset inputs are not used by Monte Carlo payoff from Mopt
	exotics::ProductOption Mopt(*(domestic_economy.underlying[0]),*(domestic_economy.underlying[0]),T(0),world.time_horizon(),strike,*(domestic_economy.initialTS),-1);
	boost::function<double (const Array<double,1>&,const Array<double,2>&)> payoff = boost::bind(&exotics::ProductOption::early_exercise_payoff,&Mopt,_1,_2);
	Array<double,1> tmp(training_paths.state_variables()(1,Range::all(),1));
	cout << "Training path for bond:\n" << tmp << endl;
	RegressionExerciseBoundary boundary(T,training_paths.state_variables(),training_paths.numeraires(),payoff,basisfunctions);
	LSExerciseStrategy<RegressionExerciseBoundary> exercise_strategy(boundary);
	MCMapping<GaussMarkovWorld,Array<double,2> > mc_mapping(exercise_strategy,world,*(domestic_economy.initialTS),numeraire_index);
	boost::function<double (Array<double,2>)> func = boost::bind(&MCMapping<GaussMarkovWorld,Array<double,2> >::mapping,&mc_mapping,_1);
	MCGeneric<Array<double,2>,
		      double,
			  RandomArray<RandomWrapper<boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >,double>,double> > 
	  mc(func,random_container);
	MCGatherer<double> mcgatherer;
	size_t n = minpaths;
	boost::math::normal normdist;
	double d = boost::math::quantile(normdist,0.95);
	cout << "Simulations,MC value,95% CI lower bound,95% CI upper bound" << endl;
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev();
	  cout << endl;
	  n = mcgatherer.number_of_simulations(); 
	  std::cerr << n << endl; }

  } // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
