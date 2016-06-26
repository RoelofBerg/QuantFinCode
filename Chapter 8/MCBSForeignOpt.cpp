/** \file  MCBSForeignOpt.cpp
    \brief Test and demonstration program.
           Copyright 2012 by Erik Schlögl

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
#include <cstdlib>
#include <boost/bind.hpp>
#include "CSV2Array.hpp"
#include "MCEngine.hpp"
#include "QFRandom.hpp"
#include "MCGeneric.hpp"
#include "GaussianEconomy.hpp"
#include "ExponentialVol.hpp"
#include "MBinaryPayoffs.hpp"

using namespace quantfin;
 
/** Test and demonstration for GBM/HJM options by Monte Carlo.

    Command-line arguments:
	  -# option parameter CSV file name.
      -# GaussMarkovWorld CSV file name.
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Maximum number of MC pricing paths.
         The default is 100.
	  -# Time line refinement.
         The default is 10.

	Usage e.g. MCBSForeignOpt.exe option.csv worldbsf.csv 100 100000000
*/
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  double mat,strike;
  try {
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
	if (argc<2) throw std::logic_error("Missing option parameter CSV file name");
	if (argc<3) throw std::logic_error("Missing GaussMarkovWorld CSV file name");
    size_t minpaths = 100;
    if (argc>3) minpaths = atoi(argv[3]);
    size_t maxpaths = 100;
    if (argc>4) maxpaths = atoi(argv[4]);
    int N = 10;
    if (argc>5) N = atoi(argv[5]);
	// Read data from files and create multicurrency term structure model
    std::ifstream is_inputs(argv[1]);
    if (!is_inputs.is_open()) throw std::logic_error("Failed to open option parameter CSV file");
    blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
    std::map<std::string,std::string> inputs_map;
    for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
    std::string str("Option expiry");
    if (inputs_map.count(str)) mat = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing option expiry"); 
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
	str = "Strike";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing strike"); 
	GaussMarkovWorld world(argv[2]);
	std::cerr << "Getting economy..." << endl;
	const std::vector<boost::shared_ptr<GaussianEconomy> >& ecv = world.get_economies();
	const GaussianEconomy& ec = *(ecv[1]);
	std::cerr << "Number of assets: " << ec.underlying.size() << endl;
	std::cerr << "Getting asset..." << endl;
	boost::shared_ptr<BlackScholesAsset> Sp(ec.underlying[0]);
	std::cerr << "...got pointer..." << endl;
	const BlackScholesAsset& S = *Sp;
	std::cerr << "...got asset." << endl;
	strike *= S.initial_value() / (*world.get_economies()[1]->initialTS)(mat);
	std::cerr << "Strike: " << strike << endl;
	double CFcall = world.get_economies()[1]->hjm->option(S,mat,strike);
	cout << "Closed form value: " << CFcall << endl;
	double fxspot = (world.get_initial_exchange_rates())(0);
	CFcall *= fxspot;
	cout << "Closed form value converted to domestic currency: " << CFcall << endl;
	BlackScholesAsset fx(&*(world.get_FXvolatilities())[0],fxspot,0.0);
	world.set_timeline(T);
	std::cerr << "Number of steps: " << world.number_of_steps() << endl;
	// Underlying is first asset in foreign economy
	world.set_reporting(1,0);
	// Report exchange rate
	world.set_reporting(1,-2,mat);
	int numeraire_index = 0;
    ranlib::NormalUnit<double> normalRNG;
	// Payoff is converted to domestic currency - thus need MBinary
	boost::shared_ptr<ForeignOption> callpayoffasset(new ForeignOption(S,fx,strike,0.0,mat,*world.get_economies()[1]->initialTS));
	boost::shared_ptr<ForeignOption> callpayoffstrike(new ForeignOption(S,fx,strike,0.0,mat,*world.get_economies()[1]->initialTS,false));
	MCPayoffList mclist;
	mclist.push_back(callpayoffasset);
	mclist.push_back(callpayoffstrike,-strike);
	MCMapping<GaussMarkovWorld,Array<double,2> > mc_mapping(mclist,world,*(world.get_economies()[0]->initialTS),numeraire_index);
	boost::function<double (Array<double,2>)> func = boost::bind(&MCMapping<GaussMarkovWorld,Array<double,2> >::mapping,&mc_mapping,_1);
	std::cerr << "Number of factors: " << world.factors() << endl;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,world.factors(),world.number_of_steps()); 
	MCGeneric<Array<double,2>,double,RandomArray<ranlib::NormalUnit<double>,double> > mc(func,random_container);
	MCGatherer<double> mcgatherer;
	size_t n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
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
