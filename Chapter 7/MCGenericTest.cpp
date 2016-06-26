/** \file  MCGenericTest.cpp
    \brief Test and demonstration program.
           Copyright 2007, 2013 by Erik Schlögl

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
#include <boost/math/distributions/normal.hpp>
#include <random/normal.h>
#include "MCGeneric.hpp"
#include "BlackScholesAsset.hpp"

using namespace quantfin;
 
inline double positivePart(double x)
{
  return (x>0.0) ? x : 0.0;      
}

class objective_class {
private:
  double initial_stock_price;
  double maturity;
  double interest_rate;
  double volatility;
  double strike;
  double discount_factor;
  double sqrt_maturity;
public:
  inline objective_class(double S,double T,double r,double sgm,double K)
	: initial_stock_price(S),maturity(T),interest_rate(r),volatility(sgm),strike(K),discount_factor(std::exp(-r*T)),sqrt_maturity(std::sqrt(maturity)) { };
  double func(double x);
};

double objective_class::func(double x)
{
  return discount_factor * positivePart(initial_stock_price*std::exp((interest_rate-0.5*volatility*volatility)*maturity+volatility*sqrt_maturity*x)-strike);
}

/** Test and demonstration for finite difference schemes.

    Command-line arguments:
      -# initial stock price.
         The default is 100.
      -# interest rate.
         The default is 5%.
      -# volatility.
         The default is 30%. If volatility is not constant, this is the volatility for the
         first time segment.
      -# maturity.
         The default is 1.5.
      -# moneyness.
         The default is 1 (at the money).
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
    size_t minpaths = 100;
    if (argc>6) minpaths = atoi(argv[6]);
    size_t maxpaths = 100;
    if (argc>7) maxpaths = atoi(argv[7]);

    ConstVol vol(sgm);
    BlackScholesAsset stock(&vol,S);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl << endl;
    double CFcall = stock.option(mat,K,r);
	ranlib::NormalUnit<double> normal_RNG;
	objective_class obj(S,mat,r,sgm,K);
	boost::function<double (double)> func = boost::bind(&objective_class::func,&obj,_1);
	MCGeneric<double,double,ranlib::NormalUnit<double> > mc(func,normal_RNG);
	MCGatherer<double> mcgatherer;
	unsigned long n = minpaths;
	cout << "Paths,Black/Scholes value,MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors\n";
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  cout << mcgatherer.number_of_simulations() << ',' << CFcall << ',' << mcgatherer.mean() << ',' << mcgatherer.mean()-d*mcgatherer.stddev() << ',' << mcgatherer.mean()+d*mcgatherer.stddev() << ',' << (mcgatherer.mean()-CFcall)/mcgatherer.stddev() << endl;
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
