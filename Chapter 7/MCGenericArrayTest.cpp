/** \file  main.cpp
    \brief Test and demonstration program.
           Copyright 2008, 2010 by Erik Schloegl 
     */

#include <iostream> 
#include <cstdlib>
#include <boost/math/distributions/normal.hpp>
#include <boost/bind.hpp>
#include <random/normal.h>
#include "MCGeneric.hpp"
#include "QFRandom.hpp"
//#include "MCGathererArray.hpp"
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
  Array<double,1> func(double x);
};

// Return discounted payoff of call and put
Array<double,1> objective_class::func(double x)
{
  Array<double,1> payoff(2);
  payoff(0) = discount_factor * positivePart(initial_stock_price*std::exp((interest_rate-0.5*volatility*volatility)*maturity+volatility*sqrt_maturity*x)-strike);
  payoff(1) = discount_factor * positivePart(strike-initial_stock_price*std::exp((interest_rate-0.5*volatility*volatility)*maturity+volatility*sqrt_maturity*x));
  return payoff;
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
      -# minimum MC relative accuracy.
         The default is 1e-2.
      -# maximum MC relative accuracy.
         The default is 1e-2.
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
    //K *= S*std::exp(r*mat);
    K *= S;
    double minacc = 1e-2;
    if (argc>6) minacc = atof(argv[6]);
    double maxacc = 1e-2;
    if (argc>7) maxacc = atof(argv[7]);

    ConstVol vol(sgm);
    BlackScholesAsset stock(&vol,S);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl << endl;
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
	ranlib::NormalUnit<double> normal_RNG;
	objective_class obj(S,mat,r,sgm,K);
	boost::function<Array<double,1> (double)> func = boost::bind(&objective_class::func,&obj,_1);
	boost::function<double (double)> antithetic = normal_antithetic<double>;
	MCGeneric<double,Array<double,1>,ranlib::NormalUnit<double> > mc(func,normal_RNG);
	MCGeneric<double,Array<double,1>,ranlib::NormalUnit<double> > mc_antithetic(func,normal_RNG,antithetic);
	MCGatherer<Array<double,1> > mcgatherer(2),mcgatherer_antithetic(2);
	double acc = minacc * CFcall;
	maxacc *= CFcall;
	cout << "Aiming for maximum accuracy of " << maxacc << endl;
	cout << "Method,Option,Required accuracy,Actual accuracy,Paths,Black/Scholes value,MC value,95% CI lower bound,95% CI upper bound,Difference in standard errors\n";
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	while (acc>=maxacc) {
	  double curr_acc = mc.simulate(mcgatherer,100,acc);
	  Array<double,1> mean(mcgatherer.mean());
	  Array<double,1> stddev(mcgatherer.stddev());
	  cout << "Standard,Call," << acc << ',' << curr_acc << ',' << mcgatherer.number_of_simulations() << ',' << CFcall << ',' << mean(0) << ',' << mean(0)-d*stddev(0) << ',' << mean(0)+d*stddev(0) << ',' << (mean(0)-CFcall)/stddev(0) << endl;
	  cout << "Standard,Put," << acc << ',' << curr_acc << ',' << mcgatherer.number_of_simulations() << ',' << CFput  << ',' << mean(1) << ',' << mean(1)-d*stddev(1) << ',' << mean(1)+d*stddev(1) << ',' << (mean(1)-CFput)/stddev(1) << endl;
	  double curr_acc_antithetic = mc_antithetic.simulate(mcgatherer_antithetic,100,acc);
	  Array<double,1> mean_antithetic(mcgatherer_antithetic.mean());
	  Array<double,1> stddev_antithetic(mcgatherer_antithetic.stddev());
	  cout << "Antithetic,Call," << acc << ',' << curr_acc_antithetic << ',' << mcgatherer_antithetic.number_of_simulations() << ',' << CFcall << ',' << mean_antithetic(0) << ',' << mean_antithetic(0)-d*stddev_antithetic(0) << ',' << mean_antithetic(0)+d*stddev_antithetic(0) << ',' << (mean_antithetic(0)-CFcall)/stddev_antithetic(0) << endl;
	  cout << "Antithetic,Put," << acc << ',' << curr_acc_antithetic << ',' << mcgatherer_antithetic.number_of_simulations() << ',' << CFput  << ',' << mean_antithetic(1) << ',' << mean_antithetic(1)-d*stddev_antithetic(1) << ',' << mean_antithetic(1)+d*stddev_antithetic(1) << ',' << (mean_antithetic(1)-CFput)/stddev_antithetic(1) << endl;
	  acc /= 2; }
    
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
