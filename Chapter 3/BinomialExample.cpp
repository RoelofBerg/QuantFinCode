/** \file  main.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2010, 2013 by Erik Schloegl 
     */

#include <iostream> 
#include <cstdlib>
#include <boost/bind.hpp>
#include "BlackScholesAsset.hpp"
#include "Payoff.hpp"
#include "Binomial.hpp"

using namespace quantfin;
 
/** Test and demonstration for finite difference schemes.

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
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
    ConstVol vol(sgm);
    BlackScholesAsset stock(&vol,S);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl;
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
    cout << "Time line refinement: " << N << endl;
    cout << "Creating BinomialLattice object" << endl;
    BinomialLattice btree(stock,r,mat,N);  
    Payoff call(K);
    Payoff put(K,-1);
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&call,_1);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (CRR): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (JR): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&put,_1);
    btree.set_CoxRossRubinstein();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (CRR): " << btree.result() << "\nDifference to closed form: " << CFput - btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (JR): " << btree.result() << "\nDifference to closed form: " << CFput - btree.result() << endl;
    btree.set_Tian();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (Tian): " << btree.result() << "\nDifference to closed form: " << CFput - btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (LR): " << btree.result() << "\nDifference to closed form: " << CFput - btree.result() << endl;
    EarlyExercise amput(put);
    boost::function<double (double,double)> g;
    g = boost::bind(boost::mem_fn(&EarlyExercise::operator()),&amput,_1,_2);
    btree.set_CoxRossRubinstein();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (CRR): " << btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (JR): " << btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (LR): " << btree.result() << endl;

    // with dividends
    stock.dividend_yield(0.03);
    CFcall = stock.option(mat,K,r);
    CFput  = stock.option(mat,K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&call,_1);
    btree.set_CoxRossRubinstein();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (CRR): " << btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (JR): " << btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (LR): " << btree.result() << endl;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&put,_1);
    btree.set_CoxRossRubinstein();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (CRR): " << btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (JR): " << btree.result() << endl;
    btree.set_Tian();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (Tian): " << btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial put (LR): " << btree.result() << endl;
    btree.set_CoxRossRubinstein();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (CRR): " << btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (JR): " << btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0,g);
    cout << "Binomial American put (LR): " << btree.result() << endl;
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
