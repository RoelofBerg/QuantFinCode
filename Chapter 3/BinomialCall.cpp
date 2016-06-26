/** \file  BinomialCall.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2010, 2013 by Erik Schlögl

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
#include "BlackScholesAsset.hpp"
#include "Payoff.hpp"
#include "Binomial.hpp"

using namespace quantfin;
 
/** Test and demonstration for binomial lattices.

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
    double dt = mat/N;
    ConstVol vol(sgm);
    BlackScholesAsset stock(&vol,S);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl;
    double CFcall = stock.option(mat,K,r);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Time line refinement: " << N << endl;
    cout << "Creating BinomialLattice object" << endl;
    BinomialLattice btree(stock,r,mat,N);  
    Payoff call(K);
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&call,_1);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (CRR): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
    btree.set_JarrowRudd();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (JR): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
    btree.set_Tian();
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (Tian): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
    btree.set_LeisenReimer(K);
    btree.apply_payoff(N-1,f);
    btree.rollback(N-1,0);
    cout << "Binomial call (LR): " << btree.result() << "\nDifference to closed form: " << CFcall - btree.result() << endl;
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
