/** \file  FDTest.cpp
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
#include "BlackScholesAsset.hpp"
#include "FiniteDifference.hpp"
#include "MCGatherer.hpp"
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
      -# state space refinement.
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
    int Nj = 10;
    if (argc>7) Nj = atoi(argv[7]);
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
	cout << "Closed form call: " << CFcall << "  delta: " << stock.delta(mat,K,r) << "  gamma: " << stock.gamma(mat,K,r) << endl;
    cout << "Closed form put: " << CFput << "  delta: " << stock.delta(mat,K,r,-1) << "  gamma: " << stock.gamma(mat,K,r,-1) << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
    cout << "Time line refinement: " << N << "\nState space refinement: " << Nj << endl;
    cout << "Creating FiniteDifference object" << endl;
    FiniteDifference fd(stock,r,mat,N,Nj);
    cout << "Creating CrankNicolson object" << endl;
    CrankNicolson cn(stock,r,mat,N,Nj);  
    Payoff call(K);
    Payoff put(K,-1);
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&call,_1);
    fd.apply_payoff(N-1,f);
    fd.rollback(N-1,0);
	cout << "FD call: " << fd.result() << "\nDifference to closed form: " << CFcall - fd.result() << endl;
	cout << "FD call delta: " << fd.delta() << "  FD call gamma: " << fd.gamma() << endl;
    cn.apply_payoff(N-1,f);
    cn.rollback(N-1,0);
    cout << "CrankNicolson call: " << cn.result() << "\nDifference to closed form: " << CFcall - cn.result() << endl;
	cout << "CrankNicolson call delta: " << cn.delta() << "  CrankNicolson call gamma: " << cn.gamma() << endl;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&put,_1);
    fd.apply_payoff(N-1,f);
    fd.rollback(N-1,0);
    cout << "FD put: " << fd.result() << "\nDifference to closed form: " << CFput - fd.result() << endl;
	cout << "FD put delta: " << fd.delta() << "  FD put gamma: " << fd.gamma() << endl;
    cn.apply_payoff(N-1,f);
    cn.rollback(N-1,0);
    cout << "CrankNicolson put: " << cn.result() << "\nDifference to closed form: " << CFput - cn.result() << endl;
	cout << "CrankNicolson put delta: " << cn.delta() << "  CrankNicolson put gamma: " << cn.gamma() << endl;
    EarlyExercise amput(put);
    boost::function<double (double,double)> g;
    g = boost::bind(boost::mem_fn(&EarlyExercise::operator()),&amput,_1,_2);
    fd.apply_payoff(N-1,f);
    fd.rollback(N-1,0,g);
    double fdamput = fd.result();
    cout << "FD American put: " << fdamput << endl;
	cout << "FD American put delta: " << fd.delta() << "  \nFD American put gamma: " << fd.gamma() << endl;
    cn.apply_payoff(N-1,f);
    cn.rollback(N-1,0,g);
    cout << "CrankNicolson American put: " << cn.result() << endl;
	cout << "CrankNicolson American put delta: " << cn.delta() << "  \nCrankNicolson American put gamma: " << cn.gamma() << endl;

    // with dividends
    stock.dividend_yield(0.03);
    FiniteDifference fd2(stock,r,mat,N,Nj);
    CFcall = stock.option(mat,K,r);
    CFput  = stock.option(mat,K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&call,_1);
    fd2.apply_payoff(N-1,f);
    fd2.rollback(N-1,0);
    cout << "FD call: " << fd2.result() << endl;
    CrankNicolson cn2(stock,r,mat,N,Nj);  
    cn2.apply_payoff(N-1,f);
    cn2.rollback(N-1,0);
    cout << "CrankNicolson call: " << cn2.result() << endl;
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
