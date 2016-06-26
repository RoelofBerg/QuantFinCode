/** \file  HestonTest.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2007, 2008, 2009, 2015 by Erik Schlögl

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
#include <random/uniform.h>
#include <random/normal.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/normal.hpp>
#include "BlackScholesAsset.hpp"
#include "PiecewiseVol.hpp"
#include "GeometricBrownianMotion.hpp"
#include "Payoff.hpp"
#include "StringForm.hpp"
#include "MCEngine.hpp"
#include "MCGeneric.hpp"
#include "QFRandom.hpp"
#include "HestonAsset.hpp"

using namespace quantfin;

/** Test and demonstration for Heston (1993) model.

    Command-line arguments:
      -# initial stock price.
         The default is 100.
      -# interest rate.
         The default is 5%.
      -# volatility of volatility.
         The default is 10%.
      -# maturity.
         The default is 1.5.
      -# moneyness.
         The default is 1 (at the money).
      -# kappa (speed of mean reversion for volatility process).
         The default is 2.
      -# theta (long term mean volatility).
         The default is 0.01.
      -# initial volatility level.
         The default is 0.01.
      -# rho.
         The default is 0.
      -# lambda.
         The default is 0.
      -# Gauss/Laguerre quadrature refinement.
         The default is 10.
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
    double sgm = 0.1;
    if (argc>3) sgm = atof(argv[3]);
    double mat = 1.5;
    if (argc>4) mat = atof(argv[4]);
    double K = 1.0;
    if (argc>5) K = atof(argv[5]);
    K *= S;
    double kappa = 2.0;
    if (argc>6) kappa = atof(argv[6]);
    double theta = 0.01;
    if (argc>7) theta = atof(argv[7]);
    double Vzero = 0.01;
    if (argc>8) Vzero = atof(argv[8]);
    double rho = 0.0;
    if (argc>9) rho = atof(argv[9]);
    double lambda = 0.0;
    if (argc>10) lambda = atof(argv[10]);
    int gausslaguerre_n = 10;
    if (argc>11) gausslaguerre_n = atoi(argv[11]);
    int N = 10;
    if (argc>12) N = atoi(argv[12]);
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
    Array<double,1> sgm1(1);
    sgm1 = sgm;
	CIRprocess vol_process(kappa,theta,Vzero,sgm1);
    HestonAsset stock(vol_process,S,rho,lambda,gausslaguerre_n);
	cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nvolvol: " << sgm << "\nrho: " << rho << endl;
	cout << "kappa (mean reversion of vol process): " << kappa << endl;
	cout << "theta (long run vol mean): " << theta << endl;
	cout << "Initial vol: " << Vzero << "\nlambda: " << lambda << endl;
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
