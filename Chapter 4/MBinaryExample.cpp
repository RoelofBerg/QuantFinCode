/** \file  MBinaryExample.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2007, 2009, 2010, 2013 by Erik Schlögl

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
#include "BlackScholesAsset.hpp"
#include "PiecewiseVol.hpp"
#include "Payoff.hpp"
#include "StringForm.hpp"
#include "MBinary.hpp"
#include "MExotics.hpp"


using namespace quantfin;

/** Test and demonstration for closed-form solutions in a Black/Scholes-type model.

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
    // read command line parameters
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
	double call_strike = K;
    int N = 10;
    if (argc>6) N = atoi(argv[6]);
	// set up timeline - N is the number of time steps, mat is the maturity of the option
    Array<double,1> T(N+1);
    firstIndex idx;
    double dt = mat/N;
    T = idx*dt;
	// create underlying asset - note that volatility is a two-dimensional vector with equal entries
    Array<double,1> sgm1(2);
    sgm1 = sgm;
    ConstVol vol(sgm1);
    BlackScholesAsset stock(&vol,S,0.03);
    cout << "S: " << S << "\nK: " << K << "\nr: " << r << "\nT: " << mat << "\nsgm: " << sgm << endl;
	// closed-form vanilla option prices
    double CFcall = stock.option(mat,K,r);
    double CFput  = stock.option(mat,K,r,-1);
    double CFiput = stock.option(T(N/2),K,r,-1);
    cout << "Closed form call: " << CFcall << endl;
    cout << "Closed form put: " << CFput << endl;
    cout << "Closed form intermediate maturity put: " << CFiput << endl;
	// initialise vector of underlying assets
    std::vector<const BlackScholesAsset*> underlying;
    underlying.push_back(&stock);
	// flat term structure for discounting
    FlatTermStructure ts(r,0.0,mat+10.0);
	// create MBinaries for European call option
    Array<int,2> index(2,1);
    index = 0, N;
    Array<double,1> alpha(1);
    alpha = 1.0;
    Array<int,2> xS(1,1);
    xS = 1;
    Array<double,2> xA(1,1);
    xA = 1.0;
    Array<double,1> xa(1);
    xa = K;
    MBinary V1(underlying,ts,T,index,alpha,xS,xA,xa);
    Array<double,1> alpha2(1);
    alpha2 = 0.0;
    MBinary V2(underlying,ts,T,index,alpha2,xS,xA,xa);
    cout << "MBinary price of call option: " << V1.price(1000000000)-K*V2.price(1000000000) << std::endl;
	// create second asset
    Array<double,1> sgm2(2);
    sgm2 = 0.1, 0.4;
    ConstVol vol2(sgm2);
    BlackScholesAsset stock2(&vol2,S,0.04);
	// option to exchange one asset for another
	// closed form price via formula implemented as member function of BlackScholesAsset - K is strike factor
    K = 1.0;
	double margrabe_price = stock.Margrabe(stock2,mat,K);
    cout << "Closed form Margrabe option: " << margrabe_price << endl;
	// add second asset to vector of underlying assets
    underlying.push_back(&stock2);
	// create MBinaries for Margrabe option
    Array<int,2> mindex(2,2);
    mindex = 0, 1, 
             N, N;
    Array<double,1> malpha(2);
    malpha = 1.0, 0.0;
    Array<double,2> mA(1,2);
    mA = 1.0, -1.0;
    Array<double,1> ma(1);
    ma = K;
    MBinary M1(underlying,ts,T,mindex,malpha,xS,mA,ma);
	double M1price = M1.price(1000000000);
    Array<double,1> malpha2(2);
    malpha2 = 0.0, 1.0;
    MBinary M2(underlying,ts,T,mindex,malpha2,xS,mA,ma);
    cout << "MBinary price of Margrabe option: " << M1price-K*M2.price(1000000000) << std::endl;
	// use the Margrabe option class from MExotics.hpp
    exotics::Margrabe margrabe(stock,stock2,0.0,mat,1.0,K,ts);
    cout << "Price of Margrabe option via MExotics: " << margrabe.price() << std::endl;

	// European call option price (on asset with piecewise constant vol)
	Array<double,2> piecewise_vols(T.extent(firstDim)-1,2);
	for (i=0;i<piecewise_vols.extent(firstDim);i++) {
	  piecewise_vols(i,0) = sgm2(0) * (T(i+1)-T(0))/(T(T.extent(firstDim)-1)-T(0)) * 1.5;
	  piecewise_vols(i,1) = sgm2(1) * (T(i+1)-T(0))/(T(T.extent(firstDim)-1)-T(0)) * 1.5; }
	PiecewiseConstVol pvol(T,piecewise_vols);
	BlackScholesAsset passet(&pvol,S,0.03);
    double CFprice = passet.option(mat,call_strike,r);

	// Geometric average option price 
	Array<double,1> gsgm(1);
	gsgm = sgm;
	ConstVol gvol(gsgm);
    BlackScholesAsset gasset(&gvol,S);
	// Fixed strike
	exotics::DiscreteGeometricMeanFixedStrike geo(gasset,T,call_strike,ts,1,S);
	double CFgeo = geo.price();
	cout << "Closed form value of geometric average option: " << CFgeo << endl;
	// Floating strike
	exotics::DiscreteGeometricMeanFloatingStrike geofloat(gasset,T,call_strike/S,ts,1,S);
	CFgeo = geofloat.price();
	cout << "Closed form value of geometric average option (floating strike): " << CFgeo << endl;
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
