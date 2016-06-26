/** \file  ESXLWProject.cpp
    \brief C++ source file implementing the example in Appendix A.1 of the book Quantitative Finance: An Object-Oriented Approach in C++.
           Copyright 2012, 2014 by Erik Schlögl

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

#pragma warning (disable : 4996)

#include "ESXLWProject.hpp"
#include "TSBootstrap.hpp"
#include "BlackScholesAsset.hpp"
#include "Payoff.hpp"
#include "Binomial.hpp"
#include <boost/bind.hpp>
#include <ctime>
#include <windows.h>

std::string // tests empty args
EriksEmptyArgFunction()
{
    return "this function is useless except for testing.";
}

MyArray // echoes an array
ErikEchoArray(const MyArray& Echoee // argument to be echoed
                  )
{
    return Echoee;
}

MyArray // log-linear interpolation given set of discount factors
ESLogLinear(const MyArray& given_times // timeline for supplied discount factors
			,const MyArray& given_discount_factors  // supplied discount factors
			,const MyArray& target_times // times for which discount factors are desired
			)
{
	int i;
	int n = given_times.size();
	if (n!=given_discount_factors.size()) throw("Timeline and discount factor arrays must have same dimension");
	blitz::Array<double,1> T(n),B(n);
	MyArray result(target_times.size());
	for (i=0;i<n;i++) {
		T(i) = given_times[i];
		B(i) = given_discount_factors[i]; }
	quantfin::TSLogLinear ts(T,B);
	for (i=0;i<target_times.size();i++) result[i] = ts(target_times[i]);
    return result;
}

MyArray // log-linear interpolation given set of discount factors
ESLogLinearMatrix(const MyMatrix& given // timeline for supplied discount factors (first column) and supplied discount factors (second column)
			,const MyArray& target_times // times for which discount factors are desired
			)
{
	int i;
	int n = given.rows();
	blitz::Array<double,1> T(n),B(n);
	MyArray result(target_times.size());
	for (i=0;i<n;i++) {
		T(i) = given[i][0];
		B(i) = given[i][1]; }
	quantfin::TSLogLinear ts(T,B);
	for (i=0;i<target_times.size();i++) result[i] = ts(target_times[i]);
    return result;
}

double // Black/Scholes European option price
ESBlackScholes(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int callput // 1 for call, -1 for put
			   )
{
    quantfin::ConstVol vol(sigma);
    quantfin::BlackScholesAsset stock(&vol,S);
    return stock.option(T,K,r,callput);
}

double // Cox/Ross/Rubinstein European option price
ESCoxRossRubinstein(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int callput // 1 for call, -1 for put
			   ,int N // number of steps in the binomial lattice
			   )
{
    quantfin::ConstVol vol(sigma);
    quantfin::BlackScholesAsset stock(&vol,S);
    quantfin::BinomialLattice btree(stock,r,T,N+1);  
    quantfin::Payoff payoff(K,callput);
	// call payoff functor
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&quantfin::Payoff::operator()),&payoff,_1);
	// apply payoff functor to final period
    btree.apply_payoff(N,f);
	// roll back through lattice to time zero
    btree.rollback(N,0);
    return btree.result();
}

double // Cox/Ross/Rubinstein American put option price
ESCoxRossRubinsteinAmerican(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int N // number of steps in the binomial lattice
			   )
{
    quantfin::ConstVol vol(sigma);
    quantfin::BlackScholesAsset stock(&vol,S);
    quantfin::BinomialLattice btree(stock,r,T,N+1);  
    quantfin::Payoff payoff(K,-1);
	// call payoff functor
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&quantfin::Payoff::operator()),&payoff,_1);
    // American put option pricing
	quantfin::EarlyExercise amput(payoff);
	// Early exercise functor
    boost::function<double (double,double)> g;
    g = boost::bind(boost::mem_fn(&quantfin::EarlyExercise::operator()),&amput,_1,_2);
	// apply payoff functor to final period
    btree.apply_payoff(N,f);
	// roll back through lattice to time zero
    btree.rollback(N,0,g);
    return btree.result();
}
