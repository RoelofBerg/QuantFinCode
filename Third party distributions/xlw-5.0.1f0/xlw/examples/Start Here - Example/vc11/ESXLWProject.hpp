/** \file  ESXLWProject.hpp
    \brief C++ header file implementing the example in Appendix A.1 of the book Quantitative Finance: An Object-Oriented Approach in C++.
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

#ifndef TEST_H
#define TEST_H

#include <xlw/MyContainers.h>
#include <xlw/CellMatrix.h>
#include <xlw/DoubleOrNothing.h>
#include <xlw/ArgList.h>
   
using namespace xlw;
  
//<xlw:libraryname=EriksTestLibrary
std::string // tests empty args
EriksEmptyArgFunction();

MyArray // echoes an array
ErikEchoArray(const MyArray& Echoee // argument to be echoed
                  );

MyArray // log-linear interpolation given set of discount factors
ESLogLinear(const MyArray& given_times // timeline for supplied discount factors
			,const MyArray& given_discount_factors  // supplied discount factors
			,const MyArray& target_times // times for which discount factors are desired
			);

MyArray // log-linear interpolation given set of discount factors
ESLogLinearMatrix(const MyMatrix& given // timeline for supplied discount factors (first column) and supplied discount factors (second column)
			,const MyArray& target_times // times for which discount factors are desired
			);

double // Black/Scholes European option price
ESBlackScholes(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int callput // 1 for call, -1 for put
			   );

double // Cox/Ross/Rubinstein European option price
ESCoxRossRubinstein(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int callput // 1 for call, -1 for put
			   ,int N // number of steps in the binomial lattice
			   );

double // Cox/Ross/Rubinstein American put option price
ESCoxRossRubinsteinAmerican(double S // current price of underlying asset
			   ,double K // exercise price
			   ,double sigma // volatility of underlying asset
			   ,double r // risk-free interest rate (continuously compounded)
			   ,double T // time to maturity (in years)
			   ,int N // number of steps in the binomial lattice
			   );


#endif
