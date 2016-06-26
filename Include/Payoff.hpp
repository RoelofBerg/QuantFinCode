/** \file  Payoff.hpp
    \brief Header file declaring classes to represent option payoffs.
           Copyright 2006 by Erik Schlögl

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

#ifndef PAYOFF_HPP
#define PAYOFF_HPP

namespace quantfin {

class Payoff {
private:
  double    K;
  int    sign;
public:
  inline Payoff(double xK,int xsign = 1) : K(xK),sign(xsign) { };
  double operator()(double S);
};

class EarlyExercise {
private:
  Payoff& payoff;
public:
  inline EarlyExercise(Payoff& xpayoff) : payoff(xpayoff) { };
  double operator()(double continuationValue,double S);
};

class KnockOut {
private:
  double   barrier;
  int    direction;
public:
  // Constructor. Down-and-out (direction = -1) is the default; setting direction = 1 results in an up-and-out option.
  inline KnockOut(double xbarrier,int xdirection = -1) : barrier(xbarrier),direction(xdirection) { };
  double operator()(double continuationValue,double S);
};


}

#endif
