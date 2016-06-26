/** \file  DeterministicCashflow.hpp
    \brief Header file defining classes to represent deterministic cashflows.
           Copyright 2005 by Erik Schlögl

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

#ifndef DETERMINISTICCASHFLOW_HPP
#define DETERMINISTICCASHFLOW_HPP
#include "TermStructure.hpp"

namespace quantfin {

class DeterministicCashflow {
private:
  Array<double,1> timeline_;
  Array<double,1> payments;
  double value;
public:
  inline DeterministicCashflow() { };
  inline DeterministicCashflow(const DeterministicCashflow& cf) : timeline_(cf.timeline_.copy()),payments(cf.payments.copy()),value(cf.value) { /* std::cout << timeline_ << std::endl; */ };
  inline DeterministicCashflow(const Array<double,1>& t,const Array<double,1>& p,double val)
            : timeline_(t.copy()),payments(p.copy()),value(val) { };
  inline DeterministicCashflow& operator=(const DeterministicCashflow& cf) { timeline_ = cf.timeline_.copy(); payments = cf.payments.copy(); value = cf.value; return *this; };
  inline const Array<double,1>& timeline() const { return timeline_; };      
  inline const Array<double,1>& cashflow() const { return payments; };      
  inline double market_value() const { return value; };
  inline void set_market_value(double v) { value = v; };
  /// Calculate net present value of present and future (not past) cashflows.
  double NPV(const TermStructure& ts) const;
  /// Predicate for sorting cashflows by length
  class longer_cashflow {
  public:
    inline bool operator()(const DeterministicCashflow& first,const DeterministicCashflow& second) 
    { 
      return (first.timeline_(first.timeline_.extent(firstDim)-1)>second.timeline_(second.timeline_.extent(firstDim)-1));
    }
  };
  /// Predicate for sorting cashflows by length
  class shorter_cashflow {
  public:
    inline bool operator()(const DeterministicCashflow& first,const DeterministicCashflow& second) 
    { 
      return (first.timeline_(first.timeline_.extent(firstDim)-1)<second.timeline_(second.timeline_.extent(firstDim)-1));
    }
  };
};


}

#endif
