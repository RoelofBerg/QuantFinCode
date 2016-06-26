/** \file  TSBootstrap.hpp
    \brief Header file defining classes for the term structure of interest rates.
           Copyright 2005, 2010, 2013 by Erik Schlögl

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

#ifndef TSBOOTSTRAP_HPP
#define TSBOOTSTRAP_HPP
#include "TermStructure.hpp"

namespace quantfin {
          
/// Abstract base class for term structures which can be bootstrapped from coupon bond prices or, equivalently, swap rates.
class TSBootstrap : public TermStructure {
private:
  class bootstrap_class {
  private:
    TSBootstrap* ts_;
    const Array<double,1>& t_;
    const Array<double,1>& p_;
    const blitz::Range& slice;
    size_t idx;
  public:
    inline bootstrap_class(const Array<double,1>& t,const Array<double,1>& p,const blitz::Range& s,size_t i,TSBootstrap* ts) : ts_(ts),t_(t),p_(p),slice(s),idx(i) { };
    double operator()(double t);            
  };
public:
  /// General constructor.
  inline TSBootstrap(const Array<double,1>& xT,             ///< Time line of maturities .
		             const Array<double,1>& xB              ///< Time T(0) forward zero coupon bond prices for these maturities - thus the first bond price is always 1 
                     ) : TermStructure(xT.copy(),xB.copy()) { };
  void bootstrap(std::vector<DeterministicCashflow> cashflows,double eps = 1E-12);
};

/// Term structure using loglinear interpolation of zero coupon bond prices.
class TSLogLinear : public TSBootstrap {
public:      
  /// General constructor.
  inline TSLogLinear(const Array<double,1>& xT,           ///< Time line of maturities .
  		             const Array<double,1>& xB            ///< Time T(0) forward zero coupon bond prices for these maturities - thus the first bond price is always 1 
                     ) : TSBootstrap(xT.copy(),xB.copy()) { };
  inline TSLogLinear(const TSLogLinear& ts) : TSBootstrap(ts.T.copy(),ts.B.copy()) { };
  /// Virtual copy constructor
  virtual TSLogLinear* pointer_to_copy() const;
  /** Function returning the interpolated value, which for term structures of
      interest rates is the maturity t, time T(0) forward zero coupon bond price. */
  virtual double operator()(double t) const;    
  /// For nested Visitor pattern.
  virtual const std::string& name() const;
};

/// Term structure using linear interpolation of zero coupon bond prices.
class TSLinear : public TSBootstrap {
public:      
  /// General constructor.
  inline TSLinear(const Array<double,1>& xT,               ///< Time line of maturities .
		          const Array<double,1>& xB                ///< Time T(0) forward zero coupon bond prices for these maturities - thus the first bond price is always 1 
                  ) : TSBootstrap(xT.copy(),xB.copy()) { };
  /// Conversion constructor.
  TSLinear(const Array<double,1>& xT,const TermStructure& xts);                
  /// Virtual copy constructor
  virtual TSLinear* pointer_to_copy() const;
  /** Function returning the interpolated value, which for term structures of
      interest rates is the maturity t, time T(0) forward zero coupon bond price. */
  virtual double operator()(double t) const;    
  /// For nested Visitor pattern.
  virtual const std::string& name() const;
};

}

#endif


