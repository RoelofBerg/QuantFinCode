/** \file  TermStructure.hpp
    \brief Header file defining classes for the term structure of interest rates.
           Copyright 2003, 2005, 2007, 2010, 2011, 2013, 2014 by Erik Schlögl

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

     Term structures are represented by classes derived from the abstract base TermStructure
     */

#ifndef TERMSTRUCTURE_HPP
#define TERMSTRUCTURE_HPP
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include "QFVisitor.hpp"

namespace quantfin {
          
using blitz::Array;
using blitz::firstDim;

class DeterministicCashflow;

class TermStructure : public QFNestedVisitable {
private:
  class fit_func {
  private:
    TermStructure& ts;
    const std::vector<DeterministicCashflow>& cashflows;
  public:
    inline fit_func(const std::vector<DeterministicCashflow>& xcashflows,TermStructure& xts) : cashflows(xcashflows),ts(xts) { };
    double operator()(Array<double,1>& zcb); 
  };
protected:
  Array<double,1>           T;   ///< Time line.
  /** Interest rate term structures are represented as zero coupon bond prices.
      Note that the time line starts at T(0), not zero, thus if T(0) is greater
      than "today", B contains actually T(0)-forward bond prices.
      Note also that the first bond price is always 1.
    */
  Array<double,1>           B;   
public:
  /// General constructor.
  inline TermStructure(const Array<double,1>& xT,         ///< Time line of maturities .
		               const Array<double,1>& xB          ///< Time T(0) forward zero coupon bond prices for these maturities - thus the first bond price is always 1 
                       ) : T(xT.copy()),B(xB.copy()) { };
  /// Constructor for term structure with n (uninitialised) data points.
  inline TermStructure(int n) : T(n),B(n) { };
  /// Virtual copy constructor
  virtual TermStructure* pointer_to_copy() const = 0;
  /// Reinitialise (if necessary) when zero coupon bond prices have changed
  virtual void reinitialise();
  /// Fit the term structure (as well as possible given the time line) to a set of market values.
  virtual void approximate_fit(std::vector<DeterministicCashflow> cashflows,double eps = 1E-09);
  /// For nested Visitor pattern.
  virtual void accept(QFNestedVisitor& visitor) const;
  inline double operator[](int idx) const;    ///< Returns B[idx].
  inline double& operator[](int idx);    ///< Returns B[idx].
  /// Query time horizon
  inline double time_horizon() const { return T(T.extent(firstDim)-1); };
  /// Number of dates on the time line.
  inline int length() const;
  /// i-th date on the timeline.
  inline double t(int idx) const; 
  /// read access to time line
  inline const Array<double,1>& timeline() const { return T; };
  /** Function returning the interpolated value, which for term structures of
      interest rates is the maturity t, time T(0) forward zero coupon bond price. */
  virtual double operator()(double t) const = 0;    
  /// Find the segment of timeline containing t.
  int find_segment(double t) const;
  /// Continuously compounded forward yield for the accrual period from u to v.
  inline double forward_yield(double u,double v) const { return (std::log((*this)(u)/(*this)(v)))/(v-u); };
  /// Continuously compounded forward yields for the accrual periods defined by the Array of dates xT.
  Array<double,1> ccfwd(const Array<double,1>& xT) const;
  /// Get all delta-compounded forward rates for the tenor structure xT.
  Array<double,1> simple_rate(const Array<double,1>& xT,double delta) const; 
  /// Get the delta-compounded forward rate, observed today, for the forward date xT. 
  double simple_rate(double xT,double delta) const;
  /** Get all delta-compounded forward swap rates for the tenor structure xT, keeping the swap end date fixed at the last element of xT.
      This is the set of forward swap rates as in section 8 of 
      Jamshidian, F. (1997) LIBOR and Swap Market Models and Measures, <I>Finance and Stochastics</I>
      1(4), pp. 293-330. */
  Array<double,1> swaps(const Array<double,1>& xT,double delta) const; 
  /// Get forward swap rate, valid today.
  double swap(double T0,    ///< Start date of the forward swap.
              int n,        ///< Length of the forward swap in number of delta periods.
              double delta  ///< Length of time between two payments of the swap.
              ) const;
  /// Get forward swap rate, valid today, for a given tenor (no assumption on length of tenor periods).
  double swap(const Array<double,1>& tenor    ///< Swap tenor structure.
              ) const;
  /** Present value of a (swap) basis point.
      This is defined here as \f[ \sum_{i=1}^n \frac{B(t,T_0+i\delta)}{B(t,T_0)} \f] where
      \f$ B(t,T) \f$ denotes the price of a zero coupon bond at time t (today), maturing at time T.
    */
  double pvbp(double T0,int n,double delta) const;
  double pvbp(const Array<double,1>& tenor) const;
  /// Find the maturity such that the zero coupon bond price for this maturity matches the target value
  double find_maturity(double target) const;
};

inline double TermStructure::operator[](int idx) const { return B(idx); }

inline double& TermStructure::operator[](int idx) { return B(idx); }

inline int TermStructure::length() const { return T.extent(firstDim); }

inline double TermStructure::t(int idx) const { return T(idx); }

/// Flat (constant yield) term structure.
class FlatTermStructure : public TermStructure {
private: 
  double level;
public:      
  /// Flat term structure constructor.
  FlatTermStructure(double lvl,    ///< Flat level of continuously compounded interest rates.
                    double tzero,  ///< Start of maturity time line.
                    double horizon ///< End of maturity time line.
                    );
  /// Virtual copy constructor
  virtual FlatTermStructure* pointer_to_copy() const;
  /// Reinitialise (if necessary) when zero coupon bond prices have changed
  virtual void reinitialise();
  /** Function returning the interpolated value, which for term structures of
      interest rates is the maturity t, time T(0) forward zero coupon bond price. */
  virtual double operator()(double t) const;    
  /// For nested Visitor pattern.
  virtual void accept(QFNestedVisitor& visitor) const;
  virtual const std::string& name() const;
};

}

#endif


