/** \file  TSInstruments.hpp
    \brief Header file declaring classes to represent interest rate instruments.
           Copyright 2006, 2015 by Erik Schlögl

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

#ifndef TSINSTRUMENTS_HPP
#define TSINSTRUMENTS_HPP

#include "TermStructure.hpp"

namespace quantfin {

using blitz::firstIndex;
using blitz::Range;
using blitz::toEnd;

class TSEuropeanInstrument {
private:
  double p; ///< Current price.
  double t; ///< Time at which this price is valid.
  double T; ///< Maturity.
public:
  inline TSEuropeanInstrument(double xp,double xt,double xT) : p(xp),t(xt),T(xT) { };
  virtual double payoff(const TermStructure& ts) = 0;
  inline double& price() { return p; };
  inline double price() const { return p; };
  inline double& maturity() { return T; };
  inline double maturity() const { return T; };
  inline double& today() { return t; };
  inline double today() const { return t; };
  /// Predicate for sorting instruments by maturity
  class longer_maturity {
  public:
    inline bool operator()(const TSEuropeanInstrument* first,const TSEuropeanInstrument* second) 
    { 
      return (first->T>second->T);
    }
  };
  /// Predicate for sorting instruments by maturity
  class shorter_maturity {
  public:
    inline bool operator()(const TSEuropeanInstrument* first,const TSEuropeanInstrument* second) 
    { 
      return (first->T<second->T);
    }
  };
};

class Caplet : public TSEuropeanInstrument {
private:
  double      lvl;
  double notional;
  double    delta;  ///< Length of accrual period.
public:
  inline Caplet(double xp,double xt,double xT,double xlvl,double xdelta,double xnotional = 1.0) 
    : TSEuropeanInstrument(xp,xt,xT),lvl(xlvl),notional(xnotional),delta(xdelta) { };
  virtual double payoff(const TermStructure& ts);
};

class Floorlet : public TSEuropeanInstrument {
private:
  double      lvl;
  double notional;
  double    delta;  ///< Length of accrual period.
public:
  inline Floorlet(double xp,double xt,double xT,double xlvl,double xdelta,double xnotional = 1.0) 
    : TSEuropeanInstrument(xp,xt,xT),lvl(xlvl),notional(xnotional),delta(xdelta) { };
  virtual double payoff(const TermStructure& ts);
};

class Swaption : public TSEuropeanInstrument {
private:
  double            lvl;
  double       notional;
  Array<double,1> tenor;  ///< Tenor structure of underlying swap
  int              sign;  ///< -1 = payer, 1 = receiver swaption
public:
  Swaption(double xp,double xt,double xT,double xlvl,double xdelta,int nperiods,int xsign = 1,double xnotional = 1.0); 
  virtual double payoff(const TermStructure& ts);
};

class TSBermudanInstrument {
private:
  double          p; ///< Current price.
  double          t; ///< Time at which this price is valid.
protected:
  Array<double,1> T; ///< Exercise dates.
public:
  inline TSBermudanInstrument(double xp,double xt,int n) : p(xp),t(xt),T(n) { };
  virtual double payoff(double continuation_value,const TermStructure& ts) = 0;
  inline double& price() { return p; };
  inline double price() const { return p; };
  inline const Array<double,1>& maturity() const { return T; };
  inline Array<double,1>& maturity() { return T; };
  inline double& today() { return t; };
  inline double today() const { return t; };
  int exercise_date(double now) const;
};

/// Bermudan swaption: the assumption is that it can be exercised at any reset date of the underlying swap.
class BermudanSwaption : public TSBermudanInstrument {
private:
  double            lvl;
  double       notional;
  Array<double,1> tenor;  ///< Tenor structure of underlying swap
  int              sign;  ///< -1 = payer, 1 = receiver swaption
public:
  BermudanSwaption(double xp,double xt,double xT,double xlvl,double xdelta,int nperiods,int xsign = 1,double xnotional = 1.0); 
  virtual double payoff(double continuation_value,const TermStructure& ts);
};

/// Knock-out barrier instrument
class BarrierInstrument : public TSBermudanInstrument {
private:
  TSEuropeanInstrument&            instrument;
  double                              barrier; 
  double                          barrier_ttm; ///< Time to maturity of barrier rate
  int                               direction; ///< Up: 1; Down: -1
public:
  BarrierInstrument(double xp,double xt,TSEuropeanInstrument& xinstrument,const Array<double,1>& barrier_monitoring_dates,double xbarrier,double xbarrier_ttm,int xdirection);
  virtual double payoff(double continuation_value,const TermStructure& ts);
};

}

#endif
