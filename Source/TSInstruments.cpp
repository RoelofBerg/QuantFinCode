/** \file  TSInstruments.cpp
    \brief C++ source file implementing classes to represent interest rate instruments.
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

#include "TSInstruments.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

double Caplet::payoff(const TermStructure& ts)
{
  double bond = ts(maturity()+delta);
  return notional * std::max(0.0,1.0-bond*(1.0+delta*lvl));   
}

double Floorlet::payoff(const TermStructure& ts)
{
  double bond = ts(maturity()+delta);
  return notional * std::max(0.0,-(1.0-bond*(1.0+delta*lvl)));   
}

Swaption::Swaption(double xp,double xt,double xT,double xlvl,double xdelta,int nperiods,int xsign,double xnotional) 
    : TSEuropeanInstrument(xp,xt,xT),lvl(xlvl),notional(xnotional),tenor(nperiods+1),sign(xsign) 
{ 
  firstIndex idx;
  tenor = idx*xdelta + xT;
}

double Swaption::payoff(const TermStructure& ts)
{
  double result = std::max(0.0,sign*(lvl-ts.swap(tenor)));
  if (result>0.0) result *= ts.pvbp(tenor);
  return result;
}

int TSBermudanInstrument::exercise_date(double now) const
{
  return find_first(now,T);
}

BermudanSwaption::BermudanSwaption(double xp,double xt,double xT,double xlvl,double xdelta,int nperiods,int xsign,double xnotional) 
    : TSBermudanInstrument(xp,xt,nperiods),lvl(xlvl),notional(xnotional),tenor(nperiods+1),sign(xsign) 
{ 
  firstIndex idx;
  tenor = idx*xdelta + xT;
  T = tenor(Range(0,nperiods-1));
}

double BermudanSwaption::payoff(double continuation_value,const TermStructure& ts)
{
  double now = ts.timeline()(0);
  double result = continuation_value;
  int exercise_idx = exercise_date(now);
  if (exercise_idx>=0) {
	Array<double,1> tmp_tenor(tenor(Range(exercise_idx,toEnd)));
	double val = sign*(lvl-ts.swap(tmp_tenor)) * ts.pvbp(tmp_tenor);
	result = std::max(result,val); }
  return result;
}

BarrierInstrument::BarrierInstrument(double xp,double xt,TSEuropeanInstrument& xinstrument,const Array<double,1>& barrier_monitoring_dates,double xbarrier,double xbarrier_ttm,int xdirection)
  : TSBermudanInstrument(xp,xt,barrier_monitoring_dates.extent(firstDim)),instrument(xinstrument),barrier(xbarrier),direction(xdirection),barrier_ttm(xbarrier_ttm)
{ 
  T = barrier_monitoring_dates;
  if (instrument.maturity()!=T(T.extent(firstDim)-1)) throw std::logic_error("BarrierInstrument: last barrier monitoring date must coincide with the maturity of the payoff");
}

double BarrierInstrument::payoff(double continuation_value,const TermStructure& ts)
{
  double now = ts.timeline()(0);
  double result = continuation_value;
  if (now==instrument.maturity()) result = instrument.payoff(ts);
  int exercise_idx = exercise_date(now);
  if (exercise_idx>=0) {  // check for knock-out
	if (0.0>=direction*(barrier-ts.simple_rate(now,barrier_ttm))) result = 0.0; }
  return result;
}
