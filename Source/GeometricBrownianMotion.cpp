/** \file  GeometricBrownianMotion.cpp
    \brief C++ source file implementing class representing a multidimensional 
           (e.g. multi-asset) geometric Brownian motion.
           Copyright 2006, 2013 by Erik Schlögl

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

#include "GeometricBrownianMotion.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

GeometricBrownianMotion::GeometricBrownianMotion(std::vector<const BlackScholesAsset*>& xunderlying)
  : underlying(xunderlying),T(NULL),timeline_(NULL),asset_values(xunderlying.size(),2),time_mapping(NULL),
    dW(xunderlying[0]->volatility_function().factors()),vol_lvl(xunderlying[0]->volatility_function().factors())
{ }

GeometricBrownianMotion::~GeometricBrownianMotion()
{
  if (T) delete T;                                                  
  if (timeline_) delete timeline_;                                                  
  if (time_mapping) delete time_mapping;                                                  
}

/// Query the dimension of the process.
int GeometricBrownianMotion::dimension() const
{
  return underlying.size();
}

/// Generate a realisation of the process under the martingale measure associated with deterministic bond prices.
void GeometricBrownianMotion::operator()(Array<double,2>& underlying_values,const Array<double,2>& x,const TermStructure& ts)
{
  int i,j;
  if (!T) throw std::logic_error("Timeline not set in GeometricBrownianMotion::operator()");
  for (j=1;j<asset_values.extent(secondDim);j++) { // Loop over time line
    double fwd = ts((*T)(j-1))/ts((*T)(j));
	dW = x(Range::all(),j-1) * std::sqrt((*T)(j)-(*T)(j-1));
    for (i=0;i<asset_values.extent(firstDim);i++) { // Loop over assets
      asset_values(i,j) = asset_values(i,j-1) * fwd * underlying[i]->dividend_discount((*T)(j-1),(*T)(j)) * underlying[i]->DoleansExp((*T)(j-1),(*T)(j),dW); }}
  // set return values
  for (i=0;i<underlying_values.extent(firstDim);i++) {
	for (j=0;j<underlying_values.extent(secondDim);j++) underlying_values(i,j) = asset_values(i,(*time_mapping)(j)); }
}

/// Generate a realisation of the process under the martingale measure associated with a given numeraire asset.
void GeometricBrownianMotion::operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const Array<double,2>& x,const TermStructure& ts,int numeraire_index)
{
  int i,j;
  if (!T) throw std::logic_error("Timeline not set in GeometricBrownianMotion::operator()");
  if (numeraire_index>=0) {
    // set return values
    for (j=1;j<asset_values.extent(secondDim);j++) { // Loop over time line
      double fwd = ts((*T)(j-1))/ts((*T)(j));
	  dW = x(Range::all(),j-1) * std::sqrt((*T)(j)-(*T)(j-1));
	  add_drift((*T)(j-1),(*T)(j),numeraire_index);
      for (i=0;i<asset_values.extent(firstDim);i++) { // Loop over assets
        asset_values(i,j) = asset_values(i,j-1) * fwd * underlying[i]->dividend_discount((*T)(j-1),(*T)(j)) * underlying[i]->DoleansExp((*T)(j-1),(*T)(j),dW); }}
    // set return values
    for (i=0;i<underlying_values.extent(firstDim);i++) {
	  for (j=0;j<underlying_values.extent(secondDim);j++) underlying_values(i,j) = asset_values(i,(*time_mapping)(j)); }
    // set numeraire values
    // Note that numeraire values must be adjusted for reinvestment of dividends.
    for (j=0;j<underlying_values.extent(secondDim);j++) 
	  numeraire_values(j) = underlying_values(numeraire_index,j) / underlying[numeraire_index]->dividend_discount((*timeline_)(0),(*timeline_)(j)); }
  else (*this)(underlying_values,x,ts);
}

void GeometricBrownianMotion::add_drift(double tstart,double tend,int numeraire_index)
{
  if (numeraire_index>=0) {
    if (!((underlying[numeraire_index]->volatility_function()).get_volatility_level(tstart,tend,vol_lvl)))
      throw std::logic_error("Unable to get volatility levels in GeometricBrownianMotion::add_drift()");
    vol_lvl *= tend - tstart;
    dW += vol_lvl; }
}

/// Set process timeline.
bool GeometricBrownianMotion::set_timeline(const Array<double,1>& timeline)
{
  int i,j;
  /* check timeline conformity with volatility of underlying assets
     all volatilities must be constant on each time segment */
  Array<double,1> segments(timeline.copy());
  for (i=0;i<underlying.size();i++) {
    Array<double,1> tmp_segments = unique_merge(segments,(underlying[i]->volatility_function()).segments(timeline(0),timeline(timeline.extent(firstDim)-1)-timeline(0)));
    segments.resize(tmp_segments.extent(firstDim));
    segments = tmp_segments; }
  if (T) delete T;
  T = new Array<double,1>(segments.copy());
  if (timeline_) delete timeline_;
  timeline_ = new Array<double,1>(timeline.copy());
  // map timeline indices to segments
  if (time_mapping) delete time_mapping;
  time_mapping = new Array<int,1>(timeline.extent(firstDim));
  j = 0;
  for (i=0;i<segments.extent(firstDim);i++) {
    if (timeline(j)==segments(i)) {
      (*time_mapping)(j) = i;
      j++; }}
  // resize scratch arrays
  if (asset_values.extent(secondDim)!=T->extent(firstDim)) asset_values.resize(underlying.size(),T->extent(firstDim));
  for (i=0;i<underlying.size();i++) asset_values(i,0) = underlying[i]->initial_value();
  return T;   
}

