/** \file  MCEngine.cpp
    \brief C++ source file implementing classes for Monte Carlo simulation.
           Copyright 2006, 2007, 2008, 2009, 2010, 2013 by Erik Schl�gl

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

#include "MCEngine.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

Array<double,1> MCPayoff::payoffArray(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values)
{
  Array<double,1> result(1);
  result(0) = (*this)(underlying_values,numeraire_values);
  return result;
}

MCEuropeanCall::MCEuropeanCall(double tzero,double tend,int asset_index,double strike)
  : MCPayoff(1,1),K(strike)
{
  timeline = tzero, tend;
  index    = asset_index, 1;     
}

double MCEuropeanCall::operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values)
{
  return numeraire_values(0)/numeraire_values(1) * std::max(0.0,underlying_values(0)-K);
}

double MCPayoffList::operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values) 
{
  int i;
  std::list<boost::shared_ptr<MCPayoff> >::const_iterator payoff_iterator;
  std::list<boost::shared_ptr<Array<int,1> > >::const_iterator time_mapping_iterator;
  std::list<boost::shared_ptr<Array<int,1> > >::const_iterator index_mapping_iterator;
  std::list<double>::const_iterator coeff_iterator;
  double result = 0.0;
  for ((payoff_iterator=payoffs.begin()),(time_mapping_iterator=time_mappings.begin()),(index_mapping_iterator=index_mappings.begin()),coeff_iterator=coefficients.begin();
  	   payoff_iterator!=payoffs.end();
  	   (++payoff_iterator),(++time_mapping_iterator),(++index_mapping_iterator),++coeff_iterator) {
    // initialise arguments and pass to function
	int n = (*payoff_iterator)->index.extent(secondDim);
	underlying_tmp.resize(n);
	numeraire_tmp.resize((*payoff_iterator)->timeline.extent(firstDim));
    for (i=0;i<numeraire_tmp.extent(firstDim);i++) numeraire_tmp(i) = numeraire_values((**time_mapping_iterator)(i));
    for (i=0;i<underlying_tmp.extent(firstDim);i++) 
	  underlying_tmp(i) = underlying_values((**index_mapping_iterator)(i));
	result += (*coeff_iterator) * (**payoff_iterator)(underlying_tmp,numeraire_tmp); }
  return result;
}

Array<double,1> MCPayoffList::payoffArray(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values) 
{
  int i,j;
  Array<double,1> result(payoffs.size());
  std::list<boost::shared_ptr<MCPayoff> >::const_iterator payoff_iterator;
  std::list<boost::shared_ptr<Array<int,1> > >::const_iterator time_mapping_iterator;
  std::list<boost::shared_ptr<Array<int,1> > >::const_iterator index_mapping_iterator;
  std::list<double>::const_iterator coeff_iterator;
  for ((j = 0),(payoff_iterator=payoffs.begin()),(time_mapping_iterator=time_mappings.begin()),(index_mapping_iterator=index_mappings.begin()),coeff_iterator=coefficients.begin();
  	   payoff_iterator!=payoffs.end();
  	   (++payoff_iterator),(++time_mapping_iterator),(++index_mapping_iterator),(++coeff_iterator),j++) {
    // initialise arguments and pass to function
	int n = (*payoff_iterator)->index.extent(secondDim);
	underlying_tmp.resize(n);
	numeraire_tmp.resize((*payoff_iterator)->timeline.extent(firstDim));
    for (i=0;i<numeraire_tmp.extent(firstDim);i++) numeraire_tmp(i) = numeraire_values((**time_mapping_iterator)(i));
    for (i=0;i<underlying_tmp.extent(firstDim);i++) 
	  underlying_tmp(i) = underlying_values((**index_mapping_iterator)(i));
	result(j) = (*coeff_iterator) * (**payoff_iterator)(underlying_tmp,numeraire_tmp); }
  return result;
}

void MCPayoffList::push_back(boost::shared_ptr<MCPayoff> xpayoff,double xcoeff)
{
  int i,j;
  std::list<boost::shared_ptr<Array<int,1> > >::iterator mapping_iterator;
  firstIndex idx;
  if (payoffs.size()==0) {
    timeline.resize(xpayoff->timeline.extent(firstDim));
	timeline = xpayoff->timeline;
	index.resize(2,xpayoff->index.extent(secondDim));
	index = xpayoff->index; 
    boost::shared_ptr<Array<int,1> > newtimemapping(new Array<int,1>(xpayoff->timeline.extent(firstDim)));
	*newtimemapping = idx; 
	time_mappings.push_back(newtimemapping); 
    boost::shared_ptr<Array<int,1> > newindexmapping(new Array<int,1>(xpayoff->index.extent(secondDim)));
	*newindexmapping = idx; 
	index_mappings.push_back(newindexmapping); }
  else {
    // merge time lines and map
    Array<double,1> merged_timeline(unique_merge(timeline,xpayoff->timeline));
	if (merged_timeline.extent(firstDim)>timeline.extent(firstDim)) {
	  // adjust existing mappings
	  Array<int,1> timeline_map(timeline.extent(firstDim));
	  j = 0;
	  for (i=0;i<timeline.extent(firstDim);i++) {
	    while (timeline(i)>merged_timeline(j)) j++;
		timeline_map(i) = j; }
      for (mapping_iterator=time_mappings.begin();mapping_iterator!=time_mappings.end();++mapping_iterator) {
		for (i=0;i<(**mapping_iterator).extent(firstDim);i++) {
		  (**mapping_iterator)(i) = timeline_map((**mapping_iterator)(i)); }}
	  for (i=0;i<index.extent(secondDim);i++) index(1,i) = timeline_map(index(1,i)); }
    boost::shared_ptr<Array<int,1> > newtimemapping(new Array<int,1>(xpayoff->timeline.extent(firstDim)));
    *newtimemapping = idx; 
	if (merged_timeline.extent(firstDim)>xpayoff->timeline.extent(firstDim)) {
	  j = 0;
	  for (i=0;i<xpayoff->timeline.extent(firstDim);i++) {
	    while (xpayoff->timeline(i)>merged_timeline(j)) j++;
		(*newtimemapping)(i) = j; }}
	time_mappings.push_back(newtimemapping); 
    timeline.resize(merged_timeline.extent(firstDim));
	timeline = merged_timeline;
    // merge indices and map
	boost::shared_ptr<Array<int,1> > newindexmapping(new Array<int,1>(xpayoff->index.extent(secondDim)));
	int oldcount = index.extent(secondDim);
	int newcount = oldcount;
	for (i=0;i<newindexmapping->extent(firstDim);i++) {
	  j = 0;
	  while (((xpayoff->index(0,i)!=index(0,j))||
		      ((*newtimemapping)(xpayoff->index(1,i))!=index(1,j)))
             &&(j<index.extent(secondDim))) j++;
	  if (j<index.extent(secondDim)) (*newindexmapping)(i) = j;
	  else {
        (*newindexmapping)(i) =  newcount; 
		newcount++; }}
	index_mappings.push_back(newindexmapping); 
    index.resizeAndPreserve(2,newcount);  
	newcount = oldcount;
	for (i=0;i<newindexmapping->extent(firstDim);i++) {
	  j = 0;
	  while (((xpayoff->index(0,i)!=index(0,j))||
		      ((*newtimemapping)(xpayoff->index(1,i))!=index(1,j)))
             &&(j<oldcount)) j++;
	  if (j>=oldcount) {
	    index(0,newcount) = xpayoff->index(0,i);
	    index(1,newcount) = (*newtimemapping)(xpayoff->index(1,i));
		newcount++; }}}
  payoffs.push_back(xpayoff);  
  coefficients.push_back(xcoeff);
}

void MCPayoffList::debug_print()
{
  using std::cout;
  using std::endl;
  cout << "Time line: " << timeline << endl;
  cout << "Index: " << index << endl;
  std::list<boost::shared_ptr<Array<int,1> > >::iterator mapping_iterator;
  cout << "Time mappings: " << endl;
  for (mapping_iterator=time_mappings.begin();mapping_iterator!=time_mappings.end();++mapping_iterator) {
	cout << **mapping_iterator << endl; }
  cout << "Index mappings: " << endl;
  for (mapping_iterator=index_mappings.begin();mapping_iterator!=index_mappings.end();++mapping_iterator) {
	cout << **mapping_iterator << endl; }
}

MCDiscreteArithmeticMeanFixedStrike::MCDiscreteArithmeticMeanFixedStrike(int asset_index,const Array<double,1>& T,double xK,
																		 int number_of_observations_in_existing_average,double existing_average)
  : MCPayoff(T,T.extent(firstDim)-1),K(xK),number_of_observations_in_existing_average_(number_of_observations_in_existing_average),existing_average_(existing_average)
{
  blitz::firstIndex  idx;
  blitz::secondIndex jdx;
  index = idx * (jdx+1);
  index(0,blitz::Range::all()) = asset_index;
}

double MCDiscreteArithmeticMeanFixedStrike::operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values)
{
  int i;
  double avg = number_of_observations_in_existing_average_ * existing_average_;
  for (i=0;i<index.extent(secondDim);i++) avg += underlying_values(i);
  avg /= i + number_of_observations_in_existing_average_;
  return numeraire_values(0)/numeraire_values(i-1) * std::max(0.0,avg-K);
}

