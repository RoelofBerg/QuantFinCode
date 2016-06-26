/** \file  LS_test.cpp
    \brief Test and demonstration program for Longstaff/Schwartz (1998).
           Copyright 2006, 2010 by Erik Schlögl

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

#include <iostream>
#include <cstdlib> 
#include <boost/bind.hpp>
#include "LongstaffSchwartz.hpp"
#include "MCGatherer.hpp"
#include "Payoff.hpp"

using namespace quantfin;

int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  try {
    Array<double,2> paths(8,4);
    paths = 1.00, 1.09, 1.08, 1.34,
            1.00, 1.16, 1.26, 1.54,
            1.00, 1.22, 1.07, 1.03,
            1.00, 0.93, 0.97, 0.92,
            1.00, 1.11, 1.56, 1.52,
            1.00, 0.76, 0.77, 0.90,
            1.00, 0.92, 0.84, 1.01,
            1.00, 0.88, 1.22, 1.34;
    Array<double,1> T(4);
    T = 0.0, 1.0, 2.0, 3.0;
    double K = 1.10;
    Payoff put(K,-1);
    boost::function<double (double)> f;
    f = boost::bind(std::mem_fun(&Payoff::operator()),&put,_1);
    FlatTermStructure ts(0.06,0.0,10.0);
	// Simplest version: one state variable, polynomial basis function
    LongstaffSchwartzExerciseBoundary1D ls(T,ts,paths,f,2);
    /* Now do the valuation as in the example in the paper - note that properly the valuation paths should
       be independent of the paths used to estimate the exercise boundary. */
    MCGatherer<double> MCestimate;
    Array<double,1> path(4);
    for (i=0;i<8;i++) {
      path = paths(i,Range::all());
	  std::cout << "Apply: " << ls.apply(path) << ' ';
      MCestimate += ls.apply(path); }
    std::cout << MCestimate.mean() << std::endl;
	// Test using general version (this version would allow for multiple state variables and arbitrary basis functions)
    Array<double,3> genpaths(8,4,1);
	genpaths(Range::all(),Range::all(),0) = paths; 
	Array<double,2> numeraire_values(8,4);
	for (i=0;i<4;i++) numeraire_values(Range::all(),i) = 1.0/ts(T(i));
	boost::function<double (double,const Array<double,1>&)> payoff = boost::bind(LSArrayAdapter,_1,_2,f,0);
	std::vector<boost::function<double (double,const Array<double,1>&)> > basisfunctions;
	int degree = 2;
	Array<int,1> p(1);
	for (i=0;i<=degree;i++) {
	  p(0) = i;
	  add_polynomial_basis_function(basisfunctions,p); }
	LongstaffSchwartzExerciseBoundary genls(T,genpaths,numeraire_values,payoff,basisfunctions);
	MCestimate.reset();
    Array<double,2> genpath(4,1);
	Array<double,1> num(4);
    for (i=0;i<8;i++) {
      genpath = genpaths(i,Range::all(),Range::all());
	  num     = numeraire_values(i,Range::all());
	  std::cout << "Apply: " << genls.apply(genpath,num) << ' ';
      MCestimate += genls.apply(genpath,num); }
    std::cout << MCestimate.mean() << std::endl;
	MCestimate.reset();
    for (i=0;i<8;i++) {
      genpath = genpaths(i,Range::all(),Range::all());
	  num     = numeraire_values(i,Range::all()) * ts(T(3));
	  std::cout << "Apply: " << genls.apply(genpath,num) << ' ';
      MCestimate += genls.apply(genpath,num); }
    std::cout << MCestimate.mean() << std::endl;
	// Most general version: would allow payoffs to depend on state variable history (e.g. barriers, lookbacks, average options)
	boost::function<double (const Array<double,1>&,const Array<double,2>&)> rebpayoff = boost::bind(REBAdapter,_1,_2,f,0);
	std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> > rebbasis_functions;
	for (i=0;i<=degree;i++) {
	  p(0) = i;
	  add_polynomial_basis_function(rebbasis_functions,p); }
	RegressionExerciseBoundary reb(T,genpaths,numeraire_values,rebpayoff,rebbasis_functions);
	MCestimate.reset();
    for (i=0;i<8;i++) {
      genpath = genpaths(i,Range::all(),Range::all());
	  num     = numeraire_values(i,Range::all());
	  std::cout << "REB apply: " << reb.apply(genpath,num) << ' ';
      MCestimate += reb.apply(genpath,num); }
    std::cout << MCestimate.mean() << std::endl;
	MCestimate.reset();
    for (i=0;i<8;i++) {
      genpath = genpaths(i,Range::all(),Range::all());
	  num     = numeraire_values(i,Range::all()) * ts(T(3));
	  std::cout << "REB apply: " << reb.apply(genpath,num) << ' ';
      MCestimate += reb.apply(genpath,num); }
    std::cout << MCestimate.mean() << std::endl;

	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
