/** \file  MCAmerican.hpp
    \brief Header file declaring classes for Monte Carlo simulation of options with early exercise features.
           Copyright 2010 by Erik Schlögl

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

#ifndef MCAMERICAN_HPP
#define MCAMERICAN_HPP
#include "MCEngine.hpp"
#include "LongstaffSchwartz.hpp"

namespace quantfin {

/** Class mapping a path asset value (or state variable) realisations and numeraire value realisations to a discounted payoff, based on an
    early exercise strategy encapsulated in the template parameter class exercise_boundary. The template parameter class exercise_boundary
	must supply the following functions:
	  -# double apply(const Array<double,2>& path,const Array<double,1>& numeraire_values) const to map a path asset value (or state variable) 
	     realisations and numeraire value realisations to a discounted payoff
	  -# const Array<double,1>& get_timeline() const returning the time points at which the early exercise condition is evaluated
	  -# int number_of_state_variables() const returning the number of assets (or state variables) required to be simulated on each path

    @see LongstaffSchwartzExerciseBoundary
	or
	@see RegressionExerciseBoundary
	*/
template <class exercise_boundary>
class LSExerciseStrategy : public MCPayoff {
private:
  exercise_boundary& boundary;
  Array<double,2>        path;
public:
  LSExerciseStrategy(exercise_boundary& xboundary);
  /// Calculate discounted payoff. 
  virtual double operator()(const Array<double,1>& underlying_values, ///< Underlying values for the (asset,time) combinations in index Array.
	                        const Array<double,1>& numeraire_values   ///< Numeraire values for the dates in timeline Array.
							);
};

/** Create a (paths,timepoints,state variables) array of state variable realisations and a (paths,timepoints) array of 
    numeraire realisations, e.g. for Longstaff/Schwartz regression.

    The template parameter class stochastic_process must implement a member function dimension(),
	which returns the number of state variables (i.e. the dimension) of the price process.
	It must also implement the member function set_timeline(const Array<double,1>& timeline),
	which sets the timeline for which state variable realisations will be generated.
	Furthermore, it must provide operator()(Array<double,2>& underlying_values,const random_variable& x,const TermStructure& ts) 
	and operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const random_variable& x,const TermStructure& ts,int numeraire_index),
	which generate the values of the price process at the dates in the required timeline using
	a draw x of the random variable (the type of which is given by the template parameter random_variable),
	under the martingale measure associated with, respectively, deterministic bond prices or a chosen numeraire.
    These requirements are for example implemented by the class GeometricBrownianMotion.

	@see GeometricBrownianMotion
	*/
template <class stochastic_process,class random_variable_generator>
class MCTrainingPaths {
private:
  Array<double,3>             paths;
  Array<double,2>   numeraire_paths;
  Array<double,2> underlying_values;
  Array<double,1>  numeraire_values;
  Array<double,1>                 T;
  stochastic_process&       process;
  random_variable_generator&    rng;
  const TermStructure&           ts;
  int               numeraire_index;
public:
  MCTrainingPaths(stochastic_process& xprocess,const Array<double,1>& xtimeline,int npaths,random_variable_generator& xrng,const TermStructure& xts,int xnumeraire_index);
  void generate_paths();
  /// Returns a paths x (time points) x assets Array of state variable paths.
  inline const Array<double,3>& state_variables() const { return paths; };
  inline const Array<double,2>& numeraires() const { return numeraire_paths; };
};

template <class stochastic_process,class random_variable_generator>
MCTrainingPaths<stochastic_process,random_variable_generator>
::MCTrainingPaths(stochastic_process& xprocess,const Array<double,1>& xtimeline,int npaths,random_variable_generator& xrng,const TermStructure& xts,int xnumeraire_index)
  : process(xprocess),paths(npaths,xtimeline.extent(firstDim),xprocess.dimension()),rng(xrng),ts(xts),numeraire_index(xnumeraire_index),
    underlying_values(xprocess.dimension(),xtimeline.extent(firstDim)),numeraire_values(xtimeline.extent(firstDim)),
	numeraire_paths(npaths,xtimeline.extent(firstDim)),T(xtimeline.copy())
{
  process.set_timeline(xtimeline);
  generate_paths();  
}

template <class stochastic_process,class random_variable_generator>
void MCTrainingPaths<stochastic_process,random_variable_generator>::generate_paths()
{
  int i,j;
  firstIndex  idx;
  secondIndex jdx;
  if (numeraire_index==-1) {
	for (j=0;j<T.extent(firstDim);j++) numeraire_values(j) = 1.0/ts(T(j)); }
  for (i=0;i<paths.extent(firstDim);i++) {
	if (numeraire_index==-1) process(underlying_values,rng.random(),ts);
	else                     process(underlying_values,numeraire_values,rng.random(),ts,numeraire_index);
	numeraire_paths(i,Range::all()) = numeraire_values;
	paths(i,Range::all(),Range::all()) = underlying_values(jdx,idx); }
}

template <class exercise_boundary>
LSExerciseStrategy<exercise_boundary>::LSExerciseStrategy(exercise_boundary& xboundary)
: MCPayoff(xboundary.get_timeline(),xboundary.get_timeline().extent(firstDim)*xboundary.number_of_state_variables()),boundary(xboundary),
  path(xboundary.get_timeline().extent(firstDim),xboundary.number_of_state_variables())
{
  // initialise index Array
  int i,j,c;
  c = 0;
  for (i=0;i<boundary.number_of_state_variables();i++) {
	for (j=0;j<timeline.extent(firstDim);j++) {
      index(0,c) = i;
	  index(1,c) = j;
	  c++; }}
}

template <class exercise_boundary>
double LSExerciseStrategy<exercise_boundary>::operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values)
{
  // reformat for LongstaffSchwartzExerciseBoundary.apply()
  int i,j,c;
  c = 0;
  for (i=0;i<boundary.number_of_state_variables();i++) {
	for (j=0;j<timeline.extent(firstDim);j++) {
	  path(j,i) = underlying_values(c);
	  c++; }}
  return boundary.apply(path,numeraire_values);
}
          
}

#endif
