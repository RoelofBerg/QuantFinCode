/** \file  MCEngine.hpp
    \brief Header file declaring classes for Monte Carlo simulation.
           Copyright 2006, 2007, 2008, 2009, 2010 by Erik Schlögl

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

#ifndef MCENGINE_HPP
#define MCENGINE_HPP
#include <vector>
#include <list>
#include "MSWarnings.hpp"
#include <boost/shared_ptr.hpp>
#include "MCGatherer.hpp"
#include "TermStructure.hpp"

namespace quantfin {

using blitz::firstDim;
using blitz::secondDim;

template <class target_price_process,class controlvariate_price_process,class cv_random_variable> class MCControlVariateMapping;
     
/** Abstract base class defining the common interface for classes mapping a path asset value (or state variable) realisations
    and numeraire value realisations to a discounted payoff.
	*/
class MCPayoff {
public:
  Array<double,1> timeline;  ///< Time line collecting all event dates.
  Array<int,2>       index;  ///< A 2 x N matrix of indices, where each column represents the 
                             ///< indices of an (asset,time) combination affecting the payoff. 
  inline MCPayoff(const Array<double,1>& xtimeline,const Array<int,2>& xindex) : timeline(xtimeline),index(xindex) { };
  inline MCPayoff(const Array<double,1>& xtimeline,int number_of_indices) : timeline(xtimeline),index(2,number_of_indices) { };
  inline MCPayoff(int number_of_event_dates,int number_of_indices) : timeline(number_of_event_dates+1),index(2,number_of_indices) { };
  /// Calculate discounted payoff. 
  virtual double operator()(const Array<double,1>& underlying_values, ///< Underlying values for the (asset,time) combinations in index Array.
	                        const Array<double,1>& numeraire_values   ///< Numeraire values for the dates in timeline Array.
							) = 0;
  /// Allow for multidimensional payoff (i.e. portfolio) with meaningful default (one-dimensional) behaviour.
  virtual Array<double,1> payoffArray(const Array<double,1>& underlying_values, ///< Underlying values for the (asset,time) combinations in index Array.
	                                  const Array<double,1>& numeraire_values   ///< Numeraire values for the dates in timeline Array.
							          );
};
          
class MCPayoffList : public MCPayoff {
private:
  std::list<boost::shared_ptr<MCPayoff> > payoffs;
  std::list<double> coefficients;
  std::list<boost::shared_ptr<Array<int,1> > > time_mappings;
  std::list<boost::shared_ptr<Array<int,1> > > index_mappings;
  Array<double,1> underlying_tmp;
  Array<double,1>  numeraire_tmp;
  void debug_print();
public:
  inline MCPayoffList() : MCPayoff(1,1),underlying_tmp(16),numeraire_tmp(16) { };
  void push_back(boost::shared_ptr<MCPayoff> xpayoff,double xcoeff = 1.0);
  /// Calculate discounted payoff. 
  virtual double operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values);
  /// Multidimensional payoff (i.e. portfolio).
  virtual Array<double,1> payoffArray(const Array<double,1>& underlying_values, ///< Underlying values for the (asset,time) combinations in index Array.
	                                  const Array<double,1>& numeraire_values   ///< Numeraire values for the dates in timeline Array.
							          );
};

/** Template that maps a set of random numbers to a stochastic process to a Monte Carlo payoff.

    The template parameter class price_process must implement a member function dimension(),
	which returns the number of assets (i.e. the dimension) of the price process.
	It must also implement the member function set_timeline(const Array<double,1>& timeline),
	which sets the timeline for which asset prices will be generated.
	Furthermore, it must provide operator()(Array<double,2>& underlying_values,const random_variable& x,const TermStructure& ts) 
	and operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const random_variable& x,const TermStructure& ts,int numeraire_index),
	which generate the values of the price process at the dates in the required timeline using
	a draw x of the random variable (the type of which is given by the template parameter random_variable),
	under the martingale measure associated with, respectively, deterministic bond prices or a chosen numeraire.
    These requirements are for example implemented by the class GeometricBrownianMotion.

	@see GeometricBrownianMotion
    */
template <class price_process,class mapped_random_variable>
class MCMapping {
private:
  /** Index of underlying asset to use as the numeraire. If less than zero, use deterministic discounting
      via the initial term structure. */
  int        numeraire_index;   
  price_process&                   process;   ///< Stochastic process driving the underlying asset dynamics.
  const TermStructure&                  ts;   ///< Initial term structure.
  MCPayoff&                         payoff;   ///< Specification of the Monte Carlo payoff.
  Array<double,2>        underlying_values;
  Array<double,1> mapped_underlying_values;
  Array<double,1>         numeraire_values;
public:
  /// Constructor.
  MCMapping(MCPayoff& xpayoff,price_process& xprocess,const TermStructure& xts,int xnumeraire_index = -1); 
  /// Choose the numeraire asset for the equivalent martingale measure under which the simulation is carried out.
  bool set_numeraire_index(int xnumeraire_index);
  /// The function mapping a realisation of the (often multidimensional) random variable x to the discounted payoff.
  double mapping(mapped_random_variable x);
  /// The function mapping a realisation of the (often multidimensional) random variable x to multiple discounted payoffs.
  Array<double,1> mappingArray(mapped_random_variable x);
  template <class target_price_process,class controlvariate_price_process,class cv_random_variable> friend class MCControlVariateMapping;
  inline void print_index() {
    std::cout << payoff.index << std::endl; };
};

template <class target_price_process,class controlvariate_price_process,class cv_random_variable>
class MCControlVariateMapping {
private:
  MCMapping<target_price_process,cv_random_variable>&                         target;
  MCMapping<controlvariate_price_process,cv_random_variable>&         controlvariate;
  const Array<double,1>&                                       controlvariate_values;
  double                                                   controlvariate_values_sum;
public:
  /// Constructor.
  inline MCControlVariateMapping(MCMapping<target_price_process,cv_random_variable>& xtarget,
	                             MCMapping<controlvariate_price_process,cv_random_variable>& xcontrolvariate,
								 const Array<double,1>& xcontrolvariate_values)
	  : target(xtarget),controlvariate(xcontrolvariate),controlvariate_values(xcontrolvariate_values) { controlvariate_values_sum = sum(controlvariate_values); }; 
  /// The function mapping a realisation of the (often multidimensional) random variable x to the discounted payoff.
  double mapping(cv_random_variable x);
  /// The function mapping a realisation of the (often multidimensional) random variable x to multiple discounted payoffs.
  Array<double,1> mappingArray(cv_random_variable x);
  inline void print_index() {
    std::cout << controlvariate.payoff.index << std::endl; };
};

template <class target_price_process,class controlvariate_price_process,class cv_random_variable>
double MCControlVariateMapping<target_price_process,controlvariate_price_process,cv_random_variable>::mapping(cv_random_variable x)
{
  return target.mapping(x) - controlvariate.mapping(x) + controlvariate_values_sum;
}

template <class target_price_process,class controlvariate_price_process,class cv_random_variable>
Array<double,1> MCControlVariateMapping<target_price_process,controlvariate_price_process,cv_random_variable>::mappingArray(cv_random_variable x)
{
  return target.mappingArray(x) - controlvariate.mappingArray(x) + controlvariate_values;
}

class MCEuropeanCall : public MCPayoff {
private:
  double K; ///< Strike.
public:
  MCEuropeanCall(double tzero,double tend,int asset_index,double strike);
  virtual double operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values);
  inline double& strike() { return K; };
};

class MCDiscreteArithmeticMeanFixedStrike : public MCPayoff {
private:
  double                                        K; ///< Strike.
  int number_of_observations_in_existing_average_;
  double                        existing_average_;
public:
  MCDiscreteArithmeticMeanFixedStrike(int asset_index,const Array<double,1>& T,double xK,int number_of_observations_in_existing_average,double existing_average);
  virtual double operator()(const Array<double,1>& underlying_values,const Array<double,1>& numeraire_values);
};

template <class price_process,class mapped_random_variable>
MCMapping<price_process,mapped_random_variable>::MCMapping(MCPayoff& xpayoff,price_process& xprocess,const TermStructure& xts,int xnumeraire_index)
    : payoff(xpayoff),process(xprocess),ts(xts),
	  underlying_values(xprocess.dimension(),xpayoff.timeline.extent(firstDim)),mapped_underlying_values(xpayoff.index.extent(secondDim)),
	  numeraire_values(xpayoff.timeline.extent(firstDim)) 
{ 
  if (!set_numeraire_index(xnumeraire_index)) throw std::logic_error("Unable to set numeraire index in MCMapping"); 
  process.set_timeline(xpayoff.timeline);
}

template <class price_process,class mapped_random_variable>
bool MCMapping<price_process,mapped_random_variable>::set_numeraire_index(int xnumeraire_index)
{
  int i;
  numeraire_index = xnumeraire_index;
  if (numeraire_index == -1) { // Deterministic discounting using initial term structure.
	for (i=0;i<numeraire_values.extent(firstDim);i++) numeraire_values(i) = 1.0/ts(payoff.timeline(i)); 
    return true; }
  else return process.set_numeraire(xnumeraire_index);
}

template <class price_process,class mapped_random_variable>
double MCMapping<price_process,mapped_random_variable>::mapping(mapped_random_variable x)
{
  int i;
  /* Using the random draw x, generate the values of the price process at the dates in the required timeline, 
     under the martingale measure associated with the chosen numeraire. */
  process(underlying_values,numeraire_values,x,ts,numeraire_index);
  // Map underlying values to the payoff.
  for (i=0;i<mapped_underlying_values.extent(firstDim);i++) {
	mapped_underlying_values(i) = underlying_values(payoff.index(0,i),payoff.index(1,i)); }
  // Calculate the discounted payoff.
  return payoff(mapped_underlying_values,numeraire_values);
}

template <class price_process,class mapped_random_variable>
Array<double,1> MCMapping<price_process,mapped_random_variable>::mappingArray(mapped_random_variable x)
{
  int i;
  /* Using the random draw x, generate the values of the price process at the dates in the required timeline, 
     under the martingale measure associated with the chosen numeraire. */
  if (numeraire_index >= 0) process(underlying_values,numeraire_values,x,ts,numeraire_index);
  else                      process(underlying_values,x,ts);
  // Map underlying values to the payoff.
  for (i=0;i<mapped_underlying_values.extent(firstDim);i++) mapped_underlying_values(i) = underlying_values(payoff.index(0,i),payoff.index(1,i));
  // Calculate the discounted payoff.
  return payoff.payoffArray(mapped_underlying_values,numeraire_values);
}

}

#endif
