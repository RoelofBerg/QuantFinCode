/** \file  MCGeneric.hpp
    \brief Header file declaring a template for generic Monte Carlo simulation.
           Copyright 2007, 2010 by Erik Schlögl

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

#ifndef MCGENERIC_HPP
#define MCGENERIC_HPP
#include "MSWarnings.hpp"
#include <boost/function.hpp>
#include <boost/math/distributions/normal.hpp>
#include "MCGatherer.hpp"

namespace quantfin {

/** Template for generic Monte Carlo simulation.
    
	The template parameter class random_number_generator_type must implement a member function random(), 
	which returns a realisation of the random variable of the type given by the template parameter class argtype.
    */
template <class argtype,class rettype,class random_number_generator_type>
class MCGeneric {
private:
  boost::math::normal                                     N;
  boost::function<rettype (argtype)>                      f;  ///< Functor mapping a draw of the random variable to the Monte Carlo payoff(s).
  boost::function<argtype (argtype)>             antithetic;  ///< Functor mapping a draw of the random variable to their antithetic values.
  random_number_generator_type      random_number_generator;
public:
  inline MCGeneric(boost::function<rettype (argtype)> func,random_number_generator_type& rng) 
	  : f(func),random_number_generator(rng) { };
  inline MCGeneric(boost::function<rettype (argtype)> func,random_number_generator_type& rng,boost::function<argtype (argtype)> antifunc) 
	  : f(func),random_number_generator(rng),antithetic(antifunc) { };
  inline void set_antithetic(boost::function<rettype (argtype)> antifunc) { antithetic = antifunc; };
  void simulate(MCGatherer<rettype>& mcgatherer,unsigned long number_of_simulations);
  double simulate(MCGatherer<rettype>& mcgatherer,unsigned long initial_number_of_simulations,double required_accuracy,double confidence_level = 0.95);
};

template <class argtype,class rettype,class random_number_generator_type>
void MCGeneric<argtype,rettype,random_number_generator_type>::simulate(MCGatherer<rettype>& mcgatherer,unsigned long number_of_simulations)
{
  unsigned long i;
  if (!antithetic) {
	for (i=0;i<number_of_simulations;i++) {
	  mcgatherer += f(random_number_generator.random()); }}
  else {
	for (i=0;i<number_of_simulations;i++) {
	  argtype rnd = random_number_generator.random();
	  rettype res = f(rnd);
	  res = 0.5 * (res + f(antithetic(rnd))); 
	  mcgatherer += res; }}
}

template <class argtype,class rettype,class random_number_generator_type>
double MCGeneric<argtype,rettype,random_number_generator_type>::simulate(MCGatherer<rettype>& mcgatherer,
																 unsigned long initial_number_of_simulations,
																 double required_accuracy,
																 double confidence_level)
{
  unsigned long n = initial_number_of_simulations;
  simulate(mcgatherer,initial_number_of_simulations);
  double d = boost::math::quantile(N,confidence_level); // was N.inverse(confidence_level);
  double current_accuracy = d * mcgatherer.max_stddev();
  while (required_accuracy<current_accuracy) {
    double q = current_accuracy/required_accuracy;
	q *= q;
    unsigned long additional_simulations = n * q - n;
	if (!additional_simulations) break;
    simulate(mcgatherer,additional_simulations);
	n += additional_simulations;
	current_accuracy = d * mcgatherer.max_stddev(); }
  return current_accuracy;
}

}

#endif
