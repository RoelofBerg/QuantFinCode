/** \file  QFRandom.hpp
    \brief Header file declaring classes for random number generators.
           Copyright 2007, 2010, 2011 by Erik Schlögl

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

#ifndef QFRANDOM_HPP
#define QFRANDOM_HPP
#include <blitz/array.h>

namespace quantfin {

using blitz::firstDim;
using blitz::secondDim;

/** Wrap random_number_generator_type which implements operator() to return a realisation of the random variable of the type given by the template parameter class rntype
    to return a realisation of the random variable via the member function random()
    */
template <class random_number_generator_type,class rntype>
class RandomWrapper {
private:
  random_number_generator_type&  rng;
public:
  inline RandomWrapper(random_number_generator_type& xrng) : rng(xrng) { };
  /// Returns a draw from random_number_generator_type.
  inline rntype random() { return rng(); };
};

/** Template for an array-valued random number generator.
    
	The template parameter class random_number_generator_type must implement a member function random(), 
	which returns a realisation of the random variable of the type given by the template parameter class rntype.
    */
template <class random_number_generator_type,class rntype>
class RandomArray {
private:
  blitz::Array<rntype,2>    contents;
  random_number_generator_type&  rng;
public:
  inline RandomArray(random_number_generator_type& xrng,int rows,int cols) : contents(rows,cols),rng(xrng) { };
  /// Returns an array of draws from random_number_generator_type.
  blitz::Array<rntype,2>& random();
  inline blitz::Array<rntype,2>& previous_draw() { return contents; };
};

template <class random_number_generator_type,class rntype>
blitz::Array<rntype,2>& RandomArray<random_number_generator_type,rntype>::random()
{
  int i,j;
  for (i=0;i<contents.extent(firstDim);i++) {
	for (j=0;j<contents.extent(secondDim);j++) contents(i,j) = rng.random(); }
  return contents;
}

/// Antithetic variables mirrored around zero.
template <class T> 
T normal_antithetic(T arg)
{
  return T(-arg);
}

/// Antithetic for uniform [0,1].
template <class T> 
T uniform_antithetic(T arg)
{
  return T(1.0 - arg);
}

}

#endif
