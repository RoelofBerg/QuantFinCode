/** \file  QFQuasiRandom.hpp
    \brief Header file declaring classes for quasi-random number generators.
           Copyright 2010, 2013 by Erik Schlögl (except as noted below)

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

#ifndef QFQUASIRANDOM_HPP
#define QFQUASIRANDOM_HPP
#include <boost/shared_ptr.hpp>
#include <boost/math/distributions/normal.hpp>
#include <blitz/array.h>

namespace quantfin {

using blitz::Array;

extern const char default_direction_number_file[];

/** Class for Sobol sequence generator based on code (sobol.cc) by Frances Y. Kuo and Stephen Joe.
    -----------------------------------------------------------------------------
    Licence pertaining to sobol.cc and the accompanying sets of direction numbers
    -----------------------------------------------------------------------------
    Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
    
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
    
        * Neither the names of the copyright holders nor the names of the
          University of New South Wales and the University of Waikato
          and its contributors may be used to endorse or promote products derived
          from this software without specific prior written permission.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
class Sobol {
private:
  Array<double,1>                current_point;
  Array<unsigned,1>                   currentX;
  boost::shared_ptr<Array<unsigned,2> >      V;
  size_t                                  maxN;
  double                               scaling;
  unsigned                            currentC;
  unsigned                            currenti;
  int                                      dim;
  void initialise(const char* direction_number_file = default_direction_number_file);
public:
  Sobol(int dimension,size_t maximum_number_of_points,const char* direction_number_file = default_direction_number_file);
  Array<double,1>& random();
  Array<double,1>& random(const Array<double,1>& shift);
  void reset();
  inline size_t number_of_points() const { return maxN; };
  inline void set_number_of_points(size_t maximum_number_of_points,const char* direction_number_file = default_direction_number_file) 
    { maxN = maximum_number_of_points; initialise(direction_number_file); };
};

class SobolArray {
private:
  Array<double,2>    contents;
  Sobol                   rng;
public:
  inline SobolArray(int rows,int cols,size_t maximum_number_of_points,const char* direction_number_file = default_direction_number_file) 
	  : contents(rows,cols),rng(rows*cols,maximum_number_of_points,direction_number_file) { };
  /// Returns an array of draws from Sobol sequence.
  inline Array<double,2>& random();
};

inline Array<double,2>& SobolArray::random()
{
  int i,j,k;
  Array<double,1> rnd(rng.random());
  k = 0;
  for (i=0;i<contents.extent(blitz::firstDim);i++) {
    for (j=0;j<contents.extent(blitz::secondDim);j++) {
	  contents(i,j) = rnd(k);
	  k++; }}
  return contents;
}

class SobolArrayNormal {
private:
  Array<double,2>    contents;
  Sobol                   rng;
  boost::math::normal  normal;
  int                       n;
  const double        epsilon;
public:
  inline SobolArrayNormal(int rows,int cols,size_t maximum_number_of_points,const char* direction_number_file = default_direction_number_file) 
	  : contents(rows,cols),rng(rows*cols,maximum_number_of_points,direction_number_file),n(rows*cols),epsilon(1e-64) {  };
  /// Returns an array of draws from Sobol sequence.
  inline Array<double,2>& random();
  inline Array<double,2>& random(Array<double,2>& shift);
  inline void reset() { rng.reset(); };
  inline size_t number_of_points() const { return rng.number_of_points(); };
  inline void set_number_of_points(size_t maximum_number_of_points,const char* direction_number_file = default_direction_number_file) 
    { rng.set_number_of_points(maximum_number_of_points,direction_number_file); };
};

inline Array<double,2>& SobolArrayNormal::random()
{
  int i,j,k;
  Array<double,1> rnd(rng.random());
  k = 0;
  for (i=0;i<contents.extent(blitz::firstDim);i++) {
    for (j=0;j<contents.extent(blitz::secondDim);j++) {
	  double x = rnd(k);
	  if (x==0) x = epsilon;
	  contents(i,j) = boost::math::quantile(normal,x);
	  k++; }}
  return contents;
}

inline Array<double,2>& SobolArrayNormal::random(Array<double,2>& shift)
{
  int i,j,k;
  Array<double,1> tmp(shift.dataFirst(),n,blitz::neverDeleteData);
  Array<double,1> rnd(rng.random(tmp));
  k = 0;
  for (i=0;i<contents.extent(blitz::firstDim);i++) {
    for (j=0;j<contents.extent(blitz::secondDim);j++) {
	  double x = rnd(k);
	  if (x==0) x = epsilon;
	  contents(i,j) = boost::math::quantile(normal,x);
	  k++; }}
  return contents;
}

}

#endif
