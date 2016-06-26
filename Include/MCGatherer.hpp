/** \file  MCGatherer.hpp
    \brief Header file defining classes to collect Monte Carlo simulations.
           Copyright 2006, 2008, 2010, 2013 by Erik Schlögl

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

#ifndef MCGATHERER_HPP
#define MCGATHERER_HPP
#include <cmath>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

namespace quantfin { 

using blitz::Array;
using blitz::firstIndex;
using blitz::secondIndex;
using blitz::firstDim;
using blitz::secondDim;

template <class T>
class MCGatherer {
private:
  boost::shared_ptr<Array<T,1> > histogram_buckets;
  boost::shared_ptr<Array<size_t,1> > histogram_count;
  T    sum;
  T   sum2;
  size_t c;    
  void add2bucket(T add);
public:
  inline MCGatherer() : sum(0.0),sum2(0.0),c(0) { };
  /// Constructor for types T which require a size argument on construction.
  inline MCGatherer(size_t array_size) : sum(array_size),sum2(array_size) { reset(); };
  void set_histogram(boost::shared_ptr<Array<T,1> > xhistogram_buckets);
  inline void reset() { sum = 0.0; sum2 = 0.0; c = 0; if (histogram_buckets) (*histogram_count) = 0.0; };
  inline void operator+=(T add) { sum += add; sum2 += add*add; c++; if (histogram_buckets) add2bucket(add); };
  inline void operator*=(T f) { sum *= f; sum2 *= f*f; };
  inline T mean() const { return T(sum/double(c)); };
  inline T stddev() const { return T(sqrt((sum2/double(c)-mean()*mean())/double(c-1))); };
  inline T variance() const { return T((sum2/double(c)-mean()*mean())/double(c-1)); };
  inline double max_stddev() const; // Must be specialised!
  inline size_t number_of_simulations() const { return c; };
  boost::shared_ptr<Array<double,2> > histogram() const;
};

template <class T>
void MCGatherer<T>::set_histogram(boost::shared_ptr<Array<T,1> > xhistogram_buckets)
{
  histogram_buckets.reset(xhistogram_buckets);
  histogram_count.reset(new Array<T,1>(xhistogram_buckets->extent(firstDim)+1));
  (*histogram_count) = 0.0;
}

template <class T>
void MCGatherer<T>::add2bucket(T add)
{
  int i = 0;
  while ((i<histogram_buckets->extent(firstDim))&&((*histogram_buckets)(i)>add)) i++;
  (*histogram_count)(i)++;
}

template <class T>
boost::shared_ptr<Array<double,2> > MCGatherer<T>::histogram() const
{
  int i;
  if (!histogram_buckets) throw std::logic_error("Histogram not generated");
  boost::shared_ptr<Array<double,2> > result(new Array<double,2>(2,histogram_count->extent(firstDim)));
  double prev = std::numeric_limits<T>::min();
  for (i=0;i<histogram_buckets->extent(firstDim);i++) {
    (*result)(0,i) = (*histogram_buckets)(i);
	(*result)(1,i) = double((*histogram_count)(i))/(double(c)*((*histogram_buckets)(i)-prev)); 
	prev           = (*histogram_buckets)(i); }
  (*result)(0,i) = std::numeric_limits<T>::max();
  (*result)(1,i) = double((*histogram_count)(i))/(double(c)*(std::numeric_limits<T>::max()-prev)); 
  return result;
}

/// Specialisation of function returning the maximum standard deviation of a set of MC estimates; specialised for double.
template <>
inline double MCGatherer<double>::max_stddev() const
{
  return stddev();
}

/// Specialisation of MCGatherer for array of values, allowing some values to be used as control variates.
template <>
class MCGatherer<Array<double,1> > {
private:
  Array<double,1>                   sum;
  Array<double,1>                  sum2;
  Array<double,1>             cv_values;
  size_t                              c;    
  Array<double,2>         covar,weights;
  Array<int,1>               cv_indices;
  firstIndex                        idx;
  secondIndex                       jdx;
  bool                             CVon;
  bool                    weights_fixed;
  int                number_of_controls;
  int                      target_index;
public:
  inline MCGatherer(int array_size) 
	: sum(array_size),sum2(array_size),covar(array_size,array_size),weights(array_size,1),cv_values(array_size),cv_indices(array_size),
	  CVon(false),weights_fixed(false) { reset(); };
  inline void reset() { sum = 0.0; sum2 = 0.0; c = 0; covar = 0.0; };
  inline void set_control_variate(bool cv) { CVon = cv; if (cv) weights_fixed = false; };
  inline void operator+=(const Array<double,1>& add);
  inline void operator*=(const Array<double,1>& f) { sum *= f; sum2 *= f*f; if (CVon) covar = covar(idx,jdx) * f(idx)*f(jdx); };
  inline Array<double,1> mean() const { return Array<double,1>(sum/double(c)); };
  inline double mean(int i) const { return (mean())(i); };
  inline Array<double,1> stddev() const { return Array<double,1>(sqrt((sum2/double(c)-mean()*mean())/double(c-1))); };
  inline double stddev(int i) const { return (stddev())(i); };
  inline Array<double,1> variance() const { return Array<double,1>((sum2/double(c)-mean()*mean())/double(c-1)); };
  inline double variance(int i) const { return (variance())(i); };
  inline size_t number_of_simulations() const { return c; };
  /// Control variate estimate when control variate weights have been fixed.
  double CVestimate(int target,int CV,double CV_expectation) const;
  /// Control variate estimate standard deviation when control variate weights have been fixed.
  double CVestimate_stddev(int target,int CV) const;
  double CVestimate(int target,const Array<int,1>& CV,const Array<double,1>& CV_expectation) const;
  double CVweight(int target,int CV) const;
  Array<double,2> CVweight(int target,const Array<int,1>& CV) const;
  inline double max_stddev() const { blitz::Array<double,1> tmp(stddev()); return blitz::max(tmp); };
  void fix_weights(int target,const Array<int,1>& CV,const Array<double,1>& CV_expectation);
  inline int dimension() const { return sum.extent(blitz::firstDim); };
};

inline void MCGatherer<Array<double,1> >::operator+=(const Array<double,1>& add) 
{ 
  int i;
  if (weights_fixed) { // Control variate estimate when control variate weights have been fixed.
    double cv = 0.0;
	for (i=0;i<number_of_controls;i++) cv += weights(i) * (add(cv_indices(i)) - cv_values(i)); 
	for (i=0;i<add.extent(blitz::firstDim);i++) {
	  double tmp;
	  if (i==target_index) tmp = add(i) - cv;
	  else                 tmp = add(i);
	  sum(i)  += tmp;
	  sum2(i) += tmp*tmp; }}
  else {
    sum  += add; 
	sum2 += add*add; }
  c++; 
  if (CVon) covar = covar(idx,jdx) + add(idx)*add(jdx); 
}

}

#endif
