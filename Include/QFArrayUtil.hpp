/** \file  QFArrayUtil.hpp
    \brief Header file defining utility functions for Blitz++ Arrays.
           Copyright 2005, 2006, 2007, 2015 by Erik Schlögl

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

#ifndef QFARRAYUTIL_HPP
#define QFARRAYUTIL_HPP
#include <stdexcept>
#include <iostream>
#include <blitz/array.h>

namespace quantfin {

using blitz::Array;
using blitz::firstIndex;
using blitz::firstDim;
using blitz::secondDim;

double const_func(const Array<double,1>& x);

template <class T,int N>
void swap(Array<T,N>& A,Array<T,N>& B)
{
  Array<T,N> tmp = A.copy();
  A = B.copy();
  B = tmp;
}

template <int N>
Array<double,N> logistic_transform(const Array<double,N>& x,double a,double b)
{
  Array<double,N> result = x.copy();
  result = a + (b-a)/(1+exp(-x));
  return result;
}

template <int N>
Array<double,N> inverse_logistic_transform(const Array<double,N>& x,double a,double b)
{
  Array<double,N> result = x.copy();
  result = -log((b-a)/(x-a)-1);
  return result;
}

/** Assuming that T is a vector sorted in ascending order, find the index i such that
    T(i) < xi <= T(i+1). */
template <class X>
int find_segment(X xi,const Array<X,1>& T) 
{
  if (xi<T(0)) throw std::out_of_range("cannot extrapolate");
  for (int i=1;i<T.extent(firstDim);i++) if (xi<=T(i)) return i-1;
  throw std::out_of_range("cannot extrapolate");
}

/** Assuming that T is a vector of doubles sorted in ascending order, find the index i such that
    T(i) < xi <= T(i+1). */
int find_segment(double xi,const Array<double,1>& T); 

template <class X>
int find_first(X xi,const Array<X,1>& T) 
{
  for (int i=0;i<T.extent(firstDim);i++) if (xi==T(i)) return i;
  return -1;
}

template <class X>
void partial_sum(Array<X,1>& T)
{
  int i;
  for (i=1;i<T.extent(firstDim);i++) T(i) += T(i-1);
}

template <class X>
Array<X,1> unique_merge(const Array<X,1>& A,const Array<X,1>& B)
{
  int apos = 0;
  int bpos = 0;
  Array<X,1> result(A.extent(firstDim)+B.extent(firstDim));
  // size of the merged arrays
  int len = 0;
  while ((apos<A.extent(firstDim))||(bpos<B.extent(firstDim))) {
    if (apos<A.extent(firstDim)) {
      if (bpos<B.extent(firstDim)) {
        if (A(apos)<=B(bpos)) {
          result(len) = A(apos);
          if (B(bpos)==A(apos)) bpos++;
          apos++; }
        else {
          result(len) = B(bpos);
          bpos++; }}
      else {
        result(len) = A(apos);
        apos++; }}
    else {
      result(len) = B(bpos);
      bpos++; }
    len++; }
  result.resizeAndPreserve(len);
  return result;
}

/// Determine whether A is a subset of B
template <class X>
bool subset(const Array<X,1>& A,const Array<X,1>& B)
{
  int i;
  for (i=0;i<A.extent(firstDim);i++) {
	if (-1==find_first(A(i),B)) return false; }
  return true;
}

void CSVwrite(std::ostream& os,const Array<double,1>& v);
void CSVwrite(std::ostream& os,const Array<double,2>& A);

}

#endif
