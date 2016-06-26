/** \file  QFArrayUtil.cpp
    \brief C++ source file implementing utility functions for Blitz++ Arrays.
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

#include "QFArrayUtil.hpp"

using namespace quantfin;

double quantfin::const_func(const Array<double,1>& x)
{
  return 1.0;
}

/** Assuming that T is a vector of doubles sorted in ascending order, find the index i such that
    T(i) < xi <= T(i+1). */
int quantfin::find_segment(double xi,const Array<double,1>& T) 
{
  if (xi<T(0)-std::max(1e-7,std::abs(T(0)*1e-7))) throw std::out_of_range("cannot extrapolate");
  for (int i=1;i<T.extent(firstDim);i++) if (xi<=T(i)+std::max(1e-7,std::abs(T(i)*1e-7))) return i-1;
  throw std::out_of_range("cannot extrapolate");
}

void quantfin::CSVwrite(std::ostream& os,const Array<double,1>& v)
{
  int i;
  for (i=0;i<v.extent(firstDim)-1;i++) os << v(i) << ',';
  os << v(v.extent(firstDim)-1) << std::endl;
}

void quantfin::CSVwrite(std::ostream& os,const Array<double,2>& A)
{
  int i,j;
  for (i=0;i<A.extent(firstDim);i++) {
    for (j=0;j<A.extent(secondDim)-1;j++) os << A(i,j) << ',';
	os << A(i,A.extent(secondDim)-1) << std::endl; }
}

