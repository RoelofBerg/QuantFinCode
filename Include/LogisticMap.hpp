/** \file  LogisticMap.hpp
    \brief Header file defining class for logistic mapping.
           Copyright 2006 by Erik Schlögl

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

#ifndef LOGISTICMAP_HPP
#define LOGISTICMAP_HPP

#include <algorithm>

class LogisticMap {
private:
  double a,b;
public:
  inline LogisticMap(double xa = 0.0,double xb = 1.0) : a(xa),b(xb) { };
  inline double min() const { return a; }; 
  inline double max() const { return b; }; 
  inline void set_min(double xa) { a = xa; }; 
  inline void set_max(double xb) { b = xb; }; 
  inline void operator()(double& x) const;
  inline void inverse(double& x) const;
  inline void operator*=(double x);
};

inline void LogisticMap::operator()(double& x) const
{
  x = a + (b-a)/(1.0+std::exp(-x));
}

inline void LogisticMap::inverse(double& x) const
{
  x = -std::log((b-a)/(x-a)-1.0);
}

inline void LogisticMap::operator*=(double x)
{
  a *= x;
  b *= x;
  if (x<0.0) std::swap(a,b);      
}

#endif
