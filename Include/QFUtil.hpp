/** \file  QFUtil.hpp
    \brief Header file defining utility functions.
           Copyright 2005, 2011 by Erik Schlögl

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

#ifndef QFUTIL_HPP
#define QFUTIL_HPP
#include <string>
#include <iostream>
#include <vector>

namespace quantfin {

bool verify_compare(double x,double y,double eps);

void load_file(std::string& s, std::istream& is);

template <class T>
bool between(const T& x,const T& a,const T& b)
{
  T lower = (a<b) ? a : b;
  T upper = (a<b) ? b : a;
  return ((x>=lower)&&(x<=upper));    
}

template <class T>
T sign(const T& x)
{
  return (x<0) ? -1 : ((x==0) ? 0 : 1);            
}

template <class T>
int get_first_index(std::vector<T> vec,T what)
{
  int i;
  for (i=0;i<vec.size();i++) if (what==vec[i]) break;
  if (i==vec.size()) i = -1; // not found
  return i;
}

}

#endif
