/** \file  CSV2Array.hpp
    \brief Header file defining function to read an Array<> from a comma-separated
           values (CSV) file.
           Copyright 2005, 2006, 2011 by Erik Schlögl

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

#ifndef CSV2ARRAY_HPP
#define CSV2ARRAY_HPP
#include <fstream>
#include <iostream>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <list>
#include <iterator>
#include <boost/regex.hpp>
#include <boost/function.hpp>

namespace quantfin {
          
void load_file(std::string& s, std::istream& is);
std::string make_string(const char* str);
std::string make_string_strip_spaces(const char* str);

/// Second argument is a function to convert a null-terminated C string into T
template <class T> blitz::Array<T,2> CSV2Array(std::istream& is,boost::function<T (const char *)> convert)
{
  static const boost::regex e("\n",boost::regex::normal | boost::regbase::icase);
  static const boost::regex ei(",",boost::regex::normal | boost::regbase::icase);
  static const char* null = "";
  int j,k,n;
  std::string s;
  std::list<std::string> l;
  s.erase();
  load_file(s, is);
  std::vector<T> data(0);
  size_t cols = 0;
  size_t rows = 0;
  boost::regex_split(std::back_inserter(l), s, e);
  while(l.size()) {
    std::list<std::string> li;
    s = *(l.begin());
    l.pop_front();
    boost::regex_split(std::back_inserter(li), s, ei);
    j = 0;
    while(li.size()) {
      s = *(li.begin());
      li.pop_front();
      s.append(null,1);
	  T d = convert(s.data());
      data.push_back(d);
      j++;           }
    cols = (j>cols) ? j : cols;
    rows++;  }
  blitz::Array<T,2> mat(rows,cols);
  n = 0;
  for (j=0;j<rows;j++) {
    for (k=0;k<cols;k++) {
      if (n<data.size()) {
        mat(j,k) = data[n];
        n++; }}}
  return mat;
}

}

#endif
