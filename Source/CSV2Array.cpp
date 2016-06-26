/** \file  CSV2Array.cpp
    \brief C++ source file implementing function to read an Array<> from a comma-separated
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

#include "CSV2Array.hpp"

using namespace quantfin;

std::string quantfin::make_string(const char* str)
{
  std::string result(str);
  return result;
}

/// Strip leading and trailing spaces.
std::string quantfin::make_string_strip_spaces(const char* str)
{
  int len = std::strlen(str);
  char* newstr = new char[len+1];
  int i = 0;
  while ((str[i]==' ')&&(str[i]!=0)) i++;
  int j = len-1;
  while ((j>i)&&(str[j]==' ')) j--;
  int k = 0;
  for(;i<=j;i++) { 
	newstr[k] = str[i];
	k++; }
  newstr[k] = 0;
  std::string result(newstr);
  delete[] newstr;
  return result;
}
