/** \file  StringForm.cpp
    \brief C++ source file implementing classes for manipulating the StringForm of output.
 		   Redistribution and use in source and binary forms, with or without modification, are permitted provided
		   that the following conditions are met:
		
           -# Redistributions of source code must retain this list of conditions and the following disclaimer.

           -# Redistributions in binary form must reproduce this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

           THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
		   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
		   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
		   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
		   OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
		   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
		   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
		   POSSIBILITY OF SUCH DAMAGE.
           */

#include "StringForm.hpp"

Bound_StringForm::Bound_StringForm(const StringForm& ff,double v) 
{ 
  int wdt = ff.wdt;
  if (v>=0.0) {
    s << ' ';
    wdt--; }
  s.precision(ff.prc);
  s.setf(ff.fmt,std::ios_base::floatfield);
  s.setf(std::ios_base::left);
  s.width(wdt);
  s.fill(ff.fl);
  s << v;
}

Bound_StringForm::Bound_StringForm(const StringForm& ff,const char* c) 
{ 
  int i;
  if (ff.str_sch == StringForm::CENTER) {
    const char* p = c;
    int len = 0;
    while (*p++) len++;
    int rest = ff.wdt - len;
    int  add = rest % 2;
    rest /= 2;
    for (i=0;i<rest;i++) s << ff.fl;
    s << c;
    for (i=0;i<rest+add;i++) s << ff.fl; }
  else {
    s.width(ff.wdt);
    s.fill(ff.fl);
    s.setf((ff.str_sch==StringForm::RIGHT) ? std::ios_base::right : std::ios_base::left,std::ios_base::adjustfield);
    s << c; }
}

std::ostream& operator<<(std::ostream& os,const Bound_StringForm& bf)
{
  return os << bf.s.str().c_str();
}
