/** \file  StringForm.hpp
    \brief Header file defining classes for manipulating the StringForm of output.

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

#ifndef STRINGFORM_HPP
#define STRINGFORM_HPP

#include <iostream>
#include <sstream>
#include <cmath>

class Bound_StringForm;     // StringForm plus value

class StringForm {
public:
  enum str_scheme { LEFT, RIGHT, CENTER };
private:
  friend class Bound_StringForm;
  int prc;      ///< precision
  int wdt;      ///< width, 0 means as wide as necessary
  std::ios_base::fmtflags fmt;      ///< general, scientific, or fixed
  bool pl;      ///< explicit plus
  bool tz;      ///< trailing zeros
  char fl;      ///< fill character
  str_scheme str_sch;
public:
  explicit StringForm(int p = 6) : prc(p) { fmt = std::ios_base::fmtflags(0); wdt = 0; pl = false; tz = false; fl = ' '; };
  /// Make a Bound_StringForm for (*this) and d
  Bound_StringForm operator()(double d) const;
  Bound_StringForm operator()(const char* c) const;
  StringForm& scientific() { fmt = std::ios_base::scientific; return *this; };
  StringForm& fixed() { fmt = std::ios_base::fixed; return *this; };
  StringForm& general() { fmt = std::ios_base::fmtflags(0); return *this; };
  StringForm& stringscheme(str_scheme s) { str_sch = s; return *this; };
  StringForm& precision(int p) { prc = p; return *this; };
  StringForm& width(int w) { wdt = w; return *this; };
  StringForm& fill(char c) { fl = c; return *this; };
};

class Bound_StringForm {
public:
  std::ostringstream s;
  Bound_StringForm(const StringForm& ff,double v);
  Bound_StringForm(const StringForm& ff,const char* v);
};

inline Bound_StringForm StringForm::operator()(double d) const { return Bound_StringForm(*this,d); }

inline Bound_StringForm StringForm::operator()(const char* c) const { return Bound_StringForm(*this,c); }

std::ostream& operator<<(std::ostream& os,const Bound_StringForm& bf);


#endif
