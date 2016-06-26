/** \file  Polynomial.hpp
    \brief Header file defining routines for polynomials.
           Copyright 2005, 2008, 2015 by Erik Schlögl

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

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP
#include <stdexcept>
#include <complex>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "QFExceptions.hpp"

namespace quantfin { 

using blitz::TinyVector;
using blitz::TinyMatrix;
using blitz::Array;
using blitz::firstDim;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;
using std::complex;

class Polynomial {
private:
  Array<complex<double>,1> c; ///< Polynomial coefficients in ascending order (i.e. coefficient of highest power last)
  Array<complex<double>,1> r; ///< Roots of the polynomial.
  bool roots_available;
  const static int MR = 8;
  const static int MT = 10;
  static double EPS;
  static double EPSS;
  int degree;
  void calc_roots();         
public:
  /// Default constructor.
  inline Polynomial() : c(1),r(0),degree(0) { roots_available = false; c(0) = 1.0; };
  /// Constructor.
  inline Polynomial(const Array<complex<double>,1>& coefficients) 
    : c(coefficients.copy()),r(coefficients.extent(firstDim)-1),degree(coefficients.extent(firstDim)-1) { roots_available = false; };
  Polynomial(const Array<double,1>& coefficients); 
  /// Copy constructor.
  Polynomial(const Polynomial& p);
  /// Function evaluation.
  complex<double> operator()(const complex<double>& x) const;
  /// Root finding.
  const Array<complex<double>,1>& roots();
  /// Assignment operations.
  Polynomial& operator=(const Polynomial& p);
  Polynomial& operator+=(const Polynomial& p);
  Polynomial& operator-=(const Polynomial& p);
  Polynomial& operator*=(const Polynomial& p);
  Polynomial& operator+=(double d);
  Polynomial& operator-=(double d);
  Polynomial& operator*=(double d);
  /// Accessors.
  int Degree() const;
  inline const Array<complex<double>,1>& coefficients() const { return c; };
};

/// Arithmetic operations.
Polynomial operator*(double c,const Polynomial& p);

}

#endif
