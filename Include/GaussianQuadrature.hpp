/** \file  GaussianQuadrature.hpp
    \brief Header file defining classes for Gaussian quadrature.
           Copyright 2004-2006, 2008, 2009 by Erik Schlögl

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

#ifndef GAUSSIANQUADRATURE_HPP
#define GAUSSIANQUADRATURE_HPP
#include <stdexcept>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <boost/function.hpp>
#include "Constants.hpp"

namespace quantfin { 

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

class GaussianQuadrature {
protected:
  int             n;
  Array<double,1> x;                     ///< abscissas
  Array<double,1> w;                     ///< weights
  double        eps;                     ///< Relative precision.
  static const int    maxit = 10;        ///< Maximum iterations.
public:
  inline GaussianQuadrature(int xn,double xeps = 3.0e-14);
  inline const Array<double,1>& abscissas() const { return x; };
  inline const Array<double,1>& weights() const { return w; };
  double integrate(boost::function<double (double t)> f) const;
  virtual double W(double t) const = 0;
};

inline GaussianQuadrature::GaussianQuadrature(int xn,double xeps) : eps(xeps),x(xn+1),w(xn+1),n(xn) { }

class GaussHermite : public GaussianQuadrature {
private:
  void gausshermite();
public:
  inline GaussHermite(int xn);
  virtual double W(double t) const;
};

inline GaussHermite::GaussHermite(int xn) 
: GaussianQuadrature(xn)
{
  gausshermite();
}

class GaussLaguerre : public GaussianQuadrature {
private:
  void initialise();
public:
  inline GaussLaguerre(int xn);
  virtual double W(double t) const;
};

inline GaussLaguerre::GaussLaguerre(int xn) 
: GaussianQuadrature(xn)
{
  initialise();
}

class GaussLobattoKronrod : public GaussianQuadrature {
private:
  void initialise();
public:
  inline GaussLobattoKronrod(int xn);
  virtual double W(double t) const;
};

inline GaussLobattoKronrod::GaussLobattoKronrod(int xn) 
: GaussianQuadrature(xn)
{
  initialise();
}

}

#endif
