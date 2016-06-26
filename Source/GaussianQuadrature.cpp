/** \file  GaussianQuadrature.cpp
    \brief C++ source file implementing classes for Gaussian quadrature.
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

#include "GaussianQuadrature.hpp"
#include "Polynomial.hpp"

using namespace quantfin;

double GaussianQuadrature::integrate(boost::function<double (double t)> f) const
{
  int i;
  double result = 0.0;
  for (i=1;i<=n;i++) result += w(i)*f(x(i))/W(x(i));
  return result;     
}

/* Given n this routine initialises the Arrays of abscissas x (sorted largest first) and weights w for an n point Gauss-Hermite quadrature. */
void GaussHermite::gausshermite()
{
  throw std::logic_error("Implementation cannot be distributed due to copyright issues - an alternative, freely distributable implementation is planned for a future release");
}

double GaussHermite::W(double t) const
{
  return std::exp(-t*t);       
}

double GaussLaguerre::W(double t) const
{
  return std::exp(-t);       
}

/// Initialise the weights w and abscissas x.
void GaussLaguerre::initialise()
{
  int i;
  // Create Laguerre polynomial L_n of order n+1 using recurrence relation - L_{n+1} is required to calculate weights below
  Array<complex<double>,1> one(2);
  one = 1.0, -1.0;
  Polynomial currL(one); // 1-x
  Polynomial prevL;      // 1
  // 2n+1-x:
  Polynomial two_n_plus_one_minus_x(one);
  for (i=2;i<=n+1;i++) {
    two_n_plus_one_minus_x += 2.0;
    Polynomial tmp(currL);
    tmp   *= two_n_plus_one_minus_x;
	prevL *= -(i-1.0);
	tmp   += prevL;
	tmp   *= 1.0/double(i);
	prevL  = currL;
	currL  = tmp; }
  // Abscissas x are the roots of the Laguerre polynomial
  const Array<complex<double>,1>& roots = prevL.roots();
  for (i=1;i<=n;i++) x(i) = roots(i-1).real(); // note that indexing starts from 1, not 0 in this particular implementation
  // Calculate weights w - using eq. 14 at http://mathworld.wolfram.com/Laguerre-GaussQuadrature.html
  for (i=1;i<=n;i++) {
    double tmp = (n+1)*currL(x(i)).real();
	w(i) = x(i) / (tmp*tmp); }
}

double GaussLobattoKronrod::W(double t) const
{
  return std::exp(-t);       
}

/// Initialise the weights w and abscissas x.
void GaussLobattoKronrod::initialise()
{
  throw std::logic_error("GaussLobattoKronrod not yet implemented");
}
