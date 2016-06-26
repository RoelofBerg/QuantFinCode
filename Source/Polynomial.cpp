/** \file  Polynomial.cpp
    \brief C++ source file implementing routines for polynomials.
           Copyright 2005, 2015 by Erik Schlögl

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

#include <iostream>
#include <cmath>
#include "Polynomial.hpp"
#include "rootfinder.h"

using namespace quantfin;

double Polynomial::EPS = 2.0e-12;
double Polynomial::EPSS = 1.0e-14;

Polynomial::Polynomial(const Polynomial& p)
  : c(p.c.copy()),r(p.r.copy()),degree(p.Degree()),roots_available(p.roots_available)
{ }

Polynomial::Polynomial(const Array<double,1>& coefficients)
  : c(coefficients.extent(firstDim)),r(coefficients.extent(firstDim)-1),degree(coefficients.extent(firstDim)-1)
{
  int idx;
  for (idx=0;idx<=degree;idx++) c(idx) = coefficients(idx);
}

complex<double> Polynomial::operator()(const complex<double>& x) const
{
  int i;
  complex<double> y;
  y = c(Degree());
  for (i=Degree()-1;i>=0;i--) {
    y *= x;
    y += c(i); }
  return y; 
}

void Polynomial::calc_roots()
{
  int j,jj;
  bool poly = (degree>0);
  if (poly) {
    bool allzero = true;
	for (j=1;j<c.extent(firstDim);j++) {
	  if (c(j)!=0.0) allzero = false; }
	poly = !allzero; }
  if (poly) {
	bool polish_roots_after = true;
	bool use_roots_as_starting_points = false;
	cmplx_roots_gen(r,c,degree,polish_roots_after,use_roots_as_starting_points); 
	roots_available = true; }
  else throw DegeneratePolynomial();
}

int Polynomial::Degree() const 
{ 
  int j;
  bool allzero = true;
  for (j=1;j<c.extent(firstDim);j++) {
    if (c(j)!=0.0) allzero = false; }
  return (allzero) ? 0 : degree; 
}

Polynomial& Polynomial::operator=(const Polynomial& p)
{
  int i;
  if (degree!=p.Degree()) {
    c.resize(p.Degree()+1);   
    degree = p.Degree();     
    r.resize(p.Degree()); }
  for (i=0;i<=p.Degree();i++) c(i) = p.c(i);
  roots_available = p.roots_available;
  if (roots_available) for (i=0;i<p.Degree();i++) r(i) = p.r(i);
  return *this;
}

Polynomial& Polynomial::operator+=(const Polynomial& p)
{
  int i;
  if (degree<p.Degree()) {
    c.resizeAndPreserve(p.Degree()+1);   
    for (i=degree+1;i<=p.Degree();i++) c(i) = 0.0;   
    degree = p.Degree();     
    r.resize(p.Degree()); }
  for (i=0;i<=p.Degree();i++) c(i) += p.c(i);
  roots_available = false;
  return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& p)
{
  int i;
  if (degree<p.Degree()) {
    c.resizeAndPreserve(p.Degree()+1);   
    for (i=degree+1;i<=p.Degree();i++) c(i) = 0.0;   
    degree = p.Degree();     
    r.resize(p.Degree()); }
  for (i=0;i<=p.Degree();i++) c(i) -= p.c(i);
  roots_available = false;
  return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& p)
{
  int i,j;
  int newdegree = degree + p.Degree();
  Array<complex<double>,1> nc(newdegree+1);   
  nc = 0.0;   
  r.resizeAndPreserve(newdegree); 
  for (i=0;i<=p.Degree();i++) {
    for (j=0;j<=degree;j++) nc(i+j) += p.c(i)*c(j); }
  c.resize(nc.shape());
  c = nc;
  if (roots_available && p.roots_available) {
    for (i=0;i<p.Degree();i++) r(i+degree) = p.r(i); }
  else roots_available = false;
  degree = newdegree;
  return *this;
}

Polynomial& Polynomial::operator+=(double d)
{
  c(0) += d;
  roots_available = false;
  return *this;
}

Polynomial& Polynomial::operator-=(double d)
{
  c(0) -= d;
  roots_available = false;
  return *this;
}

Polynomial& Polynomial::operator*=(double d)
{
  int j;
  for (j=0;j<=Degree();j++) c(j) *= d; 
  return *this;
}

Polynomial operator*(double c,const Polynomial& p)
{
  Polynomial y(p);
  y *= c;
  return y;           
}

const Array<complex<double>,1>& Polynomial::roots() 
{ 
  if (Degree()==0) throw std::logic_error("No roots in polynomial of degree zero");
  if (!roots_available) calc_roots(); 
  return r; 
}
