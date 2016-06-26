/** \file  Binomial.hpp
    \brief Header file defining classes for option pricing using binomial trees.
           Copyright 2006, 2007 by Erik Schlögl

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

#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <stdexcept>
#include <iostream>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <boost/function.hpp>
#include "BlackScholesAsset.hpp"

namespace quantfin { 

using blitz::Array;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

class BinomialLattice {
protected:
  Array<double,1>       gridslice; ///< The current time slice of the binomial lattice.
  Array<double,1>        tmpslice;
  const BlackScholesAsset&      S;
  double                        r;
  double                        T;
  int                           N;
  double                       dt;
  double                        p;
  double                        u;
  double                        d;
  double                 discount;
  double PeizerPratt(double z) const;
public:
  BinomialLattice(const BlackScholesAsset& xS,double xr,double xT,int xN);   
  void set_CoxRossRubinstein();
  void set_JarrowRudd();
  void set_Tian();
  void set_LeisenReimer(double K);
  /// Calculate the value of the underlying asset at a given grid point.
  inline double underlying(int i, ///< Time index
                           int j  ///< State index
                           ) const;
  inline void underlying(int i,Array<double,1>& un) const;
  inline void apply_payoff(int i,boost::function<double (double)> f);    
  inline void rollback(int from,int to);
  /// Rollback with early exercise (or similar) condition.
  inline void rollback(int from,int to,boost::function<double (double,double)> f);
  inline double result() const { return gridslice(0); };
};

inline double BinomialLattice::underlying(int i,int j) const 
{ 
  double tmp = 1.0;
  i -= j;
  while (j>0) {
    tmp *= u;
    j--; }
  while (i>0) {
    tmp *= d;
    i--; }
  return S.initial_value() * tmp; 
}

inline void BinomialLattice::underlying(int i,Array<double,1>& un) const
{
  int j;
  un(0) = underlying(i,0);      
  double tmp = u/d;
  for (j=1;j<=i;j++) un(j) = un(j-1) * tmp;
}

inline void BinomialLattice::apply_payoff(int i,boost::function<double (double x)> f) 
{
  int j;
  Array<double,1> un(i+1);
  underlying (i,un);
  for (j=0;j<=i;j++) gridslice(j) = f(un(j));      
}

inline void BinomialLattice::rollback(int from,int to)
{
  int i,j;
  for (i=from-1;i>=to;i--) {
    for (j=0;j<=i;j++) tmpslice(j) = discount * (p * gridslice(j+1) + (1-p) * gridslice(j));
    gridslice = tmpslice; }
}

inline void BinomialLattice::rollback(int from,int to,boost::function<double (double,double)> f)
{
  int i,j;
  Array<double,1> un(from);
  for (i=from-1;i>=to;i--) {
    for (j=0;j<=i;j++) tmpslice(j) = discount * (p * gridslice(j+1) + (1-p) * gridslice(j));
    gridslice = tmpslice;
    // apply early exercise condition
    underlying (i,un);
    for (j=0;j<=i;j++) gridslice(j) = f(gridslice(j),un(j)); }
}

class BinLattice {
private:
  Array<double,1> data;
  int              dim;
public:
  inline BinLattice(int xdim) : data((xdim*(xdim+1))/2),dim(xdim) { };
  inline double operator()(int i,int j) const { return data((i*(i+1))/2+j); };
  inline double& operator()(int i,int j) { return data((i*(i+1))/2+j); };
  inline BinLattice& operator=(double d) { data = d; return *this; };
  inline BinLattice& operator=(const Array<double,1>& d) { data = d; return *this; };
  inline int dimension() const { return dim; };
  inline const Array<double,1>& Array() const { return data; }; 
};

std::ostream& operator<<(std::ostream& os,const BinLattice& bl);

}

#endif
