/** \file  FiniteDifference.hpp
    \brief Header file defining classes for finite difference algorithms to solve partial differential equations.
           Copyright 2006, 2011 by Erik Schlögl

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

#ifndef FINITEDIFFERENCE_HPP
#define FINITEDIFFERENCE_HPP

#include <stdexcept>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include <boost/function.hpp>
#include "BlackScholesAsset.hpp"

namespace quantfin { 

using blitz::TinyVector;
using blitz::TinyMatrix;
using blitz::Array;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

class BoundaryCondition { 
public:
  virtual void operator()(Array<double,1>& x) = 0;
};

class LinearExtrapolationBC : public BoundaryCondition {
public:
  virtual void operator()(Array<double,1>& x);
};

class DirichletBC : public BoundaryCondition {
private:
  double lc;  ///< Constant lower boundary value.
  double uc;  ///< Constant upper boundary value.
public:
  inline DirichletBC(double xlc,double xuc) : lc(xlc),uc(xuc) { };
  virtual void operator()(Array<double,1>& x);
};

class FiniteDifference {
protected:
  Array<double,1>                           gridslice; ///< The current time slice of the finite difference grid.
  Array<double,1>                                 tmp;
  BoundaryCondition&               boundary_condition;
  static LinearExtrapolationBC linear_extrapolationBC;
  const BlackScholesAsset&                          S;
  double             r;
  double             T;
  int                N;
  int               Nj;
  double            dt;
  double            nu;
  double            dx;
  double            pu;
  double            pm;
  double            pd;
  double           edx;
public:
  FiniteDifference(const BlackScholesAsset& xS,double xr,double xT,int xN,int xNj);   
  /// Calculate the value of the underlying asset at a given grid point.
  inline double underlying(int i, ///< Time index
                           int j  ///< State index
                           ) const { return S.initial_value() * std::exp((j-Nj) * dx); };
  inline void underlying(int i,Array<double,1>& u) const;
  inline void apply_payoff(int i,boost::function<double (double)> f);    
  virtual void rollback(int from,int to);
  /// Rollback with early exercise (or similar) condition.
  virtual void rollback(int from,int to,boost::function<double (double,double)> f);
  inline double result() const { return gridslice(Nj); };
  double delta() const; 
  double gamma() const; 
};

inline void FiniteDifference::underlying(int i,Array<double,1>& u) const
{
  int j;
  u(0) = S.initial_value() * std::exp(-Nj * dx);      
  for (j=1;j<=2*Nj;j++) u(j) = u(j-1) * edx;
}

inline void FiniteDifference::apply_payoff(int i,boost::function<double (double x)> f) 
{
  int j;
  Array<double,1> u(2*Nj+1);
  underlying (i,u);
  for (j=0;j<=2*Nj;j++) gridslice(j) = f(u(j));      
}

class ImplicitFiniteDifference : public FiniteDifference {
private:
  Array<double,2> tridiag; ///< Coefficients of tridiagonal system of equations to be solved for each time step.
  Array<double,2>     rhs; ///< Right-hand side of tridiagonal system of equations to be solved for each time step.
  Array<double,2>     sol; ///< Space to hold the solution of tridiagonal system of equations to be solved for each time step.
  double lambda_lower,lambda_upper;
public:
  ImplicitFiniteDifference(const BlackScholesAsset& xS,double xr,double xT,int xN,int xNj);
  inline void set_derivative_at_boundary(double lower,double upper);
  virtual void rollback(int from,int to);
  /// Rollback with early exercise (or similar) condition.
  virtual void rollback(int from,int to,boost::function<double (double,double)> f);
};

inline void ImplicitFiniteDifference::set_derivative_at_boundary(double lower,double upper) 
{ 
  lambda_lower = lower * (underlying(0,1) - underlying(0,0)); 
  lambda_upper = upper * (underlying(0,2*Nj) - underlying(0,2*Nj-1)); 
}

class CrankNicolson : public FiniteDifference {
private:
  Array<double,2> tridiag; ///< Coefficients of tridiagonal system of equations to be solved for each time step.
  Array<double,2>     rhs; ///< Right-hand side of tridiagonal system of equations to be solved for each time step.
  Array<double,2>     sol; ///< Space to hold the solution of tridiagonal system of equations to be solved for each time step.
public:
  CrankNicolson(const BlackScholesAsset& xS,double xr,double xT,int xN,int xNj);
  virtual void rollback(int from,int to);
  /// Rollback with early exercise (or similar) condition.
  virtual void rollback(int from,int to,boost::function<double (double,double)> f);
};

}

#endif
