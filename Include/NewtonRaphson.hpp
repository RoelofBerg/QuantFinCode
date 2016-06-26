/** \file  NewtonRaphson.hpp
    \brief Header file defining template for multidimensional Newton-Raphson algorithm.
           Copyright 2005 by Erik Schlögl

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

#ifndef NEWTONRAPHSON_HPP
#define NEWTONRAPHSON_HPP
#include <stdexcept>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include "InterfaceCLAPACK.hpp"

namespace quantfin { 

using blitz::TinyVector;
using blitz::TinyMatrix;
using blitz::Array;
using blitz::Range;
using blitz::firstDim;
using blitz::firstIndex;
using blitz::secondIndex;

namespace opt {
          
template <class F> 
class NumericalJacobian {
private:
  /** The function. Has to provide: Array<double,1> operator()(Array<double,1> x)
                                    int argdim()
                                    int retdim()  */
  F&                           f; 
  double                       d;  ///< The perturbation of x to calculate the Jacobian by finite differences.
public:
  inline NumericalJacobian(F& xf,double delta = 1E-6) : f(xf),d(delta) { };
  inline int argdim() { return f.argdim(); };
  inline int retdim() { return f.retdim(); };
  inline Array<double,1> operator()(Array<double,1>& x) { return f(x); };
  Array<double,2> Jacobian(Array<double,1>& x);    
};
          
template <class F>
class NewtonRaphson { 
private:
  /** The function. Has to provide: Array<double,1> operator()(Array<double,1> x)
                                and Array<double,2> Jacobian(Array<double,1> x)
                                    int argdim()
                                    int retdim()  */
  F&                           f; 
  Array<double,1>     x_solution; ///< Argument yielding the target function value.
  Array<double,1>         target; ///< Target function value.
  double                     eps;
  unsigned long            maxit;
  bool                 converged;
  class done { };
public:
  inline NewtonRaphson(F& xf,double xeps,unsigned long xmaxit = 100000) 
    : f(xf),x_solution(xf.argdim()),target(xf.retdim()),eps(xeps),maxit(xmaxit),converged(false) { };
  Array<double,1> solve(Array<double,1>& currpos,const Array<double,1>& xtarget);
};

template <class F>
Array<double,1> NewtonRaphson<F>::solve(Array<double,1>& currpos,const Array<double,1>& xtarget)
{
  unsigned long n;
  converged = false;
  Array<double,1> x = currpos.copy();
  Array<double,2> dx(x.extent(firstDim),1);
  target = xtarget;
  Array<double,1> y(target.extent(firstDim)); 
  Array<double,2> y_old(target.extent(firstDim),1);
  y_old(Range::all(),0) = target - f(x); 
  for (n=0;n<maxit;n++) {
    Array<double,2> J = f.Jacobian(x);
    interfaceCLAPACK::SolveLinear(J,dx,y_old);
    x += dx(Range::all(),0);
    y  = target - f(x);
    if ((max(abs(dx))<eps)||(max(abs(y-y_old(Range::all(),0)))<eps)) {
      converged  = true;
      x_solution = x;
      return x; }
    y_old(Range::all(),0) = y; }
  if (!converged) throw std::runtime_error("Multidimensional Newton-Raphson failed to converge");
}

template <class F>
Array<double,2> NumericalJacobian<F>::Jacobian(Array<double,1>& x)
{
  int i,j;
  Array<double,2> result(f.retdim(),f.argdim()); 
  Array<double,1> fx(f(x).copy());
  Array<double,1> dx(f.argdim());
  dx = 0.0;
  for (j=0;j<f.argdim();j++) {
    dx(j) = d;
    Array<double,1> xdx(x + dx);
    Array<double,1> fxdx(f(xdx).copy());  
    for (i=0;i<f.retdim();i++) result(i,j) = (fxdx(i)-fx(i))/d;
    dx(j) = 0.0; }
  return result;
}

}}

#endif
