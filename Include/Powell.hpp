/** \file  Powell.hpp
    \brief Header file defining template for Powell's minimisation algorithm.
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

#ifndef POWELL_HPP
#define POWELL_HPP
#include <stdexcept>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

namespace quantfin { 

using blitz::TinyVector;
using blitz::TinyMatrix;
using blitz::Array;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

namespace opt {
          
template <class F,class rettype,class linesearch,int dim>
class Powell { 
private:
  F&                           f; ///< The function. Has to provide: rettype operator()(Array<rettype,1> x)
  /** Line search object. Moves the current point (currpos) to the point that minimises the function (f)
      along a given direction. Should not change currpos if the function return value changes by less than eps.
      Has to provide: rettype operator()(F& f,TinyVector<rettype,dim>& currpos,const TinyVector<rettype,dim>& direction,rettype eps) */
  linesearch&            lsearch; 
  TinyVector<rettype,dim> minarg; ///< Argument minimising the function.
  double                     eps;
  unsigned long            maxit;
  bool                 converged;
  class done { };
public:
  inline Powell(F& xf,rettype xeps,linesearch xlsearch,unsigned long xmaxit = 10000000) 
    : f(xf),eps(xeps),lsearch(xlsearch),maxit(xmaxit),converged(false) { };
  rettype solve(TinyVector<rettype,dim>& currpos);
};

template <class F,class rettype,class linesearch,int dim>
rettype Powell<F,rettype,linesearch,dim>::solve(TinyVector<rettype,dim>& currpos)
{
  unsigned long n,m,k;
  rettype                   result = 0;
  TinyVector<rettype,dim> startpos = currpos;
  int                            N = dim;
  Array<TinyVector<rettype,dim>,1> directions(dim);
  converged = false;
  try {
    for (m=0;m<std::max((long unsigned int)1,maxit/(N*N));m++) {
      // initialise directions to basis vectors
      directions = 0;  
      for (k=0;k<N;k++) (directions(k))(k) = 1;
      for (k=0;k<N;k++) {
        for (n=0;n<N;n++) result = lsearch(f,currpos,directions((n+k)%N),eps,maxit);
        directions(k) = currpos - startpos;
        // check for convergence
        if (eps>sqrt(sum(sqr(directions(k))))) {
          converged = true;
          minarg    = currpos; 
          throw done(); }
        directions(k) *= 100.0;
        result = lsearch(f,currpos,directions(k),eps,maxit);
        startpos = currpos; }}}
  catch (done) { 
    return result; }
  if (!converged) throw std::runtime_error("Powell minimisation failed to converge");
}

/********************************************************************
**      Powell minimisation general version using Array<>          **
********************************************************************/
template <class F,class rettype,class linesearch>
class GeneralPowell { 
private:
  F&                           f; ///< The function. Has to provide: rettype operator()(Array<rettype,1> x)
  /** Line search object. Moves the current point (currpos) to the point that minimises the function (f)
      along a given direction. Should not change currpos if the function return value changes by less than eps.
      Has to provide: rettype operator()(F& f,Array<rettype,1>& currpos,const Array<rettype,1>& direction,rettype eps) */
  linesearch&            lsearch; 
  Array<rettype,1>        minarg; ///< Argument minimising the function.
  double                     eps;
  int                      maxit;
  bool                 converged;
  class done { };
public:
  inline GeneralPowell(F& xf,rettype xeps,linesearch& xlsearch,unsigned long xmaxit = 10000000) 
    : f(xf),eps(xeps),lsearch(xlsearch),maxit(xmaxit),converged(false) { };
  rettype solve(Array<rettype,1>& currpos);
};

template <class F,class rettype,class linesearch>
rettype GeneralPowell<F,rettype,linesearch>::solve(Array<rettype,1>& currpos)
{
  int n,m,k;
  firstIndex i;
  secondIndex j;
  rettype                   result = 0;
  Array<rettype,1>        startpos = currpos.copy();
  rettype               prevresult = f(startpos);
  int                            N = currpos.extent(firstDim);
  Array<rettype,2> directions(N,N);
  converged = false;
  try {
    for (m=0;m<std::max(1,maxit/(N*N));m++) {
      // initialise directions to basis vectors
      directions = (i==j);  
      for (k=0;k<N;k++) {
        for (n=0;n<N;n++) result = lsearch(f,currpos,directions((n+k)%N,Range::all()),eps,std::max(10,maxit/1000));
        directions(k,Range::all()) = currpos - startpos;
        // check for convergence
        if (eps>sqrt(sum(sqr(directions(k,Range::all()))))) {
          converged = true;
          minarg    = currpos; 
          throw done(); }
        result = lsearch(f,currpos,directions(k,Range::all()),eps,std::max(10,maxit/1000));
        startpos = currpos; }
      if (eps>prevresult-result) {
        converged = true;
        minarg    = currpos;
        throw done(); }
      else prevresult = result; }}
  catch (done) { 
    return result; }
  if (!converged) {
	std::cerr << "Powell minimisation failed to converge\nCurrent position: " << currpos << "\nResult: " << result << std::endl;
	return result; // remove this!!!
	throw std::runtime_error("Powell minimisation failed to converge"); }
}

}}

#endif
