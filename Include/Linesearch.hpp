/** \file  Linesearch.hpp
    \brief Header file defining templates for line search minimisation algorithms.
           Copyright 2005, 2011 by Erik Schlögl

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

#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP
#include <stdexcept>
#include <limits>
#include <boost/function.hpp>
#include <blitz/array.h>
#include "QFUtil.hpp"
#include "QFArrayUtil.hpp"

namespace quantfin { namespace opt {
          
using blitz::Array;
using blitz::firstDim;

template <class rettype>
double ParabolicFit(double x1,double x2,double x3,rettype f1,rettype f2,rettype f3)
{
  double delta1 = x2 - x1;
  double delta3 = x2 - x3;
  double fdelta1 = f2 - f1;
  double fdelta3 = f2 - f3;
  double result = x2 - 0.5 * (delta1*delta1*fdelta3-delta3*delta3*fdelta1)/(delta1*fdelta3-delta3*fdelta1);
  return result;
}

template <class F,class rettype>
bool GeneralBracketMinimum(F& f,Array<rettype,1>& point1,Array<rettype,1>& point2,Array<rettype,1>& point3,rettype& f1,rettype& f2,rettype& f3,double lambda1,double lambda2,double lambda3,const Array<double,1>& direction,double lambda_min,double lambda_max)
{
  f1 = f(point1);   
  f2 = f(point2);   
  if (f1<f2) {
    quantfin::swap(point1,point2);
    std::swap(f1,f2); 
    std::swap(lambda1,lambda2); }
  lambda3 = std::min(lambda_max,std::max(lambda_min,lambda2+2.0*(lambda2-lambda1)));
  point3 = point2 + (lambda3-lambda2) * direction;
  f3 = f(point3);
  while ((f3<f2)&&(lambda3>lambda_min)&&(lambda3<lambda_max)) {
    point1  = point2.copy();
    f1      = f2;
    lambda1 = lambda2;
    point2  = point3.copy();
    f2      = f3;
    lambda2 = lambda3;
    lambda3 = std::min(lambda_max,std::max(lambda_min,lambda2+2.0*(lambda2-lambda1)));
    point3  = point2 + (lambda3-lambda2) * direction;
    f3      = f(point3); }
  return ((f2<=f3)&&(f2<=f1));
}

template <class F,class rettype>
class GeneralBrentLinesearch { 
private:
  boost::function<void (Array<rettype,1>&,const Array<rettype,1>&,double&,double&)> set_bounds;
  double lambda_min,lambda_max;
public:
  inline GeneralBrentLinesearch(boost::function<void (Array<rettype,1>&,const Array<rettype,1>&,double&,double&)> xset_bounds = NULL)
	  : set_bounds(xset_bounds),lambda_min(-std::numeric_limits<double>::max()),lambda_max(std::numeric_limits<double>::max()) { };
  rettype operator()(F& f,Array<rettype,1>& currpos,const Array<rettype,1>& direction,rettype eps,unsigned long maxit = 10000);
};

template <class F,class rettype>
rettype GeneralBrentLinesearch<F,rettype>::operator()(F& f,Array<rettype,1>& currpos,const Array<rettype,1>& xdirection,rettype eps,unsigned long maxit)
{
  Array<rettype,1> direction = xdirection.copy();
  direction /= sqrt(sum(sqr(xdirection))); // normalise vector
  if (set_bounds) { // set_bounds is non-NULL if search domain is bounded
	set_bounds(currpos,direction,lambda_min,lambda_max);
    if (lambda_max-lambda_min<1e-7) return f(currpos);
    if (lambda_max<1e-7) { // reverse direction
      direction *= -1.0;
      lambda_min *= -1.0; 
      lambda_max *= -1.0; }}
  unsigned long i;
  int  dim = currpos.extent(firstDim);
  int axis = 0;
  while ((direction(axis)==0)&&(axis<dim)) axis++;
  rettype fbracket1,fbracket2,fleastsofar;
  Array<rettype,1> bracket1   = currpos.copy();
  Array<rettype,1> leastsofar = currpos.copy();
  leastsofar += std::min(1.0,lambda_max/2.0) * direction;
  Array<rettype,1> bracket2   = currpos.copy();
  bool bracketed = GeneralBracketMinimum(f,bracket1,leastsofar,bracket2,fbracket1,fleastsofar,fbracket2,0.0,std::min(1.0,lambda_max/2.0),0.0,direction,lambda_min,lambda_max);
  if (!bracketed) { // return minimum at extreme end if minimum not bracketed
    if (fbracket1<fbracket2) {
      currpos = bracket1;
      return fbracket1; }
    else {
      currpos = bracket2;
      return fbracket2; }}
  Array<rettype,1> secondleastsofar      = (fbracket1>fbracket2) ? bracket2.copy() : bracket1.copy();
  rettype          fsecondleastsofar     = (fbracket1>fbracket2) ? fbracket2 : fbracket1;
  Array<rettype,1> prevsecondleastsofar  = (fbracket1>fbracket2) ? bracket1.copy() : bracket2.copy();
  rettype          fprevsecondleastsofar = (fbracket1>fbracket2) ? fbracket1 : fbracket2;
  Array<rettype,1> latesteval            = bracket2.copy();
  rettype          flatesteval           = fbracket2;
  double           stepbeforelast        = 0.0;
  double           laststep              = stepbeforelast;
  for (i=0;i<maxit;i++) {
    double delta2 = (prevsecondleastsofar(axis)-leastsofar(axis))/direction(axis);  
    double delta3 = (secondleastsofar(axis)-leastsofar(axis))/direction(axis);  
    double lambda = ParabolicFit(0.0,delta2,delta3,fleastsofar,fprevsecondleastsofar,fsecondleastsofar);
    /* To be acceptable, the parabolic step must fall within the bounding interval 
       and imply less than half the movement of the step before last.
       Otherwise, bisect the larger segment */
    delta2 = (bracket1(axis)-leastsofar(axis))/direction(axis);  
    delta3 = (bracket2(axis)-leastsofar(axis))/direction(axis);  
    if (!(between(lambda,delta2,delta3)&&(std::abs(lambda)<0.5*stepbeforelast)&&(std::abs(lambda)>eps))) lambda = 0.5 * ((std::abs(delta2)>std::abs(delta3)) ? delta2 : delta3);
    // Return if done.
    if (std::abs(lambda)<eps) {
      currpos = latesteval;
      return flatesteval; }
    latesteval     = leastsofar + lambda*direction; 
    flatesteval    = f(latesteval);
    // Housekeeping
    stepbeforelast = laststep;
    laststep       = std::abs(lambda);
    if (flatesteval<fleastsofar) {
      if (lambda*delta2>0.0) {
        bracket2  = leastsofar;
        fbracket2 = fleastsofar; }
      else {
        bracket1  = leastsofar;
        fbracket1 = fleastsofar; }}
    else {
      if (lambda*delta2>0.0) {
        bracket1  = latesteval;
        fbracket1 = flatesteval; }
      else {
        bracket2  = latesteval;
        fbracket2 = flatesteval; }}
    if (flatesteval<fprevsecondleastsofar) {
      if (flatesteval<fsecondleastsofar) {
        prevsecondleastsofar  = secondleastsofar;
        fprevsecondleastsofar = fsecondleastsofar;
        if (flatesteval<fleastsofar) {
          secondleastsofar  = leastsofar;
          fsecondleastsofar = fleastsofar;
          leastsofar        = latesteval;
          fleastsofar       = flatesteval; }
        else {
          secondleastsofar  = latesteval;
          fsecondleastsofar = flatesteval; }}
      else {
        prevsecondleastsofar  = latesteval;
        fprevsecondleastsofar = flatesteval; }}}
  if ((fsecondleastsofar-fleastsofar)/(std::abs(fleastsofar)+1e-6)<1e-6) {
	currpos = leastsofar;
	return fleastsofar; }
  else {
	currpos = leastsofar;
	return fleastsofar; }
  throw std::runtime_error("Brent minimisation failed to converge");
}

}}

#endif
