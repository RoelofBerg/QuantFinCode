/** \file  MCGatherer.cpp
    \brief Header file implementing classes to collect Monte Carlo simulations.
           Copyright 2010 by Erik Schlögl

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

#include "MCGatherer.hpp"
#include "InterfaceCLAPACK.hpp"

using namespace quantfin;

double MCGatherer<Array<double,1> >::CVestimate(int target,int CV,double CV_expectation) const
{
  double b = CVweight(target,CV);
  return mean(target) - b*(mean(CV) - CV_expectation);
}

double MCGatherer<Array<double,1> >::CVestimate_stddev(int target,int CV) const
{
  double rhosq = 0.0;
  if (CVon) {
    double cov = (covar(target,CV)/number_of_simulations() - mean(target)*mean(CV))/number_of_simulations(); 
	rhosq = cov * cov / (variance(CV)*variance(target)); }
  return std::sqrt(variance(target)*(1.0-rhosq));
}

double MCGatherer<Array<double,1> >::CVestimate(int target,const Array<int,1>& CV,const Array<double,1>& CV_expectation) const
{
  int i;
  int n = CV.extent(blitz::firstDim);
  Array<double,2> b(CVweight(target,CV));
  double result = mean(target);
  for (i=0;i<n;i++) result -= b(i,0) * (mean(CV(i)) - CV_expectation(i));
  return result;
}

double MCGatherer<Array<double,1> >::CVweight(int target,int CV) const
{
  double b = 0.0;
  if (CVon) {
	b = (covar(target,CV)/number_of_simulations() - mean(target)*mean(CV)) / (variance(CV)*number_of_simulations()); }
  return b;
}

Array<double,2> MCGatherer<Array<double,1> >::CVweight(int target,const Array<int,1>& CV) const
{
  int i,j;
  int n = CV.extent(blitz::firstDim);
  Array<double,2> covCV(n,n),cov(n,1),b(n,1);
  b = 0.0;
  if (CVon) {
    for (i=0;i<n;i++) {
      double meani = mean(CV(i));
	  for (j=0;j<i;j++) {
	    covCV(i,j) = covCV(j,i) = covar(CV(i),CV(j))/number_of_simulations() - meani*mean(CV(j)); }
	  covCV(i,i) = covar(CV(i),CV(i))/number_of_simulations() - meani*meani;
	  cov(i,0)   = covar(target,CV(i))/number_of_simulations() - mean(CV(i))*mean(target); }
	interfaceCLAPACK::SolveLinear(covCV,b,cov); }
  return b; 
}

void MCGatherer<Array<double,1> >::fix_weights(int target,const Array<int,1>& CV,const Array<double,1>& CV_expectation) 
{ 
  int i;
  number_of_controls = CV.extent(blitz::firstDim);
  Array<double,2> w(number_of_controls,1);
  w = CVweight(target,CV);
  for (i=0;i<number_of_controls;i++) {
    weights(i) = w(i,0);
    cv_indices(i) = CV(i);
	cv_values(i) = CV_expectation(i); }
  target_index = target;
  weights_fixed = true; 
  CVon = false; 
  reset();
}

