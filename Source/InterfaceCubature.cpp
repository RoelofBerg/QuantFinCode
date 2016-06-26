/** \file  InterfaceCubature.cpp
    \brief C++ source file implementing interface routines between Blitz++ classes and Johnson's cubature routine.
           Copyright 2011 by Erik Schlögl

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

#include "InterfaceCubature.hpp"
#include "cubature.h"

using namespace quantfin;

boost::function<void (const Array<double,1>&,Array<double,1>&)> global_function_for_cubature;

void integrand_f(unsigned ndim, const double *x, void *, unsigned fdim, double *fval)
{
  int i;
  Array<double,1> in(ndim);
  for (i=0;i<ndim;i++) in(i) = x[i];
  Array<double,1> out(fval, blitz::shape(fdim), blitz::neverDeleteData);
  global_function_for_cubature(in,out);
}

int quantfin::cubature(boost::function<void (const Array<double,1>&,Array<double,1>&)> f,
                       const Array<double,1>& xmin, 
	                   const Array<double,1>& xmax, 
	                   unsigned maxEval, 
                       double reqAbsError, 
			           double reqRelError, 
		               Array<double,1>& val, 
			           Array<double,1>& err)
{
  global_function_for_cubature = f;
  int retval = adapt_integrate(val.extent(firstDim),integrand_f,NULL,xmin.extent(firstDim),xmin.data(),xmax.data(),
		                       maxEval,reqAbsError,reqRelError,val.data(),err.data());
  return retval;
}
