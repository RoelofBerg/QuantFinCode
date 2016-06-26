/** \file  InterfaceCubature.hpp
    \brief Header file declaring interface routines between Blitz++ classes and Johnson's cubature routine.
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

#ifndef INTERFACECUBATURE_HPP
#define INTERFACECUBATURE_HPP
#include "MSWarnings.hpp"
#include <boost/function.hpp>
#include <blitz/array.h>

namespace quantfin { 
          
using blitz::Array;
using blitz::firstDim;

int cubature(boost::function<void (const Array<double,1>&,Array<double,1>&)> f,
             const Array<double,1>& xmin, 
	         const Array<double,1>& xmax, 
	         unsigned maxEval, 
		     double reqAbsError, 
			 double reqRelError, 
		     Array<double,1>& val, 
			 Array<double,1>& err);

}

#endif
