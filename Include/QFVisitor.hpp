/** \file  QFVisitor.hpp
    \brief Header file defining abstract base classes for the Visitor and nested Visitor patterns.
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

#ifndef QFVISITOR_HPP
#define QFVISITOR_HPP
#include "MSWarnings.hpp"
#include <blitz/array.h>

namespace quantfin {

using blitz::Array;
using blitz::firstDim;

class QFVisitor;
class QFNestedVisitor;

class QFVisitable {
public:
  virtual void accept(QFVisitor& visitor) const = 0;
};

class QFNestedVisitable {
public:
  virtual void accept(QFNestedVisitor& visitor) const = 0;
  virtual const std::string& name() const = 0;
};

class QFVisitor {
public:
  virtual void visit(const QFVisitable& visitable) = 0;
};

class QFNestedVisitor {
public:
  virtual void visit(const QFNestedVisitable& visitable) = 0;    
  virtual void visit(const std::string& name,const QFNestedVisitable& value) = 0; 
  // Visit functions for recognised data types.
  virtual void visit(double x) = 0;    
  virtual void visit(size_t x) = 0;    
  virtual void visit(const Array<double,1>& x) = 0;    
  virtual void visit(const Array<double,2>& x) = 0;    
  virtual void visit(const std::string& name,double value) = 0; 
  virtual void visit(const std::string& name,size_t value) = 0; 
  virtual void visit(const std::string& name,const Array<double,1>& value) = 0; 
  virtual void visit(const std::string& name,const Array<double,2>& value) = 0; 
};

}

#endif
