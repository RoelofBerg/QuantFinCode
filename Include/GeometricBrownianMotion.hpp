/** \file  GeometricBrownianMotion.hpp
    \brief Header file defining class representing a multidimensional 
           (e.g. multi-asset) geometric Brownian motion.
           Copyright 2006, 2007, 2013 by Erik Schlögl

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

#ifndef GEOMETRICBROWNIANMOTION_HPP
#define GEOMETRICBROWNIANMOTION_HPP
#include <stdexcept>
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include "BlackScholesAsset.hpp"

namespace quantfin { 

using blitz::Array;
using blitz::firstDim;
using blitz::secondDim;
using blitz::Range;
using blitz::toEnd;
using blitz::firstIndex;
using blitz::secondIndex;

class GeometricBrownianMotion { 
private:
  Array<double,1>*                                T; ///< Process timeline.
  Array<double,1>*                        timeline_; ///< Timeline of asset price realisations reported by operator().
  std::vector<const BlackScholesAsset*>& underlying; ///< Vector of pointers to underlying assets.
  Array<double,1>                                dW; ///< Work space to hold Brownian motion increments.
  Array<double,1>                           vol_lvl; ///< Work space to hold volatility levels.
  Array<double,2>                      asset_values; ///< Work space to hold asset values.
  Array<int,1>*                        time_mapping; ///< Mapping of indices from internal to external timeline.
  void add_drift(double tstart,double tend,int numeraire_index);
public:
  GeometricBrownianMotion(std::vector<const BlackScholesAsset*>& xunderlying);
  ~GeometricBrownianMotion();
  /// Query the dimension of the process.
  int dimension() const;
  /// Query the number of factors driving the process.
  inline int factors() const { return dW.extent(firstDim); };
  /// Set process timeline.
  bool set_timeline(const Array<double,1>& timeline);
  inline bool set_numeraire(int num) { return (num<dimension()); };
  /// Get process timeline.
  inline const Array<double,1>& get_timeline() const { return *timeline_; };
  inline int number_of_steps() const { return T->extent(firstDim)-1; };
  /// Generate a realisation of the process under the martingale measure associated with deterministic bond prices. underlying_values is an asset x (time points) Array.
  void operator()(Array<double,2>& underlying_values,const Array<double,2>& x,const TermStructure& ts);
  /// Generate a realisation of the process under the martingale measure associated with a given numeraire asset. underlying_values is an asset x (time points) Array.
  void operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const Array<double,2>& x,const TermStructure& ts,int numeraire_index);
};

}
 
#endif

