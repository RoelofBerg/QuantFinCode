/** \file  TSBinomial.hpp
    \brief Header file defining classes for interest rate option pricing using binomial trees.
           Copyright 2006, 2015 by Erik Schlögl

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

#ifndef TSBINOMIAL_HPP
#define TSBINOMIAL_HPP

#include <stdexcept>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <boost/function.hpp>
#include "TSPayoff.hpp"
#include "TSBootstrap.hpp"
#include "TSInstruments.hpp"

namespace quantfin { 

using blitz::Array;
using blitz::Range;
using blitz::firstIndex;
using blitz::secondIndex;

class TSBinomialLatticeNode {
private:  
  double               Q; ///< State price of this node as viewed from the initial node.
  double               r; ///< Short rate realisation in this node.
public:
  TSLogLinear*        ts; ///< Term structure realisation in this node.
  inline TSBinomialLatticeNode();
  inline double state_price() const { return Q; };  
  inline double& state_price() { return Q; };  
  inline double short_rate() const { return r; };  
  inline double& short_rate() { return r; };  
};

inline TSBinomialLatticeNode::TSBinomialLatticeNode()
  : r(0.0),Q(0.0),ts(NULL)
{ }

/// Binomial interest rate term structure model following Black, Derman and Toy (1990) and Sandmann and Sondermann (1989, 1991), implemented as described in Section 3.5.2 of Schlögl (2013)
class TSBinomialMethod {
private:
  const TermStructure& initial_ts;
  Array<double,1>        T; ///< Time line.
  Array<double,1>       dt; ///< Period lengths.
  int                    n; ///< Number of time steps.
  Array<double,1>      sgm; ///< Volatility parameter.
  double                 p; ///< Probability of an up move.
  Array<double,1>  rfactor;
  int       current_period;
  std::vector<std::vector<TSBinomialLatticeNode*>*> nodes;
  Array<double,1> gridslice; ///< The current time slice of the finite difference grid.
  Array<double,1> tmp;
  inline void createLattice() { createLattice(0,n); };      
  void createLattice(int from,int to);      
  void allocateLattice();
  inline const TSBinomialLatticeNode* node(int i,int j) const { return (*nodes[i])[j]; };
  inline TSBinomialLatticeNode* node(int i,int j) { return (*nodes[i])[j]; };
  double bond_r(double r);
  inline void set_current_period(int i) { current_period = i; };
  // For volatility calibration.
  TSEuropeanInstrument* curr_instrument;
  int prev_idx,idx;
  inline void set_current_instrument(TSEuropeanInstrument& instrument,int prev_i,int i)
    { curr_instrument = &instrument; prev_idx = prev_i; idx = i; };
  double price_current_instrument(double try_sgm);
  void resetTermstructures(int i);
public:
  /// Constant volatility constructor.
  TSBinomialMethod(const TermStructure& xinitial_ts,double xsgm,const Array<double,1>& xT);
  /// Constructor with time dependent volatility.
  TSBinomialMethod(const TermStructure& xinitial_ts,const Array<double,1>& xsgm,const Array<double,1>& xT);
  ~TSBinomialMethod();
  /// Calibrate to a given set of caplets or other compatible instruments. Returns remaining pricing error.
  double calibrate(std::vector<TSEuropeanInstrument*> instruments);
  /// Create full term structures in each node by backward iteration through the lattice.
  void rollbackTermstructures();
  /// Test function for lattice.
  bool verify() const;
  /// Apply payoff f in period i
  inline void apply_payoff(int i,boost::function<double (const TermStructure&)> f);    
  void rollback(int from,int to);
  /// Rollback all the way using state prices.
  void rollback(int from);
  /// Rollback with early exercise (or similar) condition.
  void rollback(int from,int to,boost::function<double (double,const TermStructure&)> f);
  /// Access the calculated (rolled-back) price in node 0;
  inline double result() const { return gridslice(0); };
  /// Price a European payoff
  double price(TSEuropeanInstrument& instrument);
  /// Price a Bermudan (or similarly path-dependent) payoff.
  double price(TSBermudanInstrument& instrument);
  /// Access short rates
  inline double short_rate(int time,int state) const { return node(time,state)->short_rate(); };
  /// Access state prices
  inline double state_price(int time,int state) const { return node(time,state)->state_price(); };
};

inline void TSBinomialMethod::apply_payoff(int i,boost::function<double (const TermStructure&)> f) 
{
  gridslice.resize(nodes[i]->size());
  std::vector<TSBinomialLatticeNode*>::iterator j;
  int jidx = 0;
  for (j=nodes[i]->begin();j!=nodes[i]->end();j++) {
    if (!((*j)->ts)) throw std::logic_error("Term structure not available in this node");
    gridslice(jidx++) = f(*((*j)->ts)); }
}

}

#endif
