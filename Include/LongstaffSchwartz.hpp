/** \file  LongstaffSchwartz.hpp
    \brief Header file declaring template for the American/Bermudan option algorithm of Longstaff/Schwartz.
           Copyright 2006, 2010, 2013 by Erik Schlögl

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
 
#ifndef LONGSTAFFSCHWARTZ_HPP
#define LONGSTAFFSCHWARTZ_HPP
#include <vector>
#include "MSWarnings.hpp"
#include <blitz/array.h>
#include <boost/function.hpp>
#include "Regression.hpp"
#include "TermStructure.hpp"

namespace quantfin {

using blitz::Array;
using blitz::Range;
using blitz::fromStart;
using blitz::firstDim;
using blitz::secondDim;
using blitz::thirdDim;
using blitz::firstIndex;
using blitz::secondIndex;

/// Longstaff/Schwartz exercise boundary, single state variable, deterministic discounting, polynomial least squares fit.
class LongstaffSchwartzExerciseBoundary1D {
private:
  std::vector<PolynomialLeastSquaresFit*> fit;  ///< Exercise boundary functions (one for each point on the time line)
  Array<double,1>                           T;  ///< Time line
  const TermStructure&                     ts;
  boost::function<double (double)>          f;
  bool                          positive_only;  ///< True if option should be exercised only if the payoff is positive.
public:
  /// Constructor using polynomial regression.
  LongstaffSchwartzExerciseBoundary1D(const Array<double,1>& timeline,
                                      const TermStructure& xts,
                                      const Array<double,2>& state_variable_paths,
                                      boost::function<double (double)> payoff,
                                      int degree,
                                      bool positive_payoff_only = true);
  /// Destructor.
  ~LongstaffSchwartzExerciseBoundary1D();
  /// Apply the exercise to a Monte Carlo path.
  double apply(const Array<double,1>& path) const;
};

/** Longstaff/Schwartz exercise boundary, arbitrary number of state variables, discounting by numeraire, 
    arbitrary set of basis functions for the least squares fit.
	*/
class LongstaffSchwartzExerciseBoundary {
private:
  const std::vector<boost::function<double (double,const Array<double,1>&)> >& basis_functions; ///< Vector of basis functions for the regression.
  std::vector<LeastSquaresFit*>                         fit; ///< Exercise boundary functions (one for each point on the time line)
  Array<double,1>                                         T; ///< Time line
  boost::function<double (double,const Array<double,1>&)> f; ///< Early exercise payoff function taking time and state variables as arguments.
  bool                                        positive_only; ///< True if option should be exercised only if the payoff is positive.
  int                                   number_of_variables;
public:
  /** Constructor. state_variable_paths is a (paths,timepoints,state variables) array of state variable realisations. 
      numeraire_values is a (paths,timepoints) array of numeraire value realisations. */
  LongstaffSchwartzExerciseBoundary(const Array<double,1>& timeline,
                                    const Array<double,3>& state_variable_paths,
									const Array<double,2>& numeraire_values,
                                    boost::function<double (double,const Array<double,1>&)> payoff,
                                    const std::vector<boost::function<double (double,const Array<double,1>&)> >& xbasis_functions,
                                    bool positive_payoff_only = true);
  /// Destructor.
  ~LongstaffSchwartzExerciseBoundary();
  /// Apply the exercise to a Monte Carlo path.
  double apply(const Array<double,2>& path,const Array<double,1>& numeraire_values) const;
  inline int number_of_state_variables() const { return number_of_variables; };
  inline const Array<double,1>& get_timeline() const { return T; };
};

void add_polynomial_basis_function(std::vector<boost::function<double (double,const Array<double,1>&)> >& basisfunctions,const Array<int,1>& p);
void add_polynomial_basis_function(std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> >& basisfunctions,const Array<int,1>& p);

/** Longstaff/Schwartz exercise boundary, arbitrary number of state variables (including path history), discounting by numeraire, 
    arbitrary set of basis functions for the least squares fit. Basis functions take as arguments observation time points and
	observed state variable values, where the latter is a (state variables) X (observation time points) Array.
	*/
class RegressionExerciseBoundary {
protected:
  const std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> >& basis_functions; ///< Vector of basis functions for the regression.
  std::vector<LeastSquaresFit*>                         fit; ///< Exercise boundary functions (one for each point on the time line)
  Array<double,1>                                         T; ///< Time line
  boost::function<double (const Array<double,1>&,const Array<double,2>&)> f; ///< Early exercise payoff function (not discounted!) taking time and state variable history as arguments.
  bool                                        positive_only; ///< True if option should be exercised only if the payoff is positive.
  int                                   number_of_variables;
  class BasisFunction {
  private:
    const boost::function<double (const Array<double,1>&,const Array<double,2>&)> f;
	Array<double,1>                                                               T;
  public:
	inline BasisFunction(const boost::function<double (const Array<double,1>&,const Array<double,2>&)> xf,const Array<double,1>& xT)
	  : f(xf),T(xT) { };
	double operator()(const Array<double,1>& x) const;
  };
  std::vector<BasisFunction*> basis; ///< Vector of BasisFunction objects for LeastSquaresFit
  void initialise(const Array<double,3>& state_variable_paths,const Array<double,2>& numeraire_values,boost::function<double (const Array<double,1>&,const Array<double,2>&)> payoff,bool positive_payoff_only);
public:
  /** Constructor. state_variable_paths is a (paths,timepoints,state variables) array of state variable realisations. 
      numeraire_values is a (paths,timepoints) array of numeraire value realisations. */
  RegressionExerciseBoundary(const Array<double,1>& timeline,
                             const Array<double,3>& state_variable_paths,
							 const Array<double,2>& numeraire_values,
                             boost::function<double (const Array<double,1>&,const Array<double,2>&)> payoff,
                             const std::vector<boost::function<double (const Array<double,1>&,const Array<double,2>&)> >& xbasis_functions,
                             bool positive_payoff_only = true,
							 bool ini = true);
  /// Destructor.
  virtual ~RegressionExerciseBoundary();
  /// Apply the exercise to a Monte Carlo path.
  virtual double apply(const Array<double,2>& path,const Array<double,1>& numeraire_values);
  inline int number_of_state_variables() const { return number_of_variables; };
  inline const Array<double,1>& get_timeline() const { return T; };
};

double LSArrayAdapter(double t,const Array<double,1>& x,boost::function<double (double)> f,int idx);
double LSArrayAdapterT(double t,const Array<double,1>& x,boost::function<double (double,double)> f,int idx);
double REBAdapter(const Array<double,1>& T,const Array<double,2>& x,boost::function<double (double)> f,int idx); 
double REBAdapterT(const Array<double,1>& T,const Array<double,2>& x,boost::function<double (double,double)> f,int idx); 

}

#endif
