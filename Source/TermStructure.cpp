/** \file  TermStructure.cpp
    \brief C++ source file implementing classes for term structures.
           Copyright 2003, 2005, 2010, 2011, 2013 by Erik Schlögl

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

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "DeterministicCashflow.hpp"
#include "TSBootstrap.hpp"
#include "Rootsearch.hpp"
#include "Powell.hpp"
#include "QFUtil.hpp"
#include "Linesearch.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

/// Flat term structure constructor.
FlatTermStructure::FlatTermStructure(double lvl,    ///< Flat level of continuously compounded interest rates.
                                     double tzero,  ///< Start of maturity time line.
                                     double horizon ///< End of maturity time line.
                                     ) : TermStructure(2),level(lvl)
{
  T(0) = tzero;
  T(1) = horizon;
  B(0) = 1.0;
  B(1) = std::exp(-lvl*(horizon-tzero));                            
}
                       
int TermStructure::find_segment(double xi) const
{
  if (xi<T(0)) throw std::out_of_range("cannot extrapolate");
  for (int i=1;i<T.extent(firstDim);i++) if (xi<=T(i)) return i-1;
  throw std::out_of_range("cannot extrapolate");
}

void TermStructure::approximate_fit(std::vector<DeterministicCashflow> cashflows,double eps)
{
  fit_func fit_f(cashflows,*this); 
  opt::GeneralBrentLinesearch<fit_func,double> ls;
  Array<double,1> currpos = B(blitz::Range(1,blitz::toEnd)).copy();
  currpos = inverse_logistic_transform(currpos,0.0,1.0);
  if (currpos.extent(firstDim)>1) {
    opt::GeneralPowell<fit_func,double,opt::GeneralBrentLinesearch<fit_func,double> > powell(fit_f,eps,ls,5000);
    powell.solve(currpos); }
  else {
    Array<double,1> direction(1);
    direction(0) = 1.0;
    ls(fit_f,currpos,direction,eps,5000); }
}

double TermStructure::fit_func::operator()(Array<double,1>& lnzcb)
{
  int i;
  ts.B(blitz::Range(1,blitz::toEnd)) = logistic_transform(lnzcb,0.0,1.0);
  ts.reinitialise();     
  double result = 0.0;
  for (i=0;i<(int)cashflows.size();i++) {
    double err = cashflows[i].market_value() - cashflows[i].NPV(ts);
    result += err*err; }
  return result;
}

void TermStructure::reinitialise()
{ }

void FlatTermStructure::reinitialise()
{ 
  level = -log(B(1))/T(1);   
}

Array<double,1> TermStructure::simple_rate(const Array<double,1>& xT,double delta) const
{
  Array<double,1> erg(xT.extent(firstDim)-1);
  double rdelta = 1.0/delta;
  for (int i=0;i<erg.extent(firstDim);i++) 
    erg(i) = rdelta * (operator()(xT(i))/operator()(xT(i+1)) - 1.0);
  return erg;
}

double TermStructure::simple_rate(double xT,double delta) const
{
  double erg(1.0/delta);
  erg *= operator()(xT)/operator()(xT+delta) - 1.0;
  return erg;
}

Array<double,1> TermStructure::swaps(const Array<double,1>& xT,double delta) const
{
  int n = xT.extent(firstDim)-1;
  Array<double,1> erg(n);
  int i = n - 1;
  erg(i) = operator()(xT(i+1));
  i--;
  while (i>=0) {
    erg(i) = erg(i+1) + operator()(xT(i+1));
    i--; }
  erg *= delta;
  for (int j=0;j<n;j++) erg(j) = (operator()(xT(j)) - operator()(xT(n)))/erg(j);
  return erg;
}

double TermStructure::swap(double T0,int n,double delta) const
{
  double erg = pvbp(T0,n,delta);
  return ((operator()(T0) - operator()(T0+n*delta)) / (delta * erg));
}

double TermStructure::swap(const Array<double,1>& tenor) const
{
  return (operator()(tenor(0)) - operator()(tenor(tenor.extent(firstDim)-1))) / pvbp(tenor);
}

double TermStructure::pvbp(double T0,int n,double delta) const
{
  double erg = 0.0;
  for (int i=1;i<=n;i++) erg += operator()(T0+i*delta);
  return erg;
}

double TermStructure::pvbp(const Array<double,1>& tenor) const
{
  double erg = 0.0;
  for (int i=1;i<tenor.extent(firstDim);i++) erg += (tenor(i)-tenor(i-1))*operator()(tenor(i));
  return erg;
}

Array<double,1> TermStructure::ccfwd(const Array<double,1>& xT) const
{
  int i;
  Array<double,1> tmpB(xT.extent(firstDim));
  for (i=0;i<tmpB.extent(firstDim);i++) tmpB(i) = operator()(xT(i));
  Array<double,1> F(xT.extent(firstDim)-1);
  for (i=0;i<F.extent(firstDim);i++) F(i) = std::log(tmpB(i)/tmpB(i+1))/(xT(i+1)-xT(i));
  return F;
}

double TermStructure::find_maturity(double target) const
{
  int n = T.extent(firstDim);
  if ((target>B(0))||(target<B(n-1))) return -1.0;  // Maturity not found
  boost::function<double (double)> objective_function = boost::bind(&TermStructure::operator(),this,_1);
  Rootsearch<boost::function<double (double)>,double,double> rs(objective_function,target,0.5*(T(0)+T(n-1)),0.25*(T(n-1)-T(0)),T(0),T(n-1));
  double result = rs.solve();
  return result;
}

FlatTermStructure* FlatTermStructure::pointer_to_copy() const
{
  return new FlatTermStructure(level,T(0),T(T.extent(firstDim)-1));
}

/** Function returning the interpolated value, which for term structures of
    interest rates is the maturity t, time T(0) forward zero coupon bond price. */
double FlatTermStructure::operator()(double t) const
{
  if (t<T(0)) throw std::logic_error("Requested maturity before beginning of term structure");
  return std::exp(-level*(t-T(0)));      
}    

TSLogLinear* TSLogLinear::pointer_to_copy() const
{
  return new TSLogLinear(T,B);
}

/** Function returning the interpolated value, which for term structures of
    interest rates is the maturity t, time T(0) forward zero coupon bond price. */
double TSLogLinear::operator()(double xi) const
{
  int idx = find_segment(xi);
  if (xi==T(idx+1)) return B(idx+1);
  else return B(idx) * std::exp((xi-T(idx))/(T(idx+1)-T(idx))*std::log(B(idx+1)/B(idx)));
}

TSLinear* TSLinear::pointer_to_copy() const
{
  return new TSLinear(T,B);
}

/** Function returning the interpolated value, which for term structures of
    interest rates is the maturity t, time T(0) forward zero coupon bond price. */
double TSLinear::operator()(double xi) const
{
  int   idx = find_segment(xi);
  double alpha = (xi-T(idx))/(T(idx+1)-T(idx));
  return (1.0-alpha)*B(idx) + alpha*B(idx+1);
}

/** Bootstrap term structure from market value of deterministic cash flows.
    Note that the vector of cashflows needs to be passed by value, so that we have a copy, which we can sort. */
void TSBootstrap::bootstrap(std::vector<DeterministicCashflow> cashflows,double eps)
{
  int i,j;
  std::sort(cashflows.begin(),cashflows.end(),DeterministicCashflow::shorter_cashflow());
  const Array<double,1>& t = cashflows[0].timeline();
  T.resize((int)cashflows.size()+1);
  T(0) = t(0);
  T(1) = t(t.extent(firstDim)-1);
  B.resize(T.extent(firstDim));
  B(0) = 1.0;
  for (i=1;i<(int)cashflows.size();i++) {
    const Array<double,1>& t = cashflows[i].timeline();
    if (!verify_compare(T(0),t(0),eps)) throw std::logic_error("Initial date mismatch"); 
    T(i+1) = t(t.extent(firstDim)-1); }
  for (i=0;i<(int)cashflows.size();i++) {
    const Array<double,1>& t = cashflows[i].timeline();
    const Array<double,1>& payments = cashflows[i].cashflow();
    double target_npv = cashflows[i].market_value();
    for (j=1;t(j)<=T(i)+eps;j++) target_npv -= payments(j-1)*operator()(t(j));
    if (j>=t.extent(firstDim)) throw std::logic_error("Too much overlap");
    blitz::Range s(j,t.extent(firstDim)-1);
    bootstrap_class b(t,payments,s,i+1,this);
    Rootsearch<bootstrap_class,double,double> rs(b,target_npv,B(i),0.1*B(i),0.0);
    rs.solve();     }
}

double TSBootstrap::bootstrap_class::operator()(double b)
{
  int i;
  ts_->B(idx) = b;
  double npv = 0.0;
  for (i=0;i<slice.length();i++) {
    npv += p_(slice(i)-1) * ts_->operator()(t_(slice(i))); }
  return npv;
}

/// For nested Visitor pattern.
void TermStructure::accept(QFNestedVisitor& visitor) const
{
  visitor.visit("Timeline",T);    
  visitor.visit("ZCB",B);    
}

/// For nested Visitor pattern.
void FlatTermStructure::accept(QFNestedVisitor& visitor) const
{
  TermStructure::accept(visitor);    
  visitor.visit("RateLevel",level);    
}

const std::string& FlatTermStructure::name() const
{
  static const std::string str("FlatTermStructure");
  return str;      
}

const std::string& TSLinear::name() const
{
  static const std::string str("TSLinear");
  return str;      
}

const std::string& TSLogLinear::name() const
{
  static const std::string str("TSLogLinear");
  return str;      
}

TSLinear::TSLinear(const Array<double,1>& xT,const TermStructure& xts)
	: TSBootstrap(xT.copy(),xT.copy())
{
  int i;
  for (i=0;i<xT.extent(firstDim);i++) B(i) = xts(xT(i));
}