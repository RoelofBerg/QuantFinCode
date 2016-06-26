/** \file  GaussianEconomy.cpp
    \brief C++ source file implementing classes for representing a set of "economies" where asset prices are driving by multidimensional 
           geometric Brownian motion.
           Copyright 2010, 2011, 2012 by Erik Schlögl

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

#include <algorithm>
#include "GaussianEconomy.hpp"
#include "QFArrayUtil.hpp"

using namespace quantfin;

GaussianEconomy::GaussianEconomy(std::vector<boost::shared_ptr<BlackScholesAsset> >& xunderlying,
	                             boost::shared_ptr<DeterministicAssetVol> xv,
								 boost::shared_ptr<TermStructure> xinitialTS)
: underlying(xunderlying),v(xv),initialTS(xinitialTS),hjm(new GaussianHJM(&*xv,&*xinitialTS))
{
  int i;
  for (i=0;i<xv->factors();i++) component_vol.push_back(xv->component_vol(i));
}

/// Constructor
GaussMarkovWorld::GaussMarkovWorld(std::vector<boost::shared_ptr<GaussianEconomy> >& xeconomies,
	                               std::vector<boost::shared_ptr<DeterministicAssetVol> >& xvols,Array<double,1> xinitial_exchange_rates)
  : economies(xeconomies),vols(xvols),initial_exchange_rates(xinitial_exchange_rates.copy()),economy_start_index(xeconomies.size()),
    T(NULL),state_variables(NULL),state_variable_drifts(NULL),termstructure_statevariables(xeconomies[0]->v->factors()), // ,timeline_(NULL),time_mapping(NULL),dW(NULL)
	terminal_fwd_exchange_rate(xinitial_exchange_rates.extent(firstDim))
{ 
  initialise();
}

void GaussMarkovWorld::initialise()
{
  int i;
  current_numeraire = 0;
  if (!check_inputs()) throw std::logic_error("Non-conforming volatility functions in GaussMarkovWorld::GaussMarkovWorld()");
  diffusion_dimension       = economies[0]->v->factors();
  number_of_state_variables = 0;
  int max_number_of_assets  = 0;
  for (i=0;i<economies.size();i++) {
	economy_start_index(i)     = number_of_state_variables;
	max_number_of_assets       = std::max(max_number_of_assets,int(economies[i]->underlying.size()));
	number_of_state_variables += diffusion_dimension + economies[i]->underlying.size(); }
  currency_start_index       = number_of_state_variables;
  number_of_state_variables += vols.size();
  terminal_fwd_asset.resize(economies.size(),max_number_of_assets);
  zsum = z2sum = 0.0;
  zn = 0;
}

bool GaussMarkovWorld::is_asset_index(int i) const 
{ 
  int j;
  bool result = (i<number_of_state_variables-vols.size());
  if (result) {
	j = 0;
	while ((j<economies.size())&&(economy_start_index(j)<=i)) j++;
	result = (i>=economy_start_index(j-1)+diffusion_dimension); }
  return result;
}

int GaussMarkovWorld::economy_index(int j) const 
{  
  int result = -1;
  if (j<currency_start_index) {
	result = economies.size()-1;
	while (result>0) {
	  if (economy_start_index(result)<=j) return result;
	  result--; }}
  return result; 
}

bool GaussMarkovWorld::check_inputs()
{
  int i,j;
  bool result = true;
  int factors = economies[0]->v->factors();
  for (i=0;i<economies.size();i++) {
	for (j=0;j<economies[i]->underlying.size();j++) {
	  if (factors!=(economies[i]->underlying[j]->volatility_function()).factors()) {
	    result = false;
		break; }}
    if (!result) break; 
	if (factors!=economies[i]->v->factors()) {
	  result = false;
	  break; }}
  if (result) {
    for (i=0;i<vols.size();i++) {
	  if (factors!=vols[i]->factors()) {
	    result = false;
		break; }}}
  return result; 
  }

/// Destructor
GaussMarkovWorld::~GaussMarkovWorld()
{
  if (T) delete T;                                                  
  if (state_variables) delete state_variables;                                                  
  if (state_variable_drifts) delete state_variable_drifts;   
}

/** Determine which asset values should be reported back by operator()
    asset_index of -1 means fixed maturity zero coupon bond 
    asset_index of -2 means terminal forward exchange rate (expressed in units of domestic currency) 
    asset_index of -3 means fixed time-to-maturity zero coupon bond 
	Returns index of reportable asset added. */
int GaussMarkovWorld::set_reporting(int currency_index,int asset_index,double maturity)
{
  reportable r;
  r.currency_index = currency_index;
  r.asset_index    = asset_index;
  r.maturity       = maturity;
  reportable_list.push_back(r);
  return reportable_list.size()-1;
}

/// Query the dimension of the process.
int GaussMarkovWorld::dimension() const
{
  return reportable_list.size();
}

/// Generating a realisation of the process under the martingale measure associated with deterministic bond prices is not applicable in a stochastic term structure model.
void GaussMarkovWorld::operator()(Array<double,2>& underlying_values,const Array<double,2>& x,const TermStructure& ts)
{
  throw std::logic_error("GaussMarkovWorld: Deterministic zero coupon bond prices as numeraire are not applicable in a stochastic term structure model");
}

/** Generate a realisation of the process under the martingale measure associated with a given numeraire asset.
    numeraire_index = 0 means simulation under the domestic rolling spot measure.
	Negative numeraire_index means simulation under a foreign rolling spot measure, i.e. -i corresponds to i-th foreign currency. In this case
	(as in all other cases), the value of the numeraire is expressed in foreign currency.
	Positive numeraire_index means simulation under the martingale measure associated with taking the i-th risky asset in domestic currency as the numeraire.
	*/
void GaussMarkovWorld::operator()(Array<double,2>& underlying_values,Array<double,1>& numeraire_values,const Array<double,2>& x,const TermStructure& ts,int numeraire_index)
{
  int i,j;
  if (!T) throw std::logic_error("Timeline not set in GaussMarkovWorld::operator()");
  if (current_numeraire!=numeraire_index) set_numeraire(numeraire_index); // to ensure that the numeraire (if a risky asset) is in the reportable_list
  if ((x.extent(firstDim)<max_rank)||(x.extent(secondDim)<T->extent(firstDim)-1)) throw std::logic_error("Insufficient random variables passed to GaussMarkovWorld::operator()");
  // propagate state variables
  propagate_state_variables(x);
  // set return values
  for (i=0;i<reportable_list.size();i++) { // Loop over reportable assets
	if ((reportable_list[i].asset_index==-1)||(reportable_list[i].asset_index==-3)) { // zero coupon bond
	  const GaussianHJM& hjm = *(economies[reportable_list[i].currency_index]->hjm);
	  double mat = reportable_list[i].maturity;
	  int start_index = economy_start_index(reportable_list[i].currency_index);
      for (j=0;j<state_variables->extent(secondDim);j++) { // Loop over time line
	    termstructure_statevariables = (*state_variables)(Range(start_index,start_index+diffusion_dimension-1),j);
		double ttm = mat;
		if (reportable_list[i].asset_index==-1) ttm -= (*T)(j);
		underlying_values(i,j) = hjm.bond(termstructure_statevariables,(*T)(j),ttm); }}
	else {
	  if (reportable_list[i].asset_index==-2) { // terminal forward exchange rate
	    int idx = reportable_list[i].currency_index-1;
        for (j=0;j<state_variables->extent(secondDim);j++) { // Loop over time line
		  underlying_values(i,j) = terminal_fwd_exchange_rate(idx) * std::exp((*state_variables)(currency_start_index+idx,j)); }}
	  else {
	    int idx = reportable_list[i].currency_index;
	    int asset_index = economy_start_index(idx) + diffusion_dimension + reportable_list[i].asset_index;
        for (j=0;j<state_variables->extent(secondDim);j++) { // Loop over time line
		  underlying_values(i,j) = terminal_fwd_asset(idx,reportable_list[i].asset_index) * std::exp((*state_variables)(asset_index,j)); }}}}
  // set numeraire values
  if (numeraire_index<=0) { // rolling spot measure
    const GaussianHJM& hjm = *(economies[-numeraire_index]->hjm);
	int start_index = economy_start_index(-numeraire_index);
    numeraire_values(0) = 1.0;
	for (j=0;j<state_variables->extent(secondDim)-1;j++) { // Loop over time line
	  termstructure_statevariables = (*state_variables)(Range(start_index,start_index+diffusion_dimension-1),j);
      double short_bond = hjm.bond(termstructure_statevariables,(*T)(j),(*T)(j+1)-(*T)(j));
	  numeraire_values(j+1) = numeraire_values(j) / short_bond; }}
  else { // domestic risky asset as numeraire
    for (j=0;j<state_variables->extent(secondDim);j++) { // Loop over time line
	  // discount terminal forward to get current asset value
	  numeraire_values(j) = underlying_values(numeraire_reportable_index,j) * underlying_values(numeraire_reportable_ZCB,j); }}
}

/// Set process timeline.
bool GaussMarkovWorld::set_timeline(const Array<double,1>& timeline)
{
  int i,j,k,covar_i,covar_j;
  if (T) delete T;
  T = new Array<double,1>(timeline.copy());
  double T_N = (*T)(T->extent(firstDim)-1);
  Array<double,2> covar(number_of_state_variables,number_of_state_variables);
  const DeterministicAssetVol* row_v  = NULL;
  const DeterministicAssetVol* row_fv = NULL;
  const DeterministicAssetVol* col_v  = NULL;
  const DeterministicAssetVol* col_fv = NULL;
  max_rank = 0;
  // delete existing multivariate distributions of state variable increments
  j = mvn.size();
  for (i=0;i<j;i++) mvn.pop_back();
  int dim1 = number_of_state_variables-vols.size();
  const DeterministicAssetVol* domesticTS_vol = &*economies[0]->v;
  const DeterministicAssetVol* TS_vol  = NULL;
  const DeterministicAssetVol* TS_vol2 = NULL;
  for (i=0;i<timeline.extent(firstDim)-1;i++) {
    double dt = timeline(i+1) - timeline(i);
    // create covar matrix for each time step in T
	// ordering: interest rate & asset state variables for each economy first, then exchange rate state variables
	for (j=0;j<number_of_state_variables;j++) {
	  row_v = get_vol(j);
	  covar(j,j) = row_v->volproduct(timeline(i),dt,*row_v);
	  if (is_asset_index(j)) {
		TS_vol = &*economies[economy_index(j)]->v;
	    covar(j,j) -= 2.0*TS_vol->bondvolproduct(timeline(i),dt,T_N,*row_v);
	    covar(j,j) += TS_vol->bondbondvolproduct(timeline(i),dt,T_N,T_N,*TS_vol); }
	  if (j>=dim1) {
	    row_fv = &*economies[j-dim1+1]->v;
	    covar(j,j) += 2.0*row_fv->bondvolproduct(timeline(i),dt,T_N,*row_v);
	    covar(j,j) -= 2.0*domesticTS_vol->bondvolproduct(timeline(i),dt,T_N,*row_v);
	    covar(j,j) += row_fv->bondbondvolproduct(timeline(i),dt,T_N,T_N,*row_fv);
	    covar(j,j) -= 2.0*row_fv->bondbondvolproduct(timeline(i),dt,T_N,T_N,*domesticTS_vol);
	    covar(j,j) += domesticTS_vol->bondbondvolproduct(timeline(i),dt,T_N,T_N,*domesticTS_vol); }
	  for (k=0;k<j;k++) {
	    col_v = get_vol(k);
		covar(k,j) = row_v->volproduct(timeline(i),dt,*col_v); 
        if (is_asset_index(k)) {
		  TS_vol2 = &*economies[economy_index(k)]->v;
		  covar(k,j) -= TS_vol2->bondvolproduct(timeline(i),dt,T_N,*row_v); }
        if (is_asset_index(j)) {
	      covar(k,j) -= TS_vol->bondvolproduct(timeline(i),dt,T_N,*col_v);
          if (is_asset_index(k)) {
	        covar(k,j) += TS_vol2->bondbondvolproduct(timeline(i),dt,T_N,T_N,*TS_vol); }}
	    if (j>=dim1) {
	      covar(k,j) += row_fv->bondvolproduct(timeline(i),dt,T_N,*col_v);
	      covar(k,j) -= domesticTS_vol->bondvolproduct(timeline(i),dt,T_N,*col_v);
          if (is_asset_index(k)) {
	        covar(k,j) -= row_fv->bondbondvolproduct(timeline(i),dt,T_N,T_N,*TS_vol2); 
			covar(k,j) += domesticTS_vol->bondbondvolproduct(timeline(i),dt,T_N,T_N,*TS_vol2); }
	      if (k>=dim1) {
	        col_fv = &*economies[k-dim1+1]->v;
	        covar(k,j) += col_fv->bondvolproduct(timeline(i),dt,T_N,*row_v);
	        covar(k,j) -= domesticTS_vol->bondvolproduct(timeline(i),dt,T_N,*row_v);
	        covar(k,j) += row_fv->bondbondvolproduct(timeline(i),dt,T_N,T_N,*col_fv);
	        covar(k,j) -= domesticTS_vol->bondbondvolproduct(timeline(i),dt,T_N,T_N,*col_fv);
	        covar(k,j) -= row_fv->bondbondvolproduct(timeline(i),dt,T_N,T_N,*domesticTS_vol);
	        covar(k,j) += domesticTS_vol->bondbondvolproduct(timeline(i),dt,T_N,T_N,*domesticTS_vol); }}
		covar(j,k) = covar(k,j); }}
    // create MultivariateNormal from each covar matrix
	boost::shared_ptr<MultivariateNormal> mv(new MultivariateNormal(covar));
	mvn.push_back(mv);
    // required dimension of random variables on each time step is maximum rank of covar matrices
	max_rank = std::max(max_rank,mv->rank()); }
  // resize scratch arrays 
  if (!state_variables) state_variables = new Array<double,2>(number_of_state_variables,T->extent(firstDim));
  if (state_variables->extent(secondDim)!=T->extent(firstDim)) state_variables->resize(number_of_state_variables,T->extent(firstDim));
  for (i=0;i<number_of_state_variables;i++) (*state_variables)(i,0) = 0.0;
  if (!state_variable_drifts) state_variable_drifts = new Array<double,2>(number_of_state_variables,T->extent(firstDim)-1);
  if (state_variable_drifts->extent(secondDim)!=T->extent(firstDim)-1) state_variable_drifts->resize(number_of_state_variables,T->extent(firstDim)-1);
  set_numeraire(current_numeraire);
  // Set terminal forwards  
  for (i=0;i<terminal_fwd_exchange_rate.extent(firstDim);i++)
	terminal_fwd_exchange_rate(i) = initial_exchange_rates(i) * (*(economies[i+1]->initialTS))(T_N) / (*(economies[0]->initialTS))(T_N);
  for (i=0;i<economies.size();i++) {
    for (j=0;j<economies[i]->underlying.size();j++) {
	  terminal_fwd_asset(i,j)  = economies[i]->underlying[j]->initial_value();
	  // added next line to account for dividends
	  terminal_fwd_asset(i,j) *= economies[i]->underlying[j]->dividend_discount((*T)(0),T_N);
	  terminal_fwd_asset(i,j) /= (*(economies[i]->initialTS))(T_N); }}
  return T;   
}

const DeterministicAssetVol* GaussMarkovWorld::get_vol(int index)
{
  int i,economy_k;
  // exchange rate state variable?
  if (index>=number_of_state_variables-vols.size()) return &*vols[index-(number_of_state_variables-vols.size())];
  else {
    // which economy?
    i = index;
	economy_k = 0;
	while ((i>=0)&&(economy_k<economies.size())) {
	  if (i<diffusion_dimension) return &*economies[economy_k]->component_vol[i];
	  else {
	    i -= diffusion_dimension;
		if (i<economies[economy_k]->underlying.size()) return &(economies[economy_k]->underlying[i]->volatility_function());
		else {
		  i -= economies[economy_k]->underlying.size();
		  economy_k++; }}}}
  throw std::logic_error("Volatility function not found in GaussMarkovWorld::get_vol()");
  return NULL;
}

/** Set state variable drifts, etc. for simulation of the process under the martingale measure associated with a given numeraire asset.
    numeraire_index = 0 means simulation under the domestic rolling spot measure.
	Negative numeraire_index means simulation under a foreign rolling spot measure, i.e. -i corresponds to i-th foreign currency. In this case
	(as in all other cases), the value of the numeraire is expressed in foreign currency.
	Positive numeraire_index means simulation under the martingale measure associated with taking the i-th risky asset in domestic currency as the numeraire.
	*/
bool GaussMarkovWorld::set_numeraire(int numeraire_index)
{
  int i,j,k,base_j;
  double T_N = (*T)(T->extent(firstDim)-1);
  if (numeraire_index>(long)(economies[0]->underlying.size())) numeraire_index = 0;
  // base case: numeraire_index==0
  const DeterministicAssetVol* adjust_vol = NULL; // for deviations from the base case
  if (numeraire_index>0) adjust_vol = &(economies[0]->underlying[numeraire_index-1]->volatility_function());
  if (numeraire_index<0) adjust_vol = &*vols[-numeraire_index-1];
  Array<double,1> svm(economies[0]->v->factors());
  for (i=0;i<T->extent(firstDim)-1;i++) {
    double dt = (*T)(i+1) - (*T)(i);
    // domestic term structure state variables
	if (numeraire_index>=0) for (j=0;j<svm.extent(firstDim);j++) {  // expected value of the state variable increment under the discrete rolling spot measure
	  svm(j) = economies[0]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*economies[0]->component_vol[j]); }
	if (numeraire_index<0) {  // adjustment from foreign discrete rolling spot measure to foreign continuous spot measure
	  for (j=0;j<svm.extent(firstDim);j++) svm(j) = economies[-numeraire_index]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*economies[0]->component_vol[j]); }
	for (j=0;j<svm.extent(firstDim);j++) {
	  (*state_variable_drifts)(j,i) = svm(j); // if random variables are being generated under the discrete rolling spot measure, the expected value under this measure must be subtracted to convert to state variables defined as Ito integrals wrt Brownian motion under the continuous spot measure
	  if (adjust_vol) (*state_variable_drifts)(j,i) = economies[0]->component_vol[j]->volproduct((*T)(i),dt,*adjust_vol); } // adjust for domestic risky asset numeraire if applicable
    // domestic asset state variables
	for (j=0;j<economies[0]->underlying.size();j++) {
	  if (numeraire_index!=0) {  
		// adjustment from foreign discrete rolling spot measure to foreign continuous spot measure
	    if (numeraire_index<0) (*state_variable_drifts)(j+svm.extent(firstDim),i) = economies[-numeraire_index]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),economies[0]->underlying[j]->volatility_function());
		else                   (*state_variable_drifts)(j+svm.extent(firstDim),i) = 0.0;
		// from foreign continuous spot measure (or domestic risky asset measure) to domestic continuous spot measure
        if (adjust_vol) {
		  (*state_variable_drifts)(j+svm.extent(firstDim),i) += (economies[0]->underlying[j]->volatility_function()).volproduct((*T)(i),dt,*adjust_vol);
		  (*state_variable_drifts)(j+svm.extent(firstDim),i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,*adjust_vol); }
		// from domestic continuous spot measure to domestic terminal forward measure
		(*state_variable_drifts)(j+svm.extent(firstDim),i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,economies[0]->underlying[j]->volatility_function()); 
		(*state_variable_drifts)(j+svm.extent(firstDim),i) += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[0]->v)); }
	  else {
		// from domestic discrete rolling spot measure to domestic continuous spot measure
        (*state_variable_drifts)(j+svm.extent(firstDim),i) = economies[0]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),economies[0]->underlying[j]->volatility_function());
        (*state_variable_drifts)(j+svm.extent(firstDim),i) -= economies[0]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*(economies[0]->v));
		// from domestic continuous spot measure to domestic terminal forward measure
		(*state_variable_drifts)(j+svm.extent(firstDim),i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,economies[0]->underlying[j]->volatility_function()); 
		(*state_variable_drifts)(j+svm.extent(firstDim),i) += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[0]->v)); }
	  // additional drift terms not dependent on measure changes
	  double tmp = (economies[0]->underlying[j]->volatility_function()).volproduct((*T)(i),dt,economies[0]->underlying[j]->volatility_function());
	  tmp -= 2.0*economies[0]->v->bondvolproduct((*T)(i),dt,T_N,economies[0]->underlying[j]->volatility_function());
	  tmp += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*economies[0]->v);
	  (*state_variable_drifts)(j+svm.extent(firstDim),i) -= 0.5*tmp; }
	base_j = j + svm.extent(firstDim);
	// foreign economies
	for (k=1;k<economies.size();k++) {
	  // foreign term structure state variables
	  if (numeraire_index==-k) { // expected value of the state variable increment under the foreign discrete rolling spot measure
        for (j=0;j<svm.extent(firstDim);j++) {  // expected value of the state variable increment under the foreign discrete rolling spot measure
		  svm(j) = economies[k]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*economies[k]->component_vol[j]); }}	  
	  else {
	    if (numeraire_index<0) {  // needs fixing: adjustment from foreign discrete rolling spot measure to foreign continuous spot measure 
	      for (j=0;j<svm.extent(firstDim);j++) svm(j) = economies[-numeraire_index]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*economies[k]->component_vol[j]); }
		else { // adjustment from foreign continuous rolling spot measure to domestic discrete rolling spot measure - note that state variable expectation under foreign continuous rolling spot measure is zero
		  for (j=0;j<svm.extent(firstDim);j++) svm(j) = -economies[k]->component_vol[j]->volproduct((*T)(i),dt,*vols[k-1]); 
		  if (numeraire_index==0) for (j=0;j<svm.extent(firstDim);j++) svm(j) += economies[0]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*economies[k]->component_vol[j]); }}
	  for (j=0;j<svm.extent(firstDim);j++) {
	    (*state_variable_drifts)(base_j+j,i) = svm(j); // if random variables are being generated under the discrete rolling spot measure, the expected value under this measure must be subtracted to convert to state variables defined as Ito integrals wrt Brownian motion under the continuous spot measure
	    if (adjust_vol) (*state_variable_drifts)(base_j+j,i) += economies[k]->component_vol[j]->volproduct((*T)(i),dt,*adjust_vol); } // adjust for domestic risky asset numeraire if appicable
	  base_j += j;
	  // foreign asset state variables
	  for (j=0;j<economies[k]->underlying.size();j++) {
		(*state_variable_drifts)(j+base_j,i) = 0.0;
	    if (numeraire_index!=-k) {  
		  // adjustment from numeraire_index discrete rolling spot measure to numeraire_index continuous spot measure
	      if (numeraire_index<0) (*state_variable_drifts)(j+base_j,i) = economies[-numeraire_index]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),economies[k]->underlying[j]->volatility_function());
		  else                   (*state_variable_drifts)(j+base_j,i) = 0.0;
		  // from numeraire_index continuous spot measure (or domestic risky asset measure) to domestic continuous spot measure
          if (adjust_vol) {
			(*state_variable_drifts)(j+base_j,i) += (economies[k]->underlying[j]->volatility_function()).volproduct((*T)(i),dt,*adjust_vol);
		    (*state_variable_drifts)(j+base_j,i) -= economies[k]->v->bondvolproduct((*T)(i),dt,T_N,*adjust_vol); }
		  else { // from domestic discrete rolling spot measure to domestic continuous spot measure
            (*state_variable_drifts)(j+base_j,i) += economies[0]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),economies[k]->underlying[j]->volatility_function());
			(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*(economies[k]->v)); }
		  // from domestic continuous spot measure to k-th foreign continuous spot measure
          (*state_variable_drifts)(j+base_j,i) -= (economies[k]->underlying[j]->volatility_function()).volproduct((*T)(i),dt,*vols[k-1]);
          (*state_variable_drifts)(j+base_j,i) += economies[k]->v->bondvolproduct((*T)(i),dt,T_N,*vols[k-1]);
		  // from k-th foreign continuous spot measure to k-th foreign terminal forward measure
		  (*state_variable_drifts)(j+base_j,i) -= economies[k]->v->bondvolproduct((*T)(i),dt,T_N,economies[k]->underlying[j]->volatility_function()); 
		  (*state_variable_drifts)(j+base_j,i) += economies[k]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[k]->v)); }
	    else {
          // from k-th foreign discrete spot measure to k-th foreign terminal forward measure
		  (*state_variable_drifts)(j+base_j,i) = economies[k]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),economies[k]->underlying[j]->volatility_function());
		  (*state_variable_drifts)(j+base_j,i) -= economies[k]->v->bondvolproduct((*T)(i),dt,T_N,economies[k]->underlying[j]->volatility_function()); }
	    // additional drift terms not dependent on measure changes
	    double tmp = (economies[k]->underlying[j]->volatility_function()).volproduct((*T)(i),dt,economies[k]->underlying[j]->volatility_function());
	    tmp -= 2.0*economies[k]->v->bondvolproduct((*T)(i),dt,T_N,economies[k]->underlying[j]->volatility_function());
	    tmp += economies[k]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*economies[k]->v);
		(*state_variable_drifts)(j+base_j,i) -= 0.5*tmp; }
	  base_j += j; }
	// exchange rates
	for (j=0;j<vols.size();j++) {
	  if (numeraire_index!=0) {  
		// adjustment from foreign discrete rolling spot measure to foreign continuous spot measure
		if (numeraire_index<0) {
		  (*state_variable_drifts)(j+base_j,i) = economies[-numeraire_index]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*vols[j]);
		  (*state_variable_drifts)(j+base_j,i) += economies[-numeraire_index]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*economies[j+1]->v);
		  (*state_variable_drifts)(j+base_j,i) -= economies[-numeraire_index]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*economies[0]->v); }
		else (*state_variable_drifts)(j+base_j,i) = 0.0;
		// from foreign continuous spot measure (or domestic risky asset measure) to domestic continuous spot measure
		if (adjust_vol) {
		  (*state_variable_drifts)(j+base_j,i) += vols[j]->volproduct((*T)(i),dt,*adjust_vol);
		  (*state_variable_drifts)(j+base_j,i) += economies[j+1]->v->bondvolproduct((*T)(i),dt,T_N,*adjust_vol);
		  (*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,*adjust_vol); }
		// from domestic continuous spot measure to domestic terminal forward measure
		(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,*vols[j]);
		(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[j+1]->v));
		(*state_variable_drifts)(j+base_j,i) += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[0]->v)); }
	  else {
        // from domestic discrete spot measure to domestic terminal forward measure
		(*state_variable_drifts)(j+base_j,i) = economies[0]->v->bondvolproduct((*T)(i),dt,(*T)(i+1),*vols[j]);
		(*state_variable_drifts)(j+base_j,i) += economies[0]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*(economies[j+1]->v));
		(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondbondvolproduct((*T)(i),dt,(*T)(i+1),T_N,*(economies[0]->v));
		(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondvolproduct((*T)(i),dt,T_N,*vols[j]);
		(*state_variable_drifts)(j+base_j,i) -= economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[j+1]->v));
		(*state_variable_drifts)(j+base_j,i) += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*(economies[0]->v)); }
	  // additional drift terms not dependent on measure changes
	  double tmp = vols[j]->volproduct((*T)(i),dt,*vols[j]);
	  tmp += 2.0*economies[j+1]->v->bondvolproduct((*T)(i),dt,T_N,*vols[j]);
	  tmp -= 2.0*economies[0]->v->bondvolproduct((*T)(i),dt,T_N,*vols[j]);
	  tmp += economies[j+1]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*economies[j+1]->v);
	  tmp -= 2.0*economies[j+1]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*economies[0]->v);
	  tmp += economies[0]->v->bondbondvolproduct((*T)(i),dt,T_N,T_N,*economies[0]->v);
	  (*state_variable_drifts)(j+base_j,i) -= 0.5*tmp; }}
  current_numeraire = numeraire_index;
  // Ensure that the numeraire (if a risky asset) and the associated terminal zero coupon bond is in the reportable_list
  if (numeraire_index>0)  {
    bool numeraire_found     = false;
    bool numeraire_ZCB_found = false;
	for (i=0;i<reportable_list.size();i++) {
	  if (reportable_list[i].currency_index==0) {
		if ((reportable_list[i].asset_index==numeraire_index-1)) {
	      numeraire_found = true;
		  numeraire_reportable_index = i;
		  if (numeraire_ZCB_found) break; }
	    if ((reportable_list[i].asset_index==-1)&&(reportable_list[i].maturity==T_N)) {
	      numeraire_ZCB_found = true;
		  numeraire_reportable_ZCB = i;
		  if (numeraire_found) break; }}}
	if (!numeraire_found) numeraire_reportable_index = set_reporting(0,numeraire_index-1,0.0); 
	if (!numeraire_ZCB_found) numeraire_reportable_ZCB = set_reporting(0,-1,T_N); }
  return true;
}

void GaussMarkovWorld::propagate_state_variables(const Array<double,2>& x)
{
  int j;
  Array<double,1> current_x(number_of_state_variables);
  // propagate state variables
  for (j=1;j<state_variables->extent(secondDim);j++) { // Loop over time line
	// ordering: interest rate & asset state variables for each economy first, then exchange rate state variables
    // term structure state variables are defined under the native currency spot measure
    // all asset and exchange rate state variables are defined under the domestic terminal measure
    current_x(Range(0,max_rank-1)) = x(Range(0,max_rank-1),j-1);
    mvn[j-1]->transform_random_variables(current_x);
	(*state_variables)(Range::all(),j) = (*state_variables)(Range::all(),j-1) + (*state_variable_drifts)(Range::all(),j-1) + current_x; }
  zsum += (*state_variables)(0,1);
  z2sum += (*state_variables)(0,1)*(*state_variables)(0,1);
  zn++;
}

/** Generate a realisation of the process under the martingale measure associated with a given numeraire asset, 
    returning the state variables. underlying_values is an asset x (time points) Array. */
const Array<double,2>& GaussMarkovWorld::generate_state_variables(Array<double,1>& numeraire_values,const Array<double,2>& x,int numeraire_index) 
{
  int j;
  if (!T) throw std::logic_error("Timeline not set in GaussMarkovWorld::generate_state_variables");
  if (numeraire_index>0) throw std::logic_error("Cannot generate risky asset as numeraire in GaussMarkovWorld::generate_state_variables. Use operator() instead.");
  if (current_numeraire!=numeraire_index) set_numeraire(numeraire_index); // to ensure that the numeraire (if a risky asset) is in the reportable_list
  if ((x.extent(firstDim)<max_rank)||(x.extent(secondDim)<T->extent(firstDim)-1)) throw std::logic_error("Insufficient random variables passed to GaussMarkovWorld::operator()");
  // propagate state variables
  propagate_state_variables(x);
  // set numeraire values
  const GaussianHJM& hjm = *(economies[-numeraire_index]->hjm);
  int start_index = economy_start_index(-numeraire_index);
  numeraire_values(0) = 1.0;
  for (j=0;j<state_variables->extent(secondDim)-1;j++) { // Loop over time line
	termstructure_statevariables = (*state_variables)(Range(start_index,start_index+diffusion_dimension-1),j);
    double short_bond = hjm.bond(termstructure_statevariables,(*T)(j),(*T)(j+1));
	numeraire_values(j+1) = numeraire_values(j) / short_bond; }
  return *state_variables;
}

GaussMarkovTermStructure GaussMarkovWorld::get_TermStructure(int i,int j) 
{
  const GaussianHJM& hjm = *(economies[j]->hjm);
  int start_index = economy_start_index(j);
  termstructure_statevariables = (*state_variables)(Range(start_index,start_index+diffusion_dimension-1),i);
  GaussMarkovTermStructure result((*T)(i),hjm,termstructure_statevariables);
  return result;
}
