/** \file  CSVconstructors.cpp
    \brief C++ source file implementing class constructors taking a CSV file name as their argument.
           Copyright 2012, 2013 by Erik Schlögl

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

#include "GaussianEconomy.hpp"
#include "CSV2Array.hpp"
#include "ExponentialVol.hpp"
#include "PiecewiseVol.hpp"

using namespace quantfin;

boost::shared_ptr<TermStructure> CSVreadtermstructure(const char* path)
{
  int i;
  std::ifstream is_inputs(path);
  if (!is_inputs.is_open()) throw std::logic_error("Failed to open input file in CSVreadtermstructure()");
  blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
  std::map<std::string,std::string> inputs_map;
  for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
  if (!inputs_map.count("Term structure type")) throw std::logic_error("Missing term structure type");
  boost::shared_ptr<TermStructure> result;
  if (inputs_map["Term structure type"]=="flat") {
	double lvl;
	std::string str("Interest rate level");
	if (inputs_map.count(str)) {
	  double start = 0.0;
	  double end   = 40.0;
	  lvl = std::atof(inputs_map[str].data());
	  str = "Start time";
	  if (inputs_map.count(str)) start = std::atof(inputs_map[str].data());
	  str = "End time";
	  if (inputs_map.count(str)) end = std::atof(inputs_map[str].data());
	  result.reset(new FlatTermStructure(lvl,start,end));
	  return result; }
	else throw std::logic_error("Missing interest rate level specification"); }
}

boost::shared_ptr<DeterministicAssetVol> CSVreadvolatility(const char* path)
{
  int i;
  std::ifstream is_inputs(path);
  if (!is_inputs.is_open()) throw std::logic_error("Failed to open input file in CSVreadvolatility()");
  blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
  std::map<std::string,std::string> inputs_map;
  for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
  int voldim = 0;
  if (inputs_map.count("Volatility dimension")) voldim = std::atoi(inputs_map["Volatility dimension"].data());
  if (voldim<1) throw std::logic_error("Misspecified volatility dimension");
  if (!inputs_map.count("Volatility function type")) throw std::logic_error("Missing volatility function type");
  Array<double,1> vollvl(voldim);
  if (inputs_map["Volatility function type"]=="constant") {
	for (i=0;i<voldim;i++) {
	  std::string str("Volatility level ");
	  str += (char)(48+i);
	  if (inputs_map.count(str)) vollvl(i) = std::atof(inputs_map[str].data());
	  else throw std::logic_error("Missing volatility level specification"); }
	boost::shared_ptr<ConstVol> vol(new ConstVol(vollvl));
	return vol; }
  if (inputs_map["Volatility function type"]=="exponential") {
    Array<double,1> mr(voldim);
	for (i=0;i<voldim;i++) {
	  std::string str("Volatility level ");
	  str += (char)(48+i);
	  if (inputs_map.count(str)) vollvl(i) = std::atof(inputs_map[str].data());
	  else throw std::logic_error("Missing volatility level specification"); 
	  std::string str1("Mean reversion speed ");
	  str1 += (char)(48+i);
	  if (inputs_map.count(str1)) mr(i) = std::atof(inputs_map[str1].data());
	  else throw std::logic_error("Missing mean reversion speed specification"); }
	boost::shared_ptr<ExponentialVol> vol(new ExponentialVol(vollvl,mr));
	return vol; }
  throw std::logic_error("Unsupported volatility function type");
}

BlackScholesAsset::BlackScholesAsset(const char* path)
{
  int i;
  std::ifstream is_inputs(path);
  if (!is_inputs.is_open()) throw std::logic_error("Failed to open input file in BlackScholesAsset constructor");
  blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
  std::map<std::string,std::string> inputs_map;
  for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
  if (inputs_map.count("CSV file for volatility")) sv = CSVreadvolatility(inputs_map["CSV file for volatility"].data());
  else throw std::logic_error("Missing volatility specification");
  v = &*sv;
  std::string str("Initial value");
  if (inputs_map.count(str)) xzero = std::atof(inputs_map[str].data());
  else throw std::logic_error("Missing initial value for BlackScholesAsset"); 
  if (inputs_map.count("CSV file for dividend yield term structure")) 
	dividend = CSVreadtermstructure(inputs_map["CSV file for dividend yield term structure"].data());
  else dividend.reset(new FlatTermStructure(0.0,0.0,100.0));
}

GaussianEconomy::GaussianEconomy(const char* path)
{
  int i;
  std::ifstream is_inputs(path);
  if (!is_inputs.is_open()) throw std::logic_error("Failed to open input file in GaussianEconomy constructor");
  blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
  std::map<std::string,std::string> inputs_map;
  for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
  int number_of_underlyings = -1;
  if (inputs_map.count("Number of Black/Scholes assets")) number_of_underlyings = std::atoi(inputs_map["Number of Black/Scholes assets"].data());
  if (number_of_underlyings<0) throw std::logic_error("Misspecified number of Black/Scholes assets");
  boost::shared_ptr<BlackScholesAsset> asset;
  for (i=0;i<number_of_underlyings;i++) {
	std::string str("CSV file for Black/Scholes asset ");
	str += (char)(48+i);
	if (inputs_map.count(str)) {
	  asset.reset(new BlackScholesAsset(inputs_map[str].data()));
	  underlying.push_back(asset); }
	else throw std::logic_error("Missing Black/Scholes asset specification"); }
  if (inputs_map.count("CSV file for interest rate volatility")) v = CSVreadvolatility(inputs_map["CSV file for interest rate volatility"].data());
  else throw std::logic_error("Missing interest rate volatility specification");
  for (i=0;i<v->factors();i++) component_vol.push_back(v->component_vol(i));
  if (inputs_map.count("CSV file for initial interest rate term structure")) 
	initialTS = CSVreadtermstructure(inputs_map["CSV file for initial interest rate term structure"].data());
  else throw std::logic_error("Missing initial interest rate term structure");
  hjm.reset(new GaussianHJM(&*v,&*initialTS));
}

GaussMarkovWorld::GaussMarkovWorld(const char* path)
  : T(NULL),state_variables(NULL),state_variable_drifts(NULL)
{ 
  int i;
  std::ifstream is_inputs(path);
  if (!is_inputs.is_open()) throw std::logic_error("Failed to open input file in GaussMarkovWorld constructor");
  blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
  std::map<std::string,std::string> inputs_map;
  for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
  int number_of_economies = 0;
  if (inputs_map.count("Number of economies")) number_of_economies = std::atoi(inputs_map["Number of economies"].data());
  if (number_of_economies<1) throw std::logic_error("Misspecified number of economies");
  boost::shared_ptr<GaussianEconomy> economy;
  for (i=0;i<number_of_economies;i++) {
	std::string str("CSV file for economy ");
	str += (char)(48+i);
	if (inputs_map.count(str)) {
	  economy.reset(new GaussianEconomy(inputs_map[str].data()));
	  economies.push_back(economy); }
	else throw std::logic_error("Missing Gaussian economy specification"); }
  initial_exchange_rates.resize(number_of_economies-1);
  terminal_fwd_exchange_rate.resize(number_of_economies-1);
  int voldim = economies[0]->v->factors();
  termstructure_statevariables.resize(voldim);
  economy_start_index.resize(economies.size());
  boost::shared_ptr<DeterministicAssetVol> vol;
  for (i=1;i<number_of_economies;i++) {
	std::string str("Initial exchange rate ");
	str += (char)(48+i);
	if (inputs_map.count(str)) initial_exchange_rates(i-1) = std::atof(inputs_map[str].data());
	else throw std::logic_error("Missing initial exchange rate specification");
    str = "CSV file for FX volatility ";
	str += (char)(48+i);
	vol = CSVreadvolatility(inputs_map[str].data());
	vols.push_back(vol); }
  initialise();
}

