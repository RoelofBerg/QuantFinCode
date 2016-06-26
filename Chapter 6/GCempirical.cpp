/** \file  GCempirical.cpp
    \brief Test and demonstration program.
           Copyright 2011, 2015 by Erik Schlögl

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

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <boost/bind.hpp>
#include <boost/math/distributions/normal.hpp>
#include <random/normal.h>
#include "BlackScholesAsset.hpp"
#include "GramCharlierAsset.hpp"
#include "CSV2Array.hpp"

using namespace quantfin;
 
/** Test and demonstration for calibration of risk-neutral distribution based on Gram/Charlier expansion.

    Command-line arguments:
      -# input data CSV file name, e.g. USDEUR10_25delta080124.csv.
      -# highest non-normal moment.
         The default is 4 (i.e. non-zero excess kurtosis).
  */
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i;
  double Spot,ttm,r_d,r_f,Delta,volATM,volRR,volBF,Delta2,volRR2,volBF2;
  Delta2 = volRR2 = volBF2 = 0.0;
  try {
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
	if (argc<2) throw std::logic_error("Missing input data CSV file name");
	int highest_moment = 4;
	if (argc>2) highest_moment = std::atoi(argv[2]);
	std::ifstream is_inputs(argv[1]);
	blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
	std::cout << "Inputs: " << inputs_matrix << std::endl;
	std::map<std::string,std::string> inputs_map;
	for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
    if (inputs_map.count("Spot")) Spot = std::atof(inputs_map["Spot"].data());
    if (inputs_map.count("Delta")) Delta = std::atof(inputs_map["Delta"].data());
    if (inputs_map.count("ATM")) volATM = std::atof(inputs_map["ATM"].data());
    if (inputs_map.count("RR")) volRR = std::atof(inputs_map["RR"].data());
    if (inputs_map.count("BF")) volBF = std::atof(inputs_map["BF"].data());
    if (inputs_map.count("ttm")) ttm = std::atof(inputs_map["ttm"].data());
	else {
	  ttm = 0.0;
	  if (inputs_map.count("TTM")) {
	    if (inputs_map["TTM"]==std::string("1M")) ttm = 1.0/12.0;
		if (inputs_map["TTM"]==std::string("3M")) ttm = 0.25; }}
	if (ttm==0.0) throw std::logic_error("Unknown time to maturity");
    if (inputs_map.count("r_d")) {
	  r_d = std::atof(inputs_map["r_d"].data());
	  r_d = std::log(1+ttm*r_d); }
    if (inputs_map.count("r_f")) {
	  r_f = std::atof(inputs_map["r_f"].data());
	  r_f = std::log(1+ttm*r_f); }
	double sgmcdelta = volBF + 0.5*volRR + volATM;
	double sgmpdelta = sgmcdelta - volRR;
	double domestic_discount = exp(-r_d*ttm);
	double foreign_discount  = exp(-r_f*ttm);
	double forwardFX = Spot * foreign_discount/domestic_discount;
	double sqrtttm = std::sqrt(ttm);
	boost::math::normal normal;
	double d = boost::math::quantile(normal,Delta);
	double cdeltastrike = Spot * std::exp(-sgmcdelta*sqrtttm*d+(r_d-r_f+0.5*sgmcdelta*sgmcdelta)*ttm);
	d = boost::math::quantile(normal,1.0-Delta);
	double pdeltastrike = Spot * std::exp(-sgmpdelta*sqrtttm*d+(r_d-r_f+0.5*sgmpdelta*sgmpdelta)*ttm);
    if (inputs_map.count("Delta2")) Delta2 = std::atof(inputs_map["Delta2"].data());
    if (inputs_map.count("RR2")) volRR2 = std::atof(inputs_map["RR2"].data());
    if (inputs_map.count("BF2")) volBF2 = std::atof(inputs_map["BF2"].data());
	double sgmcdelta2 = volBF2 + 0.5*volRR2 + volATM;
	double sgmpdelta2 = sgmcdelta2 - volRR2;
	double cdelta2strike,pdelta2strike;
	int nstrikes = 3;
	if (Delta2>0.0) {
	  nstrikes = 5;
	  d = boost::math::quantile(normal,Delta2);
	  cdelta2strike = Spot * std::exp(-sgmcdelta2*sqrtttm*d+(r_d-r_f+0.5*sgmcdelta2*sgmcdelta2)*ttm);
	  d = boost::math::quantile(normal,1.0-Delta2);
	  pdelta2strike = Spot * std::exp(-sgmpdelta2*sqrtttm*d+(r_d-r_f+0.5*sgmpdelta2*sgmpdelta2)*ttm); }

	Array<double,1> coeff(highest_moment+1);
	coeff = 0.0;
	coeff(0) = 1.0;
    GramCharlier gc(coeff);
	GramCharlierAsset gcasset(gc,volATM*sqrtttm,Spot,ttm);

	cout << "Calibrate:" << endl;
	Array<double,1> strikes(nstrikes);
	Array<double,1> vols(nstrikes);
	strikes(blitz::Range(0,2)) = pdeltastrike,forwardFX,cdeltastrike;
	vols(blitz::Range(0,2))    = sgmpdelta,volATM,sgmcdelta;
	if (nstrikes==5) {
	  strikes(3) = pdelta2strike;
	  strikes(4) = cdelta2strike;
	  vols(3)    = sgmpdelta2;
	  vols(4)    = sgmcdelta2; }
	cout << "Strikes: " << strikes << "\nVolatilities: " << vols << endl;
    double match = gcasset.calibrate(&strikes,&vols,domestic_discount,foreign_discount,highest_moment);
	cout << "Calibration result:\nObjective function: " << match << endl;
	cout << "Sigma:, " << gcasset.standard_deviation() << "\nSkewness:, " << gcasset.skewness() << ", Excess kurtosis:, " << gcasset.excess_kurtosis() << endl;
	cout << "GC coefficients: " << gc.coefficients()  << endl;
	cout << "Strike,Implied Volatility,Black/Scholes Price,Fitted Gram/Charlier Price" << endl;
	for (i=0;i<nstrikes;i++) {
	  cout << strikes(i) << ',' << vols(i) << ',';
	  ConstVol v(vols(i));
	  BlackScholesAsset asset(&v,Spot);
	  asset.dividend_yield(r_f);
	  cout << asset.option(ttm,strikes(i),r_d) << ',' << gcasset.call(strikes(i),domestic_discount,foreign_discount) << endl; }
	cout << "Generating data for a smile plot:\nStrike,Price,Implied Volatility" << endl;
	ConstVol cv(volATM);
	BlackScholesAsset stock(&cv,Spot,r_f);
	double K;
	double startK = 0.8*strikes(0);
	double endK   = 1.2 * strikes(2);
	double dK     = (endK-startK)/100.0;
	for (K=startK;K<=endK;K+=dK) {
	  double price = gcasset.call(K,domestic_discount,foreign_discount);
	  cout << K << ',' << price << ',' << stock.implied_volatility(price,ttm,K,r_d) << endl; }
      
	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
