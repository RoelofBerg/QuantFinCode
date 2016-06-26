/** \file  MBinaryGaussianTest.cpp
    \brief Test and demonstration program.
           Copyright 2013 by Erik Schlögl

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
#include "CSV2Array.hpp"
#include "MBinaryPayoffs.hpp"
#include "GaussianEconomy.hpp"
#include "ExponentialVol.hpp"
#include "MCEngine.hpp"
#include "QFRandom.hpp"
#include "MCGeneric.hpp"
#include "MExotics.hpp"

using namespace quantfin;
 
/** Test and demonstration for GBM/HJM options by MBinary.

    Command-line arguments:
	  -# option parameter CSV file name.
      -# GaussMarkovWorld CSV file name.
      -# Minimum number of MC pricing paths.
         The default is 100.
      -# Maximum number of MC pricing paths.
         The default is 100.
	  -# Numeraire index.
	     The default is 0 (domestic discrete rolling spot measure).

	Usage e.g. MBinaryGaussianTest.exe option.csv worldbsf.csv 100 1000000000 0
  */
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::cerr;

  int i;
  double mat,strike,diff;
  bool passed = true;
  double eps = 1e-9;
  int numeraire_index = 0;
  try {
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
	if (argc<2) throw std::logic_error("Missing option parameter CSV file name");
	if (argc<3) throw std::logic_error("Missing GaussMarkovWorld CSV file name");
    size_t minpaths = 100;
    if (argc>3) minpaths = atoi(argv[3]);
    size_t maxpaths = 100;
    if (argc>4) maxpaths = atoi(argv[4]);
    if (argc>5) numeraire_index = atoi(argv[5]);
	// Read data from files and create multicurrency term structure model
    std::ifstream is_inputs(argv[1]);
    if (!is_inputs.is_open()) throw std::logic_error("Failed to open option parameter CSV file");
    blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
    std::map<std::string,std::string> inputs_map;
    for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
    std::string str("Option expiry");
    if (inputs_map.count(str)) mat = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing option expiry"); 
	str = "Strike";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing strike"); 
	double strikeFX    = strike;
	double strikeQS    = strike;
	double strikeFBS   = strike;
	double strikeZCB   = strike;
	double strikeFZCB  = strike;
	double lvl         = strike;
	double strikeIS    = strike;
	double strikeIFX   = strike;
	double strikeIQS   = strike;
	double strikeIFBS  = strike;
	double strikeIZCB  = strike;
	double strikeIFZCB = strike;
	double Ilvl        = strike;
	GaussMarkovWorld world(argv[2]);
	std::cerr << "Getting economy..." << endl;
	const std::vector<boost::shared_ptr<GaussianEconomy> >& ecv = world.get_economies();
	const GaussianEconomy& ec = *(ecv[0]);
	std::cerr << "Number of assets: " << ec.underlying.size() << endl;
	std::cerr << "Getting asset..." << endl;
	boost::shared_ptr<BlackScholesAsset> Sp(ec.underlying[0]);
	std::cerr << "...got pointer..." << endl;
	const BlackScholesAsset& S = *Sp;
	std::cerr << "...got asset." << endl;
	strike *= S.initial_value() / (*world.get_economies()[0]->initialTS)(mat);
	std::cerr << "Strike: " << strike << endl;
	double CFcall = world.get_economies()[0]->hjm->option(S,mat,strike);
	cout << "Closed form value: " << CFcall << endl;
	// Underlying is first asset in domestic economy
	Array<double,1> T(2);
	T = 0.0, mat;
	world.set_timeline(T);
	world.set_reporting(0,0);
	// Price option using MBinary
	/*
	boost::shared_ptr<AssetBinary> A_payoff(new AssetBinary(S,strike,0.0,mat,*world.get_economies()[0]->initialTS));
	boost::shared_ptr<BondBinary>  B_payoff(new BondBinary(S,strike,0.0,mat,*world.get_economies()[0]->initialTS));
	*/
	boost::shared_ptr<AssetBinary> A_payoff(new AssetBinary(world,strike,0.0,mat,0));
	boost::shared_ptr<BondBinary>  B_payoff(new BondBinary(world,strike,0.0,mat,0));
	MBinary A(world,*A_payoff);
	MBinary B(world,*B_payoff);
	double Avalue = A.price();
	double Bvalue = B.price();
	double MBCFcall = Avalue - strike*Bvalue;
	cout << "Closed form value using MBinary: " << MBCFcall << "   Diff: " << (diff = CFcall-MBCFcall) << endl;
	if (std::abs(diff)>eps) passed = false;
	double imat = mat/2.0;
	strikeIS *= S.initial_value() / (*world.get_economies()[0]->initialTS)(imat);
	std::cerr << "Strike: " << strikeIS << endl;
	strikeFX  *= world.get_forward_exchange_rate(1,mat);
	strikeIFX *= world.get_forward_exchange_rate(1,imat);
	//strikeFX = 0.0;
	double CFcallFX = world.get_economies()[0]->hjm->FXoption(world.get_initial_exchange_rates()(0),
		                                                      mat,
															  strikeFX,
															  *(world.get_economies()[1]->hjm),
															  *(world.get_FXvolatilities()[0]));
	cout << "Closed form value for FX call: " << CFcallFX << endl;
	double CFcallIFX = world.get_economies()[0]->hjm->FXoption(world.get_initial_exchange_rates()(0),
		                                                       imat,
															   strikeIFX,
															   *(world.get_economies()[1]->hjm),
															   *(world.get_FXvolatilities()[0]));
	cout << "Closed form value for intermediate maturity FX call: " << CFcallIFX << endl;
	// Underlying is foreign currency 1
	int reportable_FX_index = world.set_reporting(1,-2,mat);
	// Price option using MBinary
	boost::shared_ptr<AssetBinary> A_payoffFX(new AssetBinary(world,strikeFX,0.0,mat,1));
	boost::shared_ptr<BondBinary>  B_payoffFX(new BondBinary(world,strikeFX,0.0,mat,1));
	MBinary AFX(world,*A_payoffFX);
	MBinary BFX(world,*B_payoffFX);
	double AFXvalue = AFX.price();
	double BFXvalue = BFX.price();
	double MBCFcallFX = AFXvalue - strikeFX*BFXvalue;
	cout << "FX option closed form value using MBinary: " << MBCFcallFX << "   Diff: " << (diff = CFcallFX-MBCFcallFX) << endl;
	if (std::abs(diff)>eps) passed = false;
	int iZCBindex   = world.set_reporting(0,-1,mat);
	int iFZCBindex  = world.set_reporting(1,-1,mat);
	boost::shared_ptr<MBinaryPayoff>  A_IFXpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,3,1));
	boost::shared_ptr<MBinaryPayoff>  B_IFXpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,3,1));
    A_IFXpayoff->timeline = 0.0, imat;
    A_IFXpayoff->index    = reportable_FX_index, iZCBindex, iFZCBindex,
		                    1, 1, 1;      
    A_IFXpayoff->alpha    = 1.0, 1.0, -1.0;
    A_IFXpayoff->S        = 1.0;
    A_IFXpayoff->A        = 1.0, 1.0, -1.0;
    A_IFXpayoff->a        = strikeIFX;
    B_IFXpayoff->timeline = 0.0, imat;
    B_IFXpayoff->index    = reportable_FX_index, iZCBindex, iFZCBindex,
		                    1, 1, 1;      
    B_IFXpayoff->alpha    = 0.0;
    B_IFXpayoff->S        = 1.0;
    B_IFXpayoff->A        = 1.0, 1.0, -1.0;
    B_IFXpayoff->a        = strikeIFX;
	MBinary AIFX(world,*A_IFXpayoff);
	MBinary BIFX(world,*B_IFXpayoff);
	double AIFXvalue = AIFX.price();
	double BIFXvalue = BIFX.price();
	double MBCFcallIFX = AIFXvalue - strikeIFX*BIFXvalue;
	cout << "FX intermediate maturity option closed form value using MBinary: " << MBCFcallIFX << "   Diff: " << (diff = CFcallIFX-MBCFcallIFX) << endl;
	if (std::abs(diff)>eps) passed = false;
	// First asset in foreign economy
	int reportable_asset_index = world.set_reporting(1,0);
    //Array<int,1> reportable_asset_index(2);
	//reportable_asset_index = 1, 2; // Underlying is first asset in foreign economy, multiplied by exchange rate
	const GaussianEconomy& ecF = *(ecv[1]);
	std::cerr << "Number of assets: " << ecF.underlying.size() << endl;
	std::cerr << "Getting asset..." << endl;
	boost::shared_ptr<BlackScholesAsset> SpF(ecF.underlying[0]);
	std::cerr << "...got pointer..." << endl;
	const BlackScholesAsset& SF = *SpF;
	strikeFBS  *= SF.initial_value() / (*world.get_economies()[1]->initialTS)(mat);
	strikeIFBS *= SF.initial_value() / (*world.get_economies()[1]->initialTS)(imat);
	//strikeFBS = 0.0;
	double CFcallFBS = world.get_economies()[1]->hjm->option(SF,mat,strikeFBS);
	cout << "Closed form value for foreign asset call option: " << CFcallFBS << endl;
	double fxspot = (world.get_initial_exchange_rates())(0);
	CFcallFBS *= fxspot;
	cout << "Closed form value converted to domestic currency: " << CFcallFBS << endl;
	double CFcallIFBS = world.get_economies()[1]->hjm->option(SF,imat,strikeIFBS);
	cout << "Closed form value for intermediate maturity foreign asset call option: " << CFcallIFBS << endl;
	CFcallIFBS *= fxspot;
	cout << "Closed form value converted to domestic currency: " << CFcallIFBS << endl;
	//strikeFBS *= fxspot;
	//strikeFBS *= world.get_forward_exchange_rate(1,mat);
	//Array<double,1> alpha(2);
	//alpha = 1.0, 0.0;
	//AssetProductBinary A_payoffFBS(world,strikeFBS,0.0,mat,reportable_asset_index);
	//AssetProductBinary B_payoffFBS(world,strikeFBS,0.0,mat,reportable_asset_index,alpha);
	boost::shared_ptr<ForeignOption> A_payoffFBS(new ForeignOption(world,strikeFBS,0.0,mat,reportable_asset_index,reportable_FX_index,true));
	boost::shared_ptr<ForeignOption> B_payoffFBS(new ForeignOption(world,strikeFBS,0.0,mat,reportable_asset_index,reportable_FX_index,false));
	MBinary AFBS(world,*A_payoffFBS);
	MBinary BFBS(world,*B_payoffFBS);
	double AFBSvalue = AFBS.price();
	double BFBSvalue = BFBS.price();
	double MBCFcallFBS = AFBSvalue - strikeFBS*BFBSvalue;
	cout << "Foreign asset option closed form value using MBinary: " << MBCFcallFBS << "   Diff: " << (diff = CFcallFBS-MBCFcallFBS) << endl;
	if (std::abs(diff)>eps) passed = false;
	double IFBSdiv = SF.dividend_discount(imat,mat);
	boost::shared_ptr<MBinaryPayoff>  A_IFBSpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,4,1));
	boost::shared_ptr<MBinaryPayoff>  B_IFBSpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,4,1));
    A_IFBSpayoff->timeline = 0.0, imat;
    A_IFBSpayoff->index    = reportable_asset_index, reportable_FX_index, iZCBindex, iFZCBindex,
		                    1, 1, 1, 1;      
    A_IFBSpayoff->alpha    = 1.0, 1.0, 1.0, 0.0;
    A_IFBSpayoff->S        = 1.0;
    A_IFBSpayoff->A        = 1.0, 0.0, 0.0, 1.0;
    A_IFBSpayoff->a        = IFBSdiv*strikeIFBS;
    B_IFBSpayoff->timeline = 0.0, imat;
    B_IFBSpayoff->index    = reportable_asset_index, reportable_FX_index, iZCBindex, iFZCBindex,
		                    1, 1, 1, 1;       
    B_IFBSpayoff->alpha    = 0.0, 1.0, 1.0, -1.0;
    B_IFBSpayoff->S        = 1.0;
    B_IFBSpayoff->A        = 1.0, 0.0, 0.0, 1.0;
    B_IFBSpayoff->a        = IFBSdiv*strikeIFBS;
	MBinary AIFBS(world,*A_IFBSpayoff);
	MBinary BIFBS(world,*B_IFBSpayoff);
	double AIFBSvalue = AIFBS.price();
	double BIFBSvalue = BIFBS.price();
	double MBCFcallIFBS = (AIFBSvalue - IFBSdiv*strikeIFBS*BIFBSvalue)/IFBSdiv;
	cout << "Foreign asset intermediate maturity option closed form value using MBinary: " << MBCFcallIFBS << "   Diff: " << (diff = CFcallIFBS-MBCFcallIFBS) << endl;
	if (std::abs(diff)>eps) passed = false;
	// Call option on domestic zero coupon bond
	double bondmat = 2.5*mat;
	strikeZCB  *= (*world.get_economies()[0]->initialTS)(bondmat) / (*world.get_economies()[0]->initialTS)(mat);
	strikeIZCB *= (*world.get_economies()[0]->initialTS)(bondmat) / (*world.get_economies()[0]->initialTS)(imat);
	//strikeZCB = 0.0;
    double CFZCBcall = world.get_economies()[0]->hjm->ZCBoption(mat,bondmat,strikeZCB);
    cout << "Closed form call on domestic zero coupon bond: " << CFZCBcall << endl;
    double CFIZCBcall = world.get_economies()[0]->hjm->ZCBoption(imat,bondmat,strikeIZCB);
    cout << "Closed form intermediate maturity call on domestic zero coupon bond: " << CFIZCBcall << endl;
	int zcbidx = world.set_reporting(0,-1,bondmat);
	// Price option using MBinary
	boost::shared_ptr<AssetBinary> A_ZCBpayoff(new AssetBinary(world,strikeZCB,0.0,mat,zcbidx));
	boost::shared_ptr<BondBinary>  B_ZCBpayoff(new BondBinary(world,strikeZCB,0.0,mat,zcbidx));
	MBinary AZCB(world,*A_ZCBpayoff);
	MBinary BZCB(world,*B_ZCBpayoff);
	double AZCBvalue = AZCB.price();
	double BZCBvalue = BZCB.price();
	double MBCFZCBcall = AZCBvalue - strikeZCB*BZCBvalue;
	cout << "Closed form value using MBinary: " << MBCFZCBcall << "   Diff: " << (diff = CFZCBcall-MBCFZCBcall) << endl;
	if (std::abs(diff)>eps) passed = false;
	boost::shared_ptr<AssetBinary> A_IZCBpayoff(new AssetBinary(world,strikeIZCB,0.0,imat,zcbidx));
	boost::shared_ptr<BondBinary>  B_IZCBpayoff(new BondBinary(world,strikeIZCB,0.0,imat,zcbidx));
	MBinary AIZCB(world,*A_IZCBpayoff);
	MBinary BIZCB(world,*B_IZCBpayoff);
	double AIZCBvalue = AIZCB.price();
	double BIZCBvalue = BIZCB.price();
	double MBCFIZCBcall = AIZCBvalue - strikeIZCB*BIZCBvalue;
	cout << "Closed form intermediate maturity value using MBinary: " << MBCFIZCBcall << "   Diff: " << (diff = CFIZCBcall-MBCFIZCBcall) << endl;
	if (std::abs(diff)>eps) passed = false;
	// Call option on foreign zero coupon bond
	strikeFZCB  *= (*world.get_economies()[1]->initialTS)(bondmat) / (*world.get_economies()[1]->initialTS)(mat);
	strikeIFZCB *= (*world.get_economies()[1]->initialTS)(bondmat) / (*world.get_economies()[1]->initialTS)(imat);
	//strikeFZCB = 0.0;
    double CFFZCBcall = world.get_economies()[1]->hjm->ZCBoption(mat,bondmat,strikeFZCB);
    cout << "Closed form call on foreign zero coupon bond: " << CFFZCBcall << endl;
	CFFZCBcall *= fxspot;
    cout << "Price in domestic currency: " << CFFZCBcall << endl;
    double CFIFZCBcall = world.get_economies()[1]->hjm->ZCBoption(imat,bondmat,strikeIFZCB);
    cout << "Closed form intermediate maturity call on foreign zero coupon bond: " << CFIFZCBcall << endl;
	CFIFZCBcall *= fxspot;
    cout << "Price in domestic currency: " << CFIFZCBcall << endl;
	int fzcbidx = world.set_reporting(1,-1,bondmat);
	// Price option using MBinary
	boost::shared_ptr<ForeignOption> A_FZCBpayoff(new ForeignOption(world,strikeFZCB,0.0,mat,fzcbidx,reportable_FX_index,true));
	boost::shared_ptr<ForeignOption> B_FZCBpayoff(new ForeignOption(world,strikeFZCB,0.0,mat,fzcbidx,reportable_FX_index,false));
	MBinary AFZCB(world,*A_FZCBpayoff);
	MBinary BFZCB(world,*B_FZCBpayoff);
	double AFZCBvalue = AFZCB.price();
	double BFZCBvalue = BFZCB.price();
	double MBCFcallFZCB = AFZCBvalue - strikeFZCB*BFZCBvalue;
	cout << "Foreign zero coupon bond option closed form value using MBinary: " << MBCFcallFZCB << "   Diff: " << (diff = CFFZCBcall-MBCFcallFZCB) << endl;
	if (std::abs(diff)>eps) passed = false;
	boost::shared_ptr<MBinaryPayoff>  A_IFZCBpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,4,1));
	boost::shared_ptr<MBinaryPayoff>  B_IFZCBpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,4,1));
    A_IFZCBpayoff->timeline = 0.0, imat;
    A_IFZCBpayoff->index    = fzcbidx, reportable_FX_index, iZCBindex, iFZCBindex,
		                      1, 1, 1, 1;      
    A_IFZCBpayoff->alpha    = 1.0, 1.0, 1.0, -1.0;
    A_IFZCBpayoff->S        = 1.0;
    A_IFZCBpayoff->A        = 1.0, 0.0, 0.0, 0.0;
    A_IFZCBpayoff->a        = strikeIFZCB;
    B_IFZCBpayoff->timeline = 0.0, imat;
    B_IFZCBpayoff->index    = fzcbidx, reportable_FX_index, iZCBindex, iFZCBindex,
		                      1, 1, 1, 1;       
    B_IFZCBpayoff->alpha    = 0.0, 1.0, 1.0, -1.0;
    B_IFZCBpayoff->S        = 1.0;
    B_IFZCBpayoff->A        = 1.0, 0.0, 0.0, 0.0;
    B_IFZCBpayoff->a        = strikeIFZCB;
	MBinary AIFZCB(world,*A_IFZCBpayoff);
	MBinary BIFZCB(world,*B_IFZCBpayoff);
	double AIFZCBvalue = AIFZCB.price();
	double BIFZCBvalue = BIFZCB.price();
	double MBCFcallIFZCB = AIFZCBvalue - strikeIFZCB*BIFZCBvalue;
	cout << "Foreign zero coupon bond intermediate maturity option closed form value using MBinary: " << MBCFcallIFZCB << "   Diff: " << (diff = CFIFZCBcall-MBCFcallIFZCB) << endl;
	if (std::abs(diff)>eps) passed = false;
	// Quanto caplet
	double delta = 0.75;
	lvl *= (*world.get_economies()[1]->initialTS).simple_rate(mat,delta);
	Ilvl *= (*world.get_economies()[1]->initialTS).simple_rate(imat,delta);
	double CFquanto = world.get_economies()[0]->hjm->QuantoCaplet(fxspot,mat,delta,lvl,*world.get_economies()[1]->hjm,*((world.get_FXvolatilities())[0]));
	cout << "Closed form value for quanto caplet: " << CFquanto << endl;
	double CFquanto_old = world.get_economies()[0]->hjm->QuantoCaplet_old(fxspot,mat,delta,lvl,*world.get_economies()[1]->hjm,*((world.get_FXvolatilities())[0]));
	cout << "Closed form value for quanto caplet: " << CFquanto_old << endl;
	double CFIquanto = world.get_economies()[0]->hjm->QuantoCaplet(fxspot,imat,delta,Ilvl,*world.get_economies()[1]->hjm,*((world.get_FXvolatilities())[0]));
	cout << "Closed form value for intermediate maturity quanto caplet: " << CFIquanto << endl;
	// Price option using MBinary
	int Qidx  = world.set_reporting(0,-1,mat+delta);
	int QFidx = world.set_reporting(1,-1,mat+delta);
	double strikeQ = 1.0+lvl*delta;
	boost::shared_ptr<MBinaryPayoff>  A_Qpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,2,1));
	boost::shared_ptr<MBinaryPayoff>  B_Qpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,2,1));
    A_Qpayoff->timeline = 0.0, mat;
    A_Qpayoff->index    = Qidx, QFidx,
		                  1, 1;      
    A_Qpayoff->alpha    = 1.0, -1.0;
    A_Qpayoff->S        = 1.0;
    A_Qpayoff->A        = 0.0, -1.0;
    A_Qpayoff->a        = strikeQ;
    B_Qpayoff->timeline = 0.0, mat;
    B_Qpayoff->index    = Qidx, QFidx,
		                  1, 1;      
    B_Qpayoff->alpha    = 1.0, 0.0;
    B_Qpayoff->S        = 1.0;
    B_Qpayoff->A        = 0.0, -1.0;
    B_Qpayoff->a        = strikeQ;
	MBinary Aquanto(world,*A_Qpayoff);
	MBinary Bquanto(world,*B_Qpayoff);
	double AQvalue = Aquanto.price();
	double BQvalue = Bquanto.price();
	double MBCFquanto = AQvalue - strikeQ*BQvalue;
	MBCFquanto *= fxspot;
	cout << "Closed form value using MBinary: " << MBCFquanto << "   Diff: " << (diff = CFquanto-MBCFquanto) << endl;
	diff /= CFquanto;
	if (std::abs(diff)>eps) passed = false;
	int IQidx  = world.set_reporting(0,-1,imat+delta);
	int IQFidx = world.set_reporting(1,-1,imat+delta);
	double strikeIQ = 1.0+Ilvl*delta;
	boost::shared_ptr<MBinaryPayoff>  A_IQpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,2,1));
	boost::shared_ptr<MBinaryPayoff>  B_IQpayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,2,1));
    A_IQpayoff->timeline = 0.0, imat;
    A_IQpayoff->index    = IQidx, IQFidx,
		                  1, 1;      
    A_IQpayoff->alpha    = 1.0, -1.0;
    A_IQpayoff->S        = 1.0;
    A_IQpayoff->A        = 0.0, -1.0;
    A_IQpayoff->a        = strikeIQ;
    B_IQpayoff->timeline = 0.0, imat;
    B_IQpayoff->index    = IQidx, IQFidx,
		                  1, 1;      
    B_IQpayoff->alpha    = 1.0, 0.0;
    B_IQpayoff->S        = 1.0;
    B_IQpayoff->A        = 0.0, -1.0;
    B_IQpayoff->a        = strikeIQ;
	MBinary AIquanto(world,*A_IQpayoff);
	MBinary BIquanto(world,*B_IQpayoff);
	double AIQvalue = AIquanto.price();
	double BIQvalue = BIquanto.price();
	double MBCFIquanto = AIQvalue - strikeIQ*BIQvalue;
	MBCFIquanto *= fxspot;
	cout << "Closed form intermediate maturity value using MBinary: " << MBCFIquanto << "   Diff: " << (diff = CFIquanto-MBCFIquanto) << endl;
	diff /= CFIquanto;
	if (std::abs(diff)>eps) passed = false;
	double CFcallIS = world.get_economies()[0]->hjm->option(S,imat,strikeIS);
	cout << "Closed form value of intermediate maturity option: " << CFcallIS << endl;
	// Price intermediate maturity option using MBinary
	// MC works for this one. 
	// int iZCBindex  = world.set_reporting(0,-1,mat);
	double IBSdiv = S.dividend_discount(imat,mat);
 	Array<int,1> product_index(2);
	product_index = 0, iZCBindex;
	boost::shared_ptr<AssetProductBinary> A_Ipayoff(new AssetProductBinary(world,strikeIS*IBSdiv,0.0,imat,product_index));
	boost::shared_ptr<BondProductBinary>  B_Ipayoff(new BondProductBinary(world,strikeIS*IBSdiv,0.0,imat,product_index));
	/*
	boost::shared_ptr<MBinaryPayoff>  A_Ipayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,3,1));
	boost::shared_ptr<MBinaryPayoff>  B_Ipayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,2,3,1));
 	Array<int,1> product_index(3);
	int iZCBindex  = world.set_reporting(0,-1,mat);
 	int iZCBindex2 = world.set_reporting(0,-1,imat);
	product_index = 0, iZCBindex, iZCBindex2;
    A_Ipayoff->timeline = 0.0, imat;
    A_Ipayoff->index    = 0, iZCBindex, iZCBindex2,
		                  1, 1, 1;      
    A_Ipayoff->alpha    = 1.0, 1.0, -1.0;
    A_Ipayoff->S        = 1.0;
    A_Ipayoff->A        = 1.0, 1.0, -1.0;
    A_Ipayoff->a        = strikeIS;
    B_Ipayoff->timeline = 0.0, imat;
    B_Ipayoff->index    = 0, iZCBindex, iZCBindex2,
		                  1, 1, 1;      
    B_Ipayoff->alpha    = 0.0;
    B_Ipayoff->S        = 1.0;
    B_Ipayoff->A        = 1.0, 1.0, -1.0;
    B_Ipayoff->a        = strikeIS;
	*/
	/*
	boost::shared_ptr<MBinaryPayoff>  A_Ipayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,3,2,1));
	boost::shared_ptr<MBinaryPayoff>  B_Ipayoff(new MBinaryPayoff(*world.get_economies()[0]->initialTS,3,2,1));
 	Array<int,1> product_index(2);
	int iZCBindex  = world.set_reporting(0,-1,mat);
	product_index = 0, iZCBindex;
    A_Ipayoff->timeline = 0.0, imat, mat;
    A_Ipayoff->index    = 0, iZCBindex,
		                  1, 1;      
    A_Ipayoff->alpha    = 1.0, 0.0;
    A_Ipayoff->S        = 1.0;
    A_Ipayoff->A        = 1.0, 1.0;
    A_Ipayoff->a        = strikeIS;
    B_Ipayoff->timeline = 0.0, imat, mat;
    B_Ipayoff->index    = 0, iZCBindex;
		                  1, 1;      
    B_Ipayoff->alpha    = -1.0;
    B_Ipayoff->S        = 1.0;
    B_Ipayoff->A        = 1.0, 1.0;
    B_Ipayoff->a        = strikeIS;
	*/
	/*
	boost::shared_ptr<AssetBinary> A_ITpayoff(new AssetBinary(world,strikeIS,0.0,imat,0));
	boost::shared_ptr<BondBinary>  B_ITpayoff(new BondBinary(world,strikeIS,0.0,imat,0));
	MBinary AIS_test(world,*A_ITpayoff);
	MBinary BIS_test(world,*B_ITpayoff);
	double AIS_testvalue = AIS_test.price();
	double BIS_testvalue = BIS_test.price();
	double MBCFcallIS_test = AIS_testvalue - strikeIS*BIS_testvalue;
	cout << "Closed form value using MBinary: " << MBCFcallIS_test << "   Diff: " << (diff = CFcallIS-MBCFcallIS_test) << endl;
	if (std::abs(diff)>eps) passed = false;
	*/
	MBinary AIS(world,*A_Ipayoff);
	MBinary BIS(world,*B_Ipayoff);
	double AISvalue = AIS.price();
	double BISvalue = BIS.price();
	double MBCFcallIS = (AISvalue - IBSdiv*strikeIS*BISvalue)/IBSdiv;
	cout << "Closed form value using MBinary: " << MBCFcallIS << "   Diff: " << (diff = CFcallIS-MBCFcallIS) << endl;
	if (std::abs(diff)>eps) passed = false;
	// permutate
	product_index = iZCBindex, 0;
	boost::shared_ptr<AssetProductBinary> A_IPpayoff(new AssetProductBinary(world,IBSdiv*strikeIS,0.0,imat,product_index));
	boost::shared_ptr<BondProductBinary>  B_IPpayoff(new BondProductBinary(world,IBSdiv*strikeIS,0.0,imat,product_index));
	MBinary AIPS(world,*A_IPpayoff);
	MBinary BIPS(world,*B_IPpayoff);
	AISvalue = AIPS.price();
	BISvalue = BIPS.price();
	double MBCFcallIPS = (AISvalue - IBSdiv*strikeIS*BISvalue)/IBSdiv;
	cout << "Closed form value using MBinary: " << MBCFcallIPS << "   Diff: " << (diff = CFcallIS-MBCFcallIPS) << endl;
	if (std::abs(diff)>eps) passed = false;
	if (passed) cout << "************ PASSED ************" << endl;
	else        cout << "************ FAILED ************" << endl;
	if (passed) cerr << "************ PASSED ************" << endl;
	else        cerr << "************ FAILED ************" << endl;

	// Price asset quanto option using MBinary
	int QSidx = reportable_asset_index;
	strikeQS *= SF.initial_value() / (*world.get_economies()[1]->initialTS)(mat);
	boost::shared_ptr<AssetBinary> A_QSpayoff(new AssetBinary(world,strikeQS,0.0,mat,QSidx));
	boost::shared_ptr<BondBinary>  B_QSpayoff(new BondBinary(world,strikeQS,0.0,mat,QSidx));
	MBinary ASquanto(world,*A_QSpayoff);
	MBinary BSquanto(world,*B_QSpayoff);
	double AQSvalue = ASquanto.price();
	double BQSvalue = BSquanto.price();
	double MBCFSquanto = AQSvalue - strikeQS*BQSvalue;
	MBCFSquanto *= fxspot;
	cout << "Closed form value for asset quanto option using MBinary: " << MBCFSquanto << endl;
	// Price quanto zero coupon bond option using MBinary
	boost::shared_ptr<AssetBinary> A_QZCBpayoff(new AssetBinary(world,strikeFZCB,0.0,mat,fzcbidx));
	boost::shared_ptr<BondBinary>  B_QZCBpayoff(new BondBinary(world,strikeFZCB,0.0,mat,fzcbidx));
	MBinary AZCBquanto(world,*A_QZCBpayoff);
	MBinary BZCBquanto(world,*B_QZCBpayoff);
	double AQZCBvalue = AZCBquanto.price();
	double BQZCBvalue = BZCBquanto.price();
	double MBCFZCBquanto = AQZCBvalue - strikeFZCB*BQZCBvalue;
	MBCFZCBquanto *= fxspot;
	cout << "Closed form value for quanto zero coupon bond option using MBinary: " << MBCFZCBquanto << endl;
	// Call option on the domestic geometric average
	int asian_n = 11;
	Array<double,1> asianT(asian_n);
	double dt = mat/(asian_n-1.0);
	firstIndex idx;
	asianT = idx*dt;
	double asianK = std::sqrt(strike*S.initial_value());
	exotics::DiscreteGeometricMeanFixedStrike asian_option(world,asianT,asianK,0,0,1.0);
	double MBCFasian = asian_option.price();
	cout << "Closed form value for call option on the domestic geometric average using MBinary: " << MBCFasian << endl;
	// Call option on the foreign geometric average
	double asianFK = std::sqrt(strikeFBS*SF.initial_value());
	exotics::DiscreteGeometricMeanFixedStrike asianF_option(world,asianT,asianFK,reportable_asset_index,0,1.0);
	double MBCFasianF = asianF_option.price();
	cout << "Closed form value for call option on the foreign geometric average using MBinary: " << MBCFasianF << endl;

	// Monte Carlo
	int number_of_options = 17;
	i = 0;
	Array<double,1> CFvalues(number_of_options);
	Array<double,1> MBCFvalues(number_of_options);
	Array<std::string,1> labels(number_of_options);
    ranlib::NormalUnit<double> normalRNG;
	MCPayoffList mclist;
	mclist.push_back(asian_option.get_payoff());
	CFvalues(i) = MBCFasian;
	MBCFvalues(i) = MBCFasian;
	labels(i) = "Domestic DiscreteGeometricMeanFixedStrike";
	i++;
	boost::shared_ptr<BondBinary> ZCBpayoff(new BondBinary(S,0.0,0.0,mat,*world.get_economies()[0]->initialTS));
	mclist.push_back(ZCBpayoff);
	MBCFvalues(i) = CFvalues(i) = (*world.get_economies()[0]->initialTS)(mat);
	labels(i) = "ZCB";
	i++;
	boost::shared_ptr<MCPayoffList> calloption(new MCPayoffList);
	calloption->push_back(A_payoff);
	calloption->push_back(B_payoff,-strike);
	mclist.push_back(calloption);
	CFvalues(i) = CFcall;
	MBCFvalues(i) = MBCFcall;
	labels(i) = "Call on domestic asset";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIS(new MCPayoffList);
	calloptionIS->push_back(A_Ipayoff);
	calloptionIS->push_back(B_Ipayoff,-strikeIS*IBSdiv);
	mclist.push_back(calloptionIS,1.0/IBSdiv);
	CFvalues(i) = CFcallIS;
	MBCFvalues(i) = MBCFcallIS;
	labels(i) = "Call on domestic asset - intermediate maturity";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionFX(new MCPayoffList);
	calloptionFX->push_back(A_payoffFX);
	calloptionFX->push_back(B_payoffFX,-strikeFX);
	mclist.push_back(calloptionFX);
	CFvalues(i) = CFcallFX;
	MBCFvalues(i) = MBCFcallFX;
	labels(i) = "FX call";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIFX(new MCPayoffList);
	calloptionIFX->push_back(A_IFXpayoff);
	calloptionIFX->push_back(B_IFXpayoff,-strikeIFX);
	mclist.push_back(calloptionIFX);
	CFvalues(i) = CFcallIFX;
	MBCFvalues(i) = MBCFcallIFX;
	labels(i) = "FX intermediate maturity call";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionFBS(new MCPayoffList);
	calloptionFBS->push_back(A_payoffFBS);
	calloptionFBS->push_back(B_payoffFBS,-strikeFBS);
	mclist.push_back(calloptionFBS);
	CFvalues(i) = CFcallFBS;
	MBCFvalues(i) = MBCFcallFBS;
	labels(i) = "Call on foreign asset";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIFBS(new MCPayoffList);
	calloptionIFBS->push_back(A_IFBSpayoff);
	calloptionIFBS->push_back(B_IFBSpayoff,-strikeIFBS*IFBSdiv);
	mclist.push_back(calloptionIFBS,1.0/IFBSdiv);
	CFvalues(i) = CFcallIFBS;
	MBCFvalues(i) = MBCFcallIFBS;
	labels(i) = "Foreign asset intermediate maturity call";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionZCB(new MCPayoffList);
	calloptionZCB->push_back(A_ZCBpayoff);
	calloptionZCB->push_back(B_ZCBpayoff,-strikeZCB);
	mclist.push_back(calloptionZCB);
	CFvalues(i) = CFZCBcall;
	MBCFvalues(i) = MBCFZCBcall;
	labels(i) = "Call on domestic ZCB";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIZCB(new MCPayoffList);
	calloptionIZCB->push_back(A_IZCBpayoff);
	calloptionIZCB->push_back(B_IZCBpayoff,-strikeIZCB);
	mclist.push_back(calloptionIZCB);
	CFvalues(i) = CFIZCBcall;
	MBCFvalues(i) = MBCFIZCBcall;
	labels(i) = "Intermediate maturity call on domestic ZCB";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionFZCB(new MCPayoffList);
	calloptionFZCB->push_back(A_FZCBpayoff);
	calloptionFZCB->push_back(B_FZCBpayoff,-strikeFZCB);
	mclist.push_back(calloptionFZCB);
	CFvalues(i) = CFFZCBcall;
	MBCFvalues(i) = MBCFcallFZCB;
	labels(i) = "Call on foreign ZCB";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIFZCB(new MCPayoffList);
	calloptionIFZCB->push_back(A_IFZCBpayoff);
	calloptionIFZCB->push_back(B_IFZCBpayoff,-strikeIFZCB);
	mclist.push_back(calloptionIFZCB);
	CFvalues(i) = CFIFZCBcall;
	MBCFvalues(i) = MBCFcallIFZCB;
	labels(i) = "Intermediate maturity call on foreign ZCB";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionQ(new MCPayoffList);
	calloptionQ->push_back(A_Qpayoff);
	calloptionQ->push_back(B_Qpayoff,-strikeQ);
	mclist.push_back(calloptionQ,fxspot);
	CFvalues(i) = CFquanto;
	MBCFvalues(i) = MBCFquanto;
	labels(i) = "Quanto caplet";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionIQ(new MCPayoffList);
	calloptionIQ->push_back(A_IQpayoff);
	calloptionIQ->push_back(B_IQpayoff,-strikeIQ);
	mclist.push_back(calloptionIQ,fxspot);
	CFvalues(i) = CFIquanto;
	MBCFvalues(i) = MBCFIquanto;
	labels(i) = "Intermediate maturity quanto caplet";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionQS(new MCPayoffList);
	calloptionQS->push_back(A_QSpayoff);
	calloptionQS->push_back(B_QSpayoff,-strikeQS);
	mclist.push_back(calloptionQS,fxspot);
	CFvalues(i) = MBCFSquanto;
	MBCFvalues(i) = MBCFSquanto;
	labels(i) = "Foreign asset quanto";
	i++;
	boost::shared_ptr<MCPayoffList> calloptionQZCB(new MCPayoffList);
	calloptionQZCB->push_back(A_QZCBpayoff);
	calloptionQZCB->push_back(B_QZCBpayoff,-strikeFZCB);
	mclist.push_back(calloptionQZCB,fxspot);
	CFvalues(i) = MBCFZCBquanto;
	MBCFvalues(i) = MBCFZCBquanto;
	labels(i) = "Foreign ZCB quanto";
	i++;
	mclist.push_back(asianF_option.get_payoff());
	CFvalues(i) = MBCFasianF;
	MBCFvalues(i) = MBCFasianF;
	labels(i) = "Foreign DiscreteGeometricMeanFixedStrike";
	MCMapping<GaussMarkovWorld,Array<double,2> > mc_mapping(mclist,world,*(world.get_economies()[0]->initialTS),numeraire_index);
	boost::function<Array<double,1> (Array<double,2>)> func = boost::bind(&MCMapping<GaussMarkovWorld,Array<double,2> >::mappingArray,&mc_mapping,_1);
	//std::cerr << "Number of steps: " << world.number_of_steps() << endl;
	std::cerr << "Number of factors: " << world.factors() << endl;
	std::cerr << "Number of steps: " << world.number_of_steps() << endl;
    RandomArray<ranlib::NormalUnit<double>,double> random_container(normalRNG,world.factors(),world.number_of_steps()); 
	MCGeneric<Array<double,2>,Array<double,1>,RandomArray<ranlib::NormalUnit<double>,double> > mc(func,random_container);
	//MCGatherer<double> mcgatherer;
    MCGatherer<Array<double,1> > mcgatherer(number_of_options);
	size_t n = minpaths;
	boost::math::normal normal;
	double d = boost::math::quantile(normal,0.95);
	for (i=0;i<number_of_options;i++) cout << "," << labels(i) << "," << CFvalues(i) << ",,,";
	cout << endl << "Paths";
	for (i=0;i<number_of_options;i++) cout << ",MC value,95% CI lower bound,95% CI upper bound,Std error,Std error against MBinary";
	cout << endl;
	while (mcgatherer.number_of_simulations()<maxpaths) {
	  mc.simulate(mcgatherer,n);
	  Array<double,1> mean(mcgatherer.mean());
	  Array<double,1> stddev(mcgatherer.stddev());
	  cout << mcgatherer.number_of_simulations();
	  bool MCpassed = true;
	  for (i=0;i<number_of_options;i++) {
		double stderror = (mean(i)-CFvalues(i))/stddev(i);
		double stderrorMB = (mean(i)-MBCFvalues(i))/stddev(i);
		cout << ',' << mean(i) << ',' << mean(i)-d*stddev(i) << ',' << mean(i)+d*stddev(i);
		cout << ',' << std::abs(stderror) << ',' << std::abs(stderrorMB);
		if (std::max(std::abs(stderror),std::abs(stderrorMB))>2.0) MCpassed = false; }
	  if (MCpassed) cout << ",PASSED";
	  else          cout << ",FAILED";
	  cout << endl; 
	  n = mcgatherer.number_of_simulations(); 
	  std::cerr << n << ' ' << std::flush;
	  if (MCpassed) std::cerr << "  PASSED" << endl;
	  else          std::cerr << "  FAILED" << endl; }

  } // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
