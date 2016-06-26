/** \file  TSBinomialExample.cpp
    \brief Test and demonstration program.
           Copyright 2006, 2015 by Erik Schloegl 
     */

#include <iostream>
#include <cstdlib>
#include <boost/bind.hpp>
#include "GaussianHJM.hpp"
#include "TSBinomial.hpp"
#include "QFArrayUtil.hpp"
#include "TSPayoff.hpp"
#include "TSInstruments.hpp"
#include "ExponentialVol.hpp"

using namespace quantfin; 
 
inline double positivePart(double x)
{
  return (x>0.0) ? x : 0.0;      
}

/** Test and demonstration for HJM implementation.

    Command-line arguments:
      -# interest rate. 
      -# volatility.
      -# number of time steps.
      -# maturity.
      -# moneyness.
  */
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;

  int i,j;
  try {
    double r = 0.05;
    if (argc>1) r = atof(argv[1]);
    double sgm = 0.3;
    if (argc>2) sgm = atof(argv[2]);
    int N = 10;
    if (argc>3) N = atoi(argv[3]);
    double mat = 1.5;
    if (argc>4) mat = atof(argv[4]);
    double K = 1.0;
    if (argc>5) K = atof(argv[5]);
    FlatTermStructure flat_ts(r,0.0,mat+10.0);
    Array<double,1> T(9);
    double T1 = mat;
    double T2 = T1 + 5.0;
    firstIndex idx;
    Array<double,1> SaSoT(N+1);
    double dtSaSo = T2/N;
    SaSoT = idx * dtSaSo;
    int iseg = find_segment(T1,SaSoT);
    if (std::abs(T1-SaSoT(iseg))>std::abs(T1-SaSoT(iseg+1))) iseg++;
    T1 = SaSoT(iseg);
    double strike = K * flat_ts(T2)/flat_ts(T1);
    double ZCBstrike = strike;
    /// Test SaSo binomial method
    cout << "Creating SaSo lattice... " << SaSoT << endl;
    TSBinomialMethod SaSo(flat_ts,sgm,SaSoT); 
    cout << "Rolling back term structures... " << endl;
    SaSo.rollbackTermstructures();
    cout << "SaSo verified: " << SaSo.verify() << endl;
    ZCBoption zcbpayoff(ZCBstrike,T2);
    boost::function<double (const TermStructure& ts)> f;
    f = boost::bind(std::mem_fun(&ZCBoption::operator()),&zcbpayoff,_1);
    SaSo.apply_payoff(iseg,f);
    SaSo.rollback(iseg,0);
    cout << "SaSo lattice ZCB call: " << SaSo.result() << endl;
    SaSo.apply_payoff(iseg,f);
    SaSo.rollback(iseg);
    cout << "SaSo lattice ZCB call using state prices: " << SaSo.result() << endl;
    cout << "Attempting to reproduce Figure 8.4 of Clewlow/Strickland..." << endl;
    Array<double,1> CST(5);
    CST = idx;
    cout << "Time line: " << CST << endl;
    Array<double,1> PDB(5);
    PDB(0) = 1.0;
    for (i=1;i<5;i++) PDB(i) = PDB(i-1)/1.05;
    cout << "Initial bonds: " << PDB << endl;
    TSLogLinear CSts(CST,PDB);
    TSBinomialMethod CSmodel(CSts,0.1,CST);
	cout << "Short rates:" << endl;
	for (i=0;i<4;i++) {
	  for (j=0;j<=i;j++) {
		cout << CSmodel.short_rate(i,j) << ','; }
	  cout << endl; }
	cout << "State prices:" << endl;
	for (i=0;i<4;i++) {
	  for (j=0;j<=i;j++) {
		cout << CSmodel.state_price(i,j) << ','; }
	  cout << endl; }
    cout << "Testing calibration to caplets..." << endl;
    double mr = 0.1;
    ExponentialVol evol(sgm/10.0,mr);
    GaussianHJM emodel(&evol,&CSts);
    std::vector<TSEuropeanInstrument*> caplets,floorlets;    
    double delta    = CST(2)-CST(1);
    double lvl      = CSts.simple_rate(CST(1),delta);
	double floorlvl = lvl;
    Caplet caplet1(emodel.caplet(CST(1),delta,lvl),CST(0),CST(1),lvl,delta);
    cout << caplet1.price() << ' ' << flush;
    caplets.push_back(&caplet1);
    Floorlet floorlet1(-1.0,CST(0),CST(1),floorlvl,delta);
    floorlets.push_back(&floorlet1);
    delta = CST(3)-CST(2);
    lvl   = CSts.simple_rate(CST(2),delta);
    Caplet caplet2(emodel.caplet(CST(2),delta,lvl),CST(0),CST(2),lvl,delta);
    cout << caplet2.price() << ' ' << flush;
    caplets.push_back(&caplet2);
    Floorlet floorlet2(-1.0,CST(0),CST(2),floorlvl,delta);
    floorlets.push_back(&floorlet2);
    delta = CST(4)-CST(3);
    lvl   = CSts.simple_rate(CST(3),delta);
    Caplet caplet3(emodel.caplet(CST(3),delta,lvl),CST(0),CST(3),lvl,delta);
    cout << caplet3.price() << ' ' << flush;
    caplets.push_back(&caplet3);
    Floorlet floorlet3(-1.0,CST(0),CST(3),floorlvl,delta);
    floorlets.push_back(&floorlet3);
    cout << endl;
    CSmodel.calibrate(caplets);
	cout << "Compare model price with input price" << endl;
    std::vector<TSEuropeanInstrument*>::iterator iter;
    for (iter=caplets.begin();iter!=caplets.end();iter++) {
      cout << CSmodel.price(**iter) << ' ' << (*iter)->price() << endl; }
	// Price an interest rate floor
	double floor = 0.0;
    for (iter=floorlets.begin();iter!=floorlets.end();iter++) {
      floor += CSmodel.price(**iter); }
	cout << "Price of interest rate floor with floor level " << floorlvl << ": " << floor << endl;
	// Instantiate larger lattice
	int N2 = 130;
    Array<double,1> SaSoT2(N2+1);
    double dtSaSo2 = T2/N2;
    SaSoT2 = idx * dtSaSo2;
    cout << "Creating SaSo lattice... " << SaSoT2 << endl;
    TSBinomialMethod SaSo2(flat_ts,sgm,SaSoT2); 
    cout << "Rolling back term structures... " << endl;
    SaSo2.rollbackTermstructures();
    cout << "SaSo verified: " << SaSo2.verify() << endl;
	// Price a European swaption
	Swaption swaption(-1.0,0.0,SaSoT2(25),r,0.5,4);
	cout << "Maturity: " << SaSoT2(25) << endl;
	cout << "Swaption price: " << SaSo2.price(swaption) << endl;
	// Price a Bermudan swaption
    BermudanSwaption bermudan(-1.0,0.0,SaSoT2(25),r,0.5,4); 
	if (!subset(bermudan.maturity(),SaSoT2)) throw std::logic_error("Timeline mismatch");
	cout << "Maturity: " << SaSoT2(25) << endl;
	cout << "Bermudan Swaption price: " << SaSo2.price(bermudan) << endl;
	// Price a barrier caplet
	Caplet capletB(-1.0,0.0,SaSoT2(50),lvl,delta);
	BarrierInstrument barrier_caplet(-1.0,0.0,capletB,SaSoT2(Range(0,50)),0.8*lvl,delta,-1);
	cout << "Caplet price: " << SaSo2.price(capletB) << endl;
	cout << "Barrier caplet price: " << SaSo2.price(barrier_caplet) << endl;

	} // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}
