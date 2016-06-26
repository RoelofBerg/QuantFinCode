/** \file  MBinaryPayoffs.hpp
    \brief Header file declaring classes to represent payoffs for the "Quintessential Option Pricing Formula" 
           of Skipper and Buchen (2003).
           Copyright 2006, 2012, 2013 by Erik Schlögl

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

#ifndef MBINARYPAYOFFS_HPP
#define MBINARYPAYOFFS_HPP
#include "MBinary.hpp"

namespace quantfin {

class AssetBinary : public MBinaryPayoff {
public:
  inline AssetBinary(const BlackScholesAsset& asset,double strike,double now,double mat,const TermStructure& xts);
  inline AssetBinary(const GaussMarkovWorld& world,double strike,double now,double mat,int reportable_asset_index,int callput = 1);
};

class BondBinary : public MBinaryPayoff {
public:
  inline BondBinary(const BlackScholesAsset& asset,double strike,double now,double mat,const TermStructure& xts);
  inline BondBinary(const GaussMarkovWorld& world,double strike,double now,double mat,int reportable_asset_index,int callput = 1);
};

/** European option on a foreign asset, struck in foreign currency (but price expressed in domestic currency).
    Exchange rate is interpreted as a BlackScholesAsset paying a continuous dividend yield. 
	If interest rates are stochastic, volatilities of the asset and the exchange rate must be the volatilities of the forward prices. */
class ForeignOption : public MBinaryPayoff {
public:
  inline ForeignOption(const BlackScholesAsset& asset,
	                   const BlackScholesAsset& fx,
					   double strike,
					   double now,
					   double mat,
					   const TermStructure& xts,
					   bool assetpart = true,
					   int callput = 1);
  inline ForeignOption(const GaussMarkovWorld& world,
	                   double strike,
					   double now,
					   double mat,
					   int reportable_asset_index,
					   int reportable_FX_index,
					   bool assetpart = true,
					   int callput = 1);
};

inline ForeignOption::ForeignOption(const BlackScholesAsset& asset,const BlackScholesAsset& fx,double strike,double now,double mat,const TermStructure& xts,bool assetpart,int callput)
  : MBinaryPayoff(xts,2,2,1) 
{      
  underlying.push_back(&asset);
  underlying.push_back(&fx);
  timeline = now, mat;
  index    = 0, 1,
	         1, 1;      
  alpha    = 1.0;
  if (!assetpart) alpha(0) = 0.0;
  S        = callput;
  A        = 1.0, 0.0;
  a        = strike;
}

inline ForeignOption::ForeignOption(const GaussMarkovWorld& world,double strike,double now,double mat,int reportable_asset_index,int reportable_FX_index,bool assetpart,int callput)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,2,1) 
{      
  timeline = now, mat;
  index    = reportable_asset_index, 1,
	         reportable_FX_index, 1;      
  alpha    = 1.0;
  if (!assetpart) alpha(0) = 0.0;
  S        = callput;
  A        = 1.0, 0.0;
  a        = strike;
}

inline AssetBinary::AssetBinary(const BlackScholesAsset& asset,double strike,double now,double mat,const TermStructure& xts)
  : MBinaryPayoff(xts,2,1,1) 
{      
  underlying.push_back(&asset);
  timeline = now, mat;
  index    = 0, 1;      
  alpha    = 1.0;
  S        = 1;
  A        = 1.0;
  a        = strike;
}

inline AssetBinary::AssetBinary(const GaussMarkovWorld& world,double strike,double now,double mat,int reportable_asset_index,int callput)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,1,1) 
{      
  timeline = now, mat;
  index    = reportable_asset_index, 1;      
  alpha    = 1.0;
  S        = callput;
  A        = 1.0;
  a        = strike;
}

inline BondBinary::BondBinary(const BlackScholesAsset& asset,double strike,double now,double mat,const TermStructure& xts)
  : MBinaryPayoff(xts,2,1,1) 
{      
  underlying.push_back(&asset);
  timeline = now, mat;
  index    = 0, 1;      
  alpha    = 0.0;
  S        = 1;
  A        = 1.0;
  a        = strike;
}

inline BondBinary::BondBinary(const GaussMarkovWorld& world,double strike,double now,double mat,int reportable_asset_index,int callput)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,1,1) 
{      
  timeline = now, mat;
  index    = reportable_asset_index, 1;      
  alpha    = 0.0;
  S        = callput;
  A        = 1.0;
  a        = strike;
}

class AssetProductBinary : public MBinaryPayoff {
public:
  inline AssetProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index);
  inline AssetProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index,const Array<double,1>& xalpha);
};

inline AssetProductBinary::AssetProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,reportable_asset_index.extent(firstDim),1) 
{      
  int i;
  timeline = now, mat;
  for (i=0;i<reportable_asset_index.extent(firstDim);i++) {
	index(0,i)  = reportable_asset_index(i);      
	index(1,i)  = 1; }  
  alpha    = 1.0;
  S        = 1;
  A        = 1.0;
  a        = strike;
}

inline AssetProductBinary::AssetProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index,const Array<double,1>& xalpha)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,reportable_asset_index.extent(firstDim),1) 
{      
  int i;
  timeline = now, mat;
  for (i=0;i<reportable_asset_index.extent(firstDim);i++) {
	index(0,i)  = reportable_asset_index(i);      
	index(1,i)  = 1; }  
  alpha    = xalpha;
  S        = 1;
  A        = 1.0;
  a        = strike;
}

class BondProductBinary : public MBinaryPayoff {
public:
  inline BondProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index);
};

inline BondProductBinary::BondProductBinary(const GaussMarkovWorld& world,double strike,double now,double mat,const Array<int,1>& reportable_asset_index)
  : MBinaryPayoff(*world.get_economies()[0]->initialTS,2,reportable_asset_index.extent(firstDim),1) 
{      
  int i;
  timeline = now, mat;
  for (i=0;i<reportable_asset_index.extent(firstDim);i++) {
	index(0,i)  = reportable_asset_index(i);      
	index(1,i)  = 1; }  
  alpha    = 0.0;
  S        = 1;
  A        = 1.0;
  a        = strike;
}

}

#endif


