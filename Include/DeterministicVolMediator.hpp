/** \file  DeterministicVolMediator.hpp
    \brief Header file declaring class for a Mediator pattern for deterministic volatility function classes.
           Copyright 2006, 2011 by Erik Schlögl

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

#ifndef DETERMINISTICVOLMEDIATOR_HPP
#define DETERMINISTICVOLMEDIATOR_HPP
#include "DeterministicVol.hpp"
#include "ExponentialVol.hpp"

namespace quantfin {

class DeterministicVolMediator {
public:
  enum type { FLAT, EXPONENTIAL };
  static double volproduct_DeterministicAssetVolDiff(double t,double dt,
                                                     const DeterministicAssetVol& xv,
                                                     const DeterministicAssetVol& v1,
                                                     const DeterministicAssetVol& v2);  
  static double volproduct_ConstVol(double t,double dt,const Array<double,1>& lvl,const DeterministicAssetVol& v2);   
  static double volproduct_ExponentialVol(double t,double dt,
                                          const Array<double,1>& lvl,
                                          const Array<double,1>& decay,
                                          const ExponentialVol& v1,
                                          const DeterministicAssetVol& v2);
  static double bondvolproduct_ConstVol(double t,double dt,double bondmat,const Array<double,1>& lvl,const DeterministicAssetVol& xv);
  static double bondvolproduct_ExponentialVol(double t,double dt,double bondmat,const Array<double,1>& lvl,const Array<double,1>& decay,const ExponentialVol& this_vol,const DeterministicAssetVol& xv);
  static double bondbondvolproduct_ConstVol(double t,double dt,double T1,double T2,const Array<double,1>& lvl,const ConstVol& v1,const DeterministicAssetVol& xv);
  static double bondbondvolproduct_ExponentialVol(double t,double dt,double T1,double T2,const Array<double,1>& lvl,const Array<double,1>& decay,const ExponentialVol& this_vol,const DeterministicAssetVol& xv);
  static Array<double,1> z_bondintegral_ExponentialVol(double t,double dt,double T,const Array<double,1>& lvl,const Array<double,1>& decay,const ExponentialVol& v1,const DeterministicAssetVol& xmodel);
  static Array<double,1> z_volintegral_ExponentialVol(double t,double dt,double T,const Array<double,1>& lvl,const Array<double,1>& decay,const ExponentialVol& v1,const DeterministicAssetVol& fxvol);
  static double FwdFXexponential_ExponentialVol(double t,double dt,const Array<double,1>& dWj,const Array<double,1>& dZj0,const Array<double,1>& dZjk,const Array<double,1>& nu,const Array<double,1>& a,const ExponentialVol* this_vol,DeterministicAssetVol* xv,DeterministicAssetVol* fxvol);
  static double covar_ExponentialVol(int j,double t,double dt,double lvl,double decay,const ExponentialVol& this_vol,const DeterministicAssetVol& xv);
};

}

#endif
