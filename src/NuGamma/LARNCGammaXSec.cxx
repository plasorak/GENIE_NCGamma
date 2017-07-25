//____________________________________________________________________________
/*
  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Pierre Lasorak <p.j.j.lasorak \at qmul.ac.uk>
  Queen Mary University, London

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________


#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "NuGamma/LARNCGammaXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "NuGamma/TensorInc.h"
#include "BaryonResonance/BaryonResonance.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace TensorUtils;

//____________________________________________________________________________
LARNCGammaXSec::LARNCGammaXSec() :
  XSecAlgorithmI("genie::LARNCGammaXSec")
{

}
//____________________________________________________________________________
LARNCGammaXSec::LARNCGammaXSec(string config) :
  XSecAlgorithmI("genie::LARNCGammaXSec", config)
{

}
//____________________________________________________________________________
LARNCGammaXSec::~LARNCGammaXSec()
{

}
//____________________________________________________________________________
double LARNCGammaXSec::XSec(const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;
  LOG("LARNCGammaXSec", pWARN)
    << "*** Calculating the cross section";
  //  this->CalcAmplitude();

  InitialState* initialstate =  interaction->InitStatePtr();
  Kinematics* kine = interaction->KinePtr();  
  (void)kine;
  double W =        interaction->KinePtr()->GetKV(kKVW);
  double Q2 =       interaction->KinePtr()->GetKV(kKVQ2);
  double EGamma =   interaction->KinePtr()->GetKV(kKVEGamma);
  double PhiGamma = interaction->KinePtr()->GetKV(kKVPhiGamma);
  LOG("LARNCGammaXSec", pWARN)  << "W        " << W;
  LOG("LARNCGammaXSec", pWARN)  << "Q2       " << Q2;
  LOG("LARNCGammaXSec", pWARN)  << "EGamma   " << EGamma;
  LOG("LARNCGammaXSec", pWARN)  << "PhiGamma " << PhiGamma;
  
  
  gNeutrinoInit = (TLorentzVector*)initialstate->GetProbeP4()->Clone();
  gTargetInit   = (TLorentzVector*)initialstate->GetTgtP4()->Clone();
 

  //TensorOrder4& output1, output2;
  //AmplitudeNum(interaction, output1, output2);
  

  return 1; // For now return 1
}
//_____________________________________________________________________________
double LARNCGammaXSec::Integral(const Interaction * interaction) const
{
  // Compute the cross section for a free fNucleon target
  const InitialState & init_state = interaction -> InitState();
  const Target &       target     = init_state.Tgt();
  (void)target;
  double Ev   = init_state.ProbeE(kRfHitNucRest);
  (void)Ev;
  double xsec = 1;

  /*
    Integrate over all the four fold cross section
   */

  return xsec;
}

//____________________________________________________________________________
bool LARNCGammaXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  //if(interaction->ProcInfo().IsRESNCGamma()) return true;
  return false;
}

//____________________________________________________________________________
bool LARNCGammaXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void LARNCGammaXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LARNCGammaXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LARNCGammaXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  (void)gc;
  (void)confp;
  //fGw = fConfig->GetDoubleDef("Gw", gc->GetDouble("RESNCGamma-Gw"));
}


//____________________________________________________________________________
void LARNCGammaXSec::ElectroMagneticFormFactor(double t, int fNucleon, double &f1, double &f2)
{
  double gd = TMath::Power(1 - t/fgVMassSq, -2);
  double gmp = fgPMagMom * gd;
  double gmn = fgNMagMom * gd;

  double tao = t / (4 * fgMn * fgMn);
  double gep = gd;
  double gen = fgNMagMom * (fgAp * tao / (1 - fgBp * tao)) * gd;
  double ge, gm;

  if(fNucleon == 1){
    ge = gep;
    gm = gmp;
  }else{
    ge = gen;
    gm = gmn;
  }

  f1 = (ge - tao * gm) / (1 - tao);
  f2 = (gm - ge)       / (1 - tao);
}
//____________________________________________________________________________
void LARNCGammaXSec::NeutralCurrentFormFactor(double t, int fNucleon, double &fv1, double &fv2, double &fva)
{
  double tao = t / (4 * fgMn * fgMn);

  double gd = TMath::Power(1 - t/fgVMassSq, -2);
  double gmp = fgPMagMom * gd;
  double gmn = fgNMagMom * gd;

  double gep = gd;
  double gen = fgNMagMom * (fgAp * tao/(1 - fgBp*tao))*gd;

  double f1_p = (gep - tao*gmp) / (1-tao);
  double f2_p = (gmp - gep)     / (1-tao);

  double f1_n = (gen - tao*gmn) / (1-tao);
  double f2_n = (gmn - gen)     / (1-tao);

  double fa = fgGa / TMath::Power(1 - t/fgNuclMa, 2);

  double f1_s = -fgF1s_0 * t / (1 - tao) / TMath::Power(1 - t / fgMv2, 2);
  double f2_s =  fgF2s_0     / (1 - tao) / TMath::Power(1 - t / fgMv2, 2);
  double fa_s =  fgDeltaS / TMath::Power(1 - t / fgNuclMa, 2);

  double stheta = 1 - 4 * fgSinThetaW2;

  if(fNucleon == 1){
    fv1=(stheta * f1_p - f1_n - f1_s)/2;
    fv2=(stheta * f2_p - f2_n - f2_s)/2;
    fva=(fa + fa_s)/2;
  }else if (fNucleon == -1){
    fv1=(stheta * f1_n - f1_p - f1_s)/2;
    fv2=(stheta * f2_n - f2_p - f2_s)/2; 
    fva=(-fa + fa_s)/2;
  }
}
//____________________________________________________________________________
void LARNCGammaXSec::P33_1232FormFactor(double q2, double &c3v, double &c4v, double &c5v, double &c3a, double &c4a, double &c5a, double &c6a)
{
  // using MAID2007 in the vector sector
  // IMPORTANT: these are the CC form factors: CiV and CiA in Tables 5.5 and 5.6 in Tina's thesis 
  // Rules to obtain EM and NC form factors can be found in these tables

  // Q2                  4-momentum transfer 
  // c3v,c4v,c5v        vector form factors (dimensionless)
  // c3a,c4a,c5a,c6a    axial form factors  (dimensionless)

  // real*8,parameter::fgAlpha=1./137.     fine structure constant
  // real*8,parameter::md=1.232           delta mass
  // real*8,parameter::mn=.939           fNucleon mass
  // real*8,parameter::mpi=.139          pion mass

  // real*8,parameter::Ma=0.93     N-Delta axial mass

  // real*8 ::a12,a32,s12          helicity amplitudes
  // real*8 ::r,fgr,mpw2,mmw2

  // real*8 ::Fd   dipole

  // -------------
  // vector form factors          
  double a12, a32, s12;
  HelicityAmplitude(kP33_1232,"p", q2, a12, a32, s12); //   fNucleon label 'p' irrelevant in this case
  //  a12=a12*1.e-3
  //  a32=a32*1.e-3
  //  s12=s12*1.e-3      (-) gives the solution like in Tina's thesis

  double kr = (TMath::Power(fgMdelta, 2) - TMath::Power(fgMnucl, 2)) / 2 / fgMdelta;
  double mpw2 = TMath::Power(fgMnucl + fgMdelta, 2);
  double mmw2 = TMath::Power(fgMnucl - fgMdelta, 2);

  double r = TMath::Sqrt(6 * kr / kPi / fgAlpha * fgMnucl * fgMdelta / (mmw2 - q2));

  c3v = -r * fgMnucl * fgMdelta / (mpw2 - q2) * (a12 + a32 / kSqrt3);
  //c3v=-r * mn     * md      / (mpw2 - q2) * (a12 + a32 / sqrt(3.))
  c4v = -2 * r * TMath::Power(fgMnucl, 2) / (mpw2 - q2) / (mmw2 - q2) * (fgMnucl * fgMdelta * a12
									- (TMath::Power(fgMdelta, 2) + TMath::Power(fgMnucl, 2) - fgMnucl * fgMdelta - q2) * a32 / kSqrt3
									+ kSqrt2 * TMath::Power(fgMdelta, 2) * (TMath::Power(fgMdelta, 2) - TMath::Power(fgMnucl, 2) - q2) / TMath::Sqrt((mpw2 - q2) * (mmw2 - q2)) * s12);
  //c4v=-2.*r *mn**2                     / (mpw2 - q2) / (mmw2 - q2) * (mn     * md      * a12 -
  //                                                                   (md**2                     + mn**2                   - mn     * md      - q2) * a32 / sqrt(3.) +
  //                                                                   sqrt(2.)* md**2                    * (md**2             - mn**2                   - q2) / sqrt(       (mpw2 - q2) * (mmw2 - q2)) * s12)

  c5v = -2. * r * TMath::Power(fgMdelta*fgMnucl, 2) / (mpw2 - q2) / (mmw2 - q2) * (-a12 + a32 / kSqrt3 - kSqrt2 * (TMath::Power(fgMdelta, 2) - TMath::Power(fgMnucl, 2) + q2) / TMath::Sqrt((mpw2 - q2) * (mmw2 - q2)) * s12);
  //c5v = -2* r * (md*mn)**2                      / (mpw2 - q2) / (mmw2 - q2) * (-a12 + a32 / sqrt(3.)-sqrt(2.)*(md**2                    - mn**2                   +q2 ) / sqrt(       (mpw2 - q2) * (mmw2 - q2)) * s12)
  
  //        axial form factors
  //        simple dipole for c5a and adler
  
  double Fd = TMath::Power(1 - q2 / TMath::Power(fgNDeltaMa, 2), -2);
  c5a = 1.2 * Fd;
  c4a = -c5a / 4.;
  c3a = 0.;
  c6a = c5a * TMath::Power(fgMnucl, 2)/(TMath::Power(fgMpi, 2) - q2);
}
//____________________________________________________________________________
void LARNCGammaXSec::P11_1440FormFactor(double t1, TensorOrder2& fem, TensorOrder2& fvem, TensorOrder2& fvnc, TensorOrder1& fanc)
{
  //Form Factor For P11(1440)   
  //real*8,dimension(2,6) :: fem,fvem,fvnc
  //real*8,dimension(2) :: fanc

  //xma=1    
  double smr = fgMP11;  //1.462
  double fa0 = fgFa0P11; //-0.52

  double rk = TMath::Sqrt(fgMn * (TMath::Power(smr, 2) - fgMnSq) / kPi / fgAlpha);
  double xmp2 = TMath::Power(smr + fgMn, 2);
  double xmm2 = TMath::Power(smr - fgMn, 2);
  double Q2 = -t1;
  double xmpq2 = TMath::Power(smr + fgMn, 2) + Q2;
  double xmmq2 = TMath::Power(smr - fgMn, 2) + Q2;

  //   a12=a12*1.e-3
  //   a32=a32*1.e-3
  //   s12=s12*1.e-3      (-) gives the solution like in Tina's thesis
  double a12, a32, s12;
  //   EM form factor 
  double tem=0;
  HelicityAmplitude(kP11_1440, "p", tem, a12, a32, s12);
  fem(0,0) = kSqrt8 * fgMnSq * rk / xmp2 / TMath::Sqrt(xmm2) * (a12 - TMath::Sqrt(8 / xmp2 / xmm2) * smr * (smr + fgMn) * s12);
  fem(0,1) = kSqrt2 * fgMn   * rk / xmp2 / TMath::Sqrt(xmm2) * (smr + fgMn) * a12;

  HelicityAmplitude(kP11_1440, "n", tem, a12, a32, s12);
  fem(1,0) = kSqrt8 * fgMnSq * rk / xmp2 / TMath::Sqrt(xmm2)* (a12 - TMath::Sqrt(8 / xmp2 / xmm2) * smr * (smr + fgMn) * s12);
  fem(1,1) = kSqrt2 * fgMn   * rk / xmp2 / TMath::Sqrt(xmm2)* (smr + fgMn) * a12;

  //  NC form factor 
  HelicityAmplitude(kP11_1440,"p", t1, a12, a32, s12);
  fvem(0,0) = kSqrt8 * fgMnSq * rk / xmpq2 / TMath::Sqrt(xmmq2) * (a12 - TMath::Sqrt(8 / xmpq2 / xmmq2) * smr * (smr + fgMn) * s12);
  fvem(0,1) = kSqrt2 * fgMn   * rk / xmpq2 / TMath::Sqrt(xmmq2) * ((smr + fgMn) * a12 + TMath::Sqrt(8/xmpq2/xmmq2)*smr*Q2*s12 );
  
  HelicityAmplitude(kP11_1440,"n", t1, a12, a32, s12);
  fvem(1,0) = kSqrt8 * fgMnSq * rk / xmpq2 / TMath::Sqrt(xmmq2) * (a12 - TMath::Sqrt(8 / xmpq2 / xmmq2) * smr * (smr + fgMn) * s12);
  fvem(1,1) = kSqrt2 * fgMn   * rk / xmpq2 / TMath::Sqrt(xmmq2) * ((smr + fgMn) * a12 + TMath::Sqrt(8 / xmpq2 / xmmq2) * smr * Q2 * s12);
  
  EMtoNCFormFactor(2, fvem, fvnc);

  fanc(0)= fa0 * TMath::Power(1 + Q2 / fg2Ma2, -2) * 0.5;
  fanc(1) = -fanc(0);
}
//____________________________________________________________________________
void LARNCGammaXSec::S11_1535FormFactor(double t1, TensorOrder2& fem, TensorOrder2& fvem, TensorOrder2& fvnc, TensorOrder1& fanc)
{
  //Form Factor For S11(1535)   
  double smr = fgMS11;  // 1.535
  double fa0 = fgFa0S11; // -0.23
  
  double rk = TMath::Sqrt(fgMn * (TMath::Power(smr, 2) - fgMnSq) / kPi / fgAlpha);
  double xmp2 = TMath::Power(smr + fgMn, 2);
  double xmm2 = TMath::Power(smr - fgMn, 2);
  double Q2 = -t1;
  double xmpq2 = TMath::Power(smr + fgMn, 2) + Q2;
  double xmmq2 = TMath::Power(smr - fgMn, 2) + Q2;
  
  
  //   EM form factor 
  double tem = 0, a12, a32, s12;
  HelicityAmplitude(kS11_1535, "p", tem, a12, a32, s12);
  fem(0,0) = kSqrt8 * fgMnSq * rk / xmm2 / TMath::Sqrt(xmp2) * (a12 + TMath::Sqrt(8. / xmp2 / xmm2) * smr * (smr - fgMn) * s12);
  fem(0,1) = kSqrt2 * fgMn   * rk / xmm2 / TMath::Sqrt(xmp2) * (smr - fgMn) * a12;

  HelicityAmplitude(kS11_1535, "n", tem, a12, a32, s12);
  fem(1,0) = kSqrt8 * fgMnSq * rk / xmm2 / TMath::Sqrt(xmp2) * (a12 + TMath::Sqrt(8. / xmp2 / xmm2) * smr * (smr - fgMn) * s12);
  fem(1,1) = kSqrt2 * fgMn   * rk / xmm2 / TMath::Sqrt(xmp2) * (smr - fgMn) * a12;

  //   NC form factor 
  HelicityAmplitude(kS11_1535, "p", t1, a12, a32, s12);
  fvem(0,0) = kSqrt8 * fgMnSq * rk / xmmq2 / TMath::Sqrt(xmpq2) * (a12 + TMath::Sqrt(8. / xmpq2 / xmmq2) * smr * (smr - fgMn) * s12);
  fvem(0,1) = kSqrt2 * fgMn   * rk / xmmq2 / TMath::Sqrt(xmpq2) * ((smr - fgMn) * a12 - TMath::Sqrt(8. / xmpq2 / xmmq2) * smr * Q2 * s12);

  HelicityAmplitude(kS11_1535, "n", t1, a12, a32, s12);
  fvem(1,0) = kSqrt8 * fgMnSq * rk / xmpq2 / TMath::Sqrt(xmmq2) * (a12 + TMath::Sqrt(8. / xmpq2 / xmmq2) * smr * (smr - fgMn) * s12);
  fvem(1,1) = kSqrt2 * fgMn   * rk / xmpq2 / TMath::Sqrt(xmmq2) * ((smr - fgMn) * a12 - TMath::Sqrt(8. / xmpq2 / xmmq2) * smr * Q2 * s12);
  
  EMtoNCFormFactor(2, fvem, fvnc);

  fanc(0) = fa0 * TMath::Power(1 + Q2 / fg2Ma2, -2) * 0.5;
  fanc(1) = -fanc(0);
}
//____________________________________________________________________________
void LARNCGammaXSec::D13_1520FormFactor(double t1, TensorOrder2& fem, TensorOrder2& fvem, TensorOrder2& fvnc, TensorOrder2& fanc)
{
  //Form Factor For D13(1520)   
  //  real*8,dimension(2,6) :: fem,fvem,fvnc,fanc
  
  double fca50 = fgFca5D13;// //-2.15
  double smr = fgMD13; //1.524
  //    fg2Ma2=1

  double rk = TMath::Sqrt(fgMn * (TMath::Power(smr, 2) - fgMnSq) / kPi / fgAlpha);
  double xmp2 = TMath::Power(smr + fgMn, 2);
  double xmm2 = TMath::Power(smr - fgMn, 2);
  double Q2 = -t1;
  double xmpq2 = TMath::Power(smr + fgMn, 2) + Q2;
  double xmmq2 = TMath::Power(smr - fgMn, 2) + Q2;

  // EM form factor 
  double tem = 0;
  double a12, a32, s12;
  HelicityAmplitude(kD13_1520, "p", tem, a12, a32, s12);

  fem(0,2) = rk/TMath::Sqrt(xmp2)*fgMn*smr/xmm2*(kSqrt3*a12-a32);
  fem(0,3) = 2*rk/xmp2/TMath::Sqrt(xmp2)*fgMnSq/xmm2 *(-kSqrt3*smr*fgMn*a12 +(TMath::Power(smr, 2)+smr*fgMn+fgMnSq)*a32 - kSqrt6*TMath::Power(smr,2)*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12 );
  fem(0,4)=2*rk/xmp2/TMath::Sqrt(xmp2)*TMath::Power(smr,2)*fgMnSq/xmm2 *( -kSqrt3*a12 - a32 + kSqrt6*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12 );

  HelicityAmplitude(kD13_1520,"n",tem,a12,a32,s12);

  fem(1,2)=rk/TMath::Sqrt(xmp2)*fgMn*smr/xmm2*(kSqrt3*a12-a32);
  fem(1,3)=2*rk/xmp2/TMath::Sqrt(xmp2)*fgMnSq/xmm2 *( -kSqrt3*smr*fgMn*a12 +(TMath::Power(smr,2)+smr*fgMn+fgMnSq)*a32 - kSqrt6*TMath::Power(smr,2)*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12);
  fem(1,4)=2*rk/xmp2/TMath::Sqrt(xmp2)*TMath::Power(smr,2)*fgMnSq/xmm2 *( -kSqrt3*a12 - a32 + kSqrt6*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12);

  //  NC form factor
  tem=t1;
  HelicityAmplitude(kD13_1520,"p",tem,a12,a32,s12);
  fvem(0,2)=rk / TMath::Sqrt(xmpq2)*fgMn*smr/xmmq2*(kSqrt3*a12-a32);
  fvem(0,3)=2 * rk / xmpq2 / TMath::Sqrt(xmpq2) * fgMnSq / xmmq2 *(-kSqrt3 * smr * fgMn * a12 + (TMath::Power(smr, 2) + smr * fgMn + fgMnSq + Q2) * a32
								   -kSqrt6 * TMath::Power(smr, 2) * (TMath::Power(smr, 2) - fgMnSq + Q2) / TMath::Sqrt(xmpq2 * xmmq2 ) * s12);
  fvem(0,4)=2 * rk / xmpq2 / TMath::Sqrt(xmpq2) * TMath::Power(smr, 2) * fgMnSq / xmmq2 *(-kSqrt3 * a12 - a32 + kSqrt6 * (TMath::Power(smr, 2) - fgMnSq - Q2) / TMath::Sqrt(xmpq2 * xmmq2) * s12);

  HelicityAmplitude(kD13_1520,"n",tem,a12,a32,s12);
  fvem(1,2) =     rk / TMath::Sqrt(xmpq2) * fgMn * smr / xmmq2 * (kSqrt3 * a12 - a32);
  fvem(1,3) = 2 * rk / xmpq2 / TMath::Sqrt(xmpq2) * fgMnSq / xmmq2 * (-kSqrt3 * smr * fgMn * a12 +
								      (TMath::Power(smr, 2) + smr * fgMn + fgMnSq + Q2) * a32
								      -kSqrt6*TMath::Power(smr, 2) * (TMath::Power(smr, 2) - fgMnSq + Q2) / TMath::Sqrt(xmpq2 * xmmq2) * s12);
  fvem(1,4) = 2 * rk / xmpq2 / TMath::Sqrt(xmpq2) * TMath::Power(smr, 2) * fgMnSq / xmmq2 * (-kSqrt3 * a12 - a32
											     + kSqrt6 * (TMath::Power(smr, 2) - fgMnSq - Q2) / TMath::Sqrt(xmpq2 * xmmq2) * s12);

  EMtoNCFormFactor(3,fvem,fvnc);
  fanc(0,2) = 0;
  fanc(0,4) = fca50 * TMath::Power(1 + Q2 / fg2Ma2,-2) * 0.5;
  fanc(0,5) = fgMnSq / (Q2 + TMath::Power(fgMpi2,2)) * fanc(0,4);
  fanc(0,3) = -0.25 * fanc(0,4);

  for(int i =3; i < 7; i++){
    fanc(1,i) = -fanc(0,i);
  }
  // Luis choice
  fanc(0,3) = 0;
  fanc(1,3) = 0;
}
//____________________________________________________________________________
void LARNCGammaXSec::P33_1232FormFactor(double t1, TensorOrder2& fem, TensorOrder2& fvem, TensorOrder2& fvnc, TensorOrder2& fanc)
{
  // Form Factor For P33(1232)   
  //real*8,dimension[1][6) :: fem,fvem,fvnc,fanc
  
  double fca50 = fgFca5P33; //  1.2
  double smr = fgMdelta;   // 1.232

  //double sw = 1 - 2 * fgSinThetaW2;

  double rk = TMath::Sqrt(fgMn*(TMath::Power(smr,2)-fgMnSq)/kPi/fgAlpha);
  double xmp2 = TMath::Power(smr+fgMn,2);
  double xmm2 = TMath::Power(smr-fgMn,2);
  double Q2 = -t1;
  double xmpq2 = TMath::Power(smr+fgMn,2) + Q2;
  double xmmq2 = TMath::Power(smr-fgMn,2) + Q2;

  // EM form factor 
  double tem=0;
  double a12,a32,s12;
  HelicityAmplitude(kP33_1232,"p",tem,a12,a32,s12);

  fem(0,2)=rk/TMath::Sqrt(xmm2)*fgMn*smr/xmp2*(kSqrt3*a12+a32);
  fem(0,3)=2*rk/xmm2/TMath::Sqrt(xmm2)*fgMnSq/xmp2 *( kSqrt3*smr*fgMn*a12 -(TMath::Power(smr, 2) - smr*fgMn+fgMnSq)*a32 + kSqrt6*TMath::Power(smr,2)*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12);
  fem(0,4)=2*rk/xmm2/TMath::Sqrt(xmm2)*TMath::Power(smr,2)*fgMnSq/xmp2 *(-kSqrt3*a12 + a32 - kSqrt6*(TMath::Power(smr,2)-fgMnSq)/TMath::Sqrt(xmp2*xmm2)*s12);

  // NC form factor
  tem=t1;
  HelicityAmplitude(kP33_1232,"p",tem,a12,a32,s12);

  fvem(0,2)=rk/TMath::Sqrt(xmmq2)*fgMn*smr/xmpq2*(kSqrt3*a12+a32);
  fvem(0,3)=2*rk/xmmq2/TMath::Sqrt(xmmq2)*fgMnSq/xmpq2 *( kSqrt3*smr*fgMn*a12 -(TMath::Power(smr,2)-smr*fgMn+fgMnSq+Q2)*a32 + kSqrt6*TMath::Power(smr,2)*(TMath::Power(smr,2)-fgMnSq+Q2)/TMath::Sqrt(xmpq2*xmmq2)*s12);
  fvem(0,4)=2*rk/xmmq2/TMath::Sqrt(xmmq2)*TMath::Power(smr,2)*fgMnSq/xmpq2 *( -kSqrt3*a12 + a32 - kSqrt6*(TMath::Power(smr,2)-fgMnSq-Q2)/TMath::Sqrt(xmpq2*xmmq2)*s12);

  for(int i = 3; i < 6; i++)
    fvnc(0,i)=(1-2*fgSinThetaW2)*fvem(0,i);
  
  fanc(0,2) = 0;
  fanc(0,4) = fca50 * TMath::Power(1 + Q2 / fgMaSq, -2) * (-1);//      Fnc = -Fcc
  fanc(0,6) = fgMnSq / (Q2 + TMath::Power(fgMpi2, 2)) * fanc(0,4);
  fanc(0,3) = -0.25 * fanc(0,4);
}
//____________________________________________________________________________
void LARNCGammaXSec::HelicityAmplitude(Resonance_t res, string nucl, double q2, double &A12, double &A32, double &S12)
{
  // Helicity amplitude
  // Helicity amplitides for fNucleon-resonance em transitions according to MAID2007
  // character(len=9),intent(in)::res   P33(1232),P11(1440), etc
  // character(len=1),intent(in)::nucl  p or n 
  // real*8,intent(in)::q2   in GeV2
  // real*8,intent(out)::A12,A32,S12   helicity ampls. in units of 10-3 GeV^(-1/2)

  // real*8,parameter::pi=3.1416
  // real*8,parameter::mn=0.9382723

  // real*8 ::Q2G    -q2
  // real*8 ::mres    resonance Breit-Wigner mass
  // real*8 ::fggcm0   energy of an equivalent real photon in cm
  // real*8 ::egcm,qcm     photon energy and momentum in cm

  // real*8 ::AE,AM,AS,Fq,tau

  double Q2G = -q2;
  switch(res){
  case(kP33_1232) :
    {
      double mres  = fgMdelta;
      double kgcm0 = (TMath::Power(mres, 2) - TMath::Power(fgMnucl, 2)) / 2. / mres;
      double egcm  = (TMath::Power(mres, 2) - Q2G - TMath::Power(fgMnucl, 2)) / 2. / mres;
      double qcm   = TMath::Sqrt(TMath::Power(egcm, 2) + Q2G);
      double tau   = Q2G / 4. / TMath::Power(fgMnucl, 2);
      double Fq    = 1./TMath::Power(1 + Q2G / 0.71, 2) * (qcm / kgcm0);
      
      double AM = 300    * (1. + 0.01  * Q2G) * TMath::Exp(-0.23 * Q2G) * Fq;
      double AE = -6.37  * (1. - 0.021 * Q2G) * TMath::Exp(-0.16 * Q2G) * Fq;
      double AS = -12.40 * (1. + 0.12  * Q2G) / (1. + 4.9 * tau) * ( qcm / kgcm0) * TMath::Exp(-0.23 * Q2G) * Fq;
      
      A12 = -.5 * (3 * AE + AM);
      A32 = kSqrt3 / 2. * (AE - AM);
      S12 = -kSqrt2 * AS;
    }
  case(kP11_1440) :
    {
      if(nucl == "p"){
	A12 = -61.4 *(1. - 1.22 * Q2G - 0.55 * TMath::Power(Q2G,4)) * TMath::Exp(-1.51 * Q2G);
	S12 = 4.2   *(1. + 40.  * Q2G + 1.5  * TMath::Power(Q2G,4)) * TMath::Exp(-1.75 * Q2G);
	A32 = 0.;
      }else{
	A12 = 54.1  * (1. + 0.95 * Q2G) * TMath::Exp(-1.77 * Q2G);
	S12 = -41.5 * (1. + 2.98 * Q2G) * TMath::Exp(-1.55 * Q2G);
	A32 = 0.;
      }
    }
  case(kD13_1520) :
    {
      if(nucl == "p"){
	A12 = -27.  * (1. + 7.77 * Q2G) * TMath::Exp(-1.09 * Q2G);
	A32 = 161.  * (1. + 0.69 * Q2G) * TMath::Exp(-2.1  * Q2G);
	S12 = -63.6 * (1. + 4.19 * Q2G) * TMath::Exp(-3.40 * Q2G);
      }else{
	A12 = -77.  * (1. - 0.53 * Q2G) * TMath::Exp(-1.55 * Q2G);
	A32 = -154. * (1. + 0.58 * Q2G) * TMath::Exp(-1.75 * Q2G);
	S12 = 13.6  * (1. + 15.7 * Q2G) * TMath::Exp(-1.57 * Q2G);
      }
    }
  case(kS11_1535) :
    {
      if(nucl =="p"){
	A12 = 66. * (1. + 1.61 * Q2G) * TMath::Exp(-0.70 * Q2G);
	S12 = -2. * (1. + 23.9 * Q2G) * TMath::Exp(-0.81 * Q2G);
	A32 = 0.;
      }else{
	A12 = -51. * (1. + 4.75 * Q2G) * TMath::Exp(-1.69 * Q2G);
	S12 = 28.5 * (1. + 0.36 * Q2G) * TMath::Exp(-1.55 * Q2G);
	A32 = 0.;
      }
    }
  default :
    LOG("LARNCGammaXSec", pFATAL) << "Resonance isn't handled";
  }
  A12 = A12 * 1.E-3;
  A32 = A32 * 1.E-3;
  S12 = S12 * 1.E-3; 
}
//____________________________________________________________________________
void LARNCGammaXSec::EMtoNCFormFactor(double n, TensorOrder2&fvem, TensorOrder2& fvnc)
{
  // Transfer from EM FF to NC FF
  //  real*8,dimension[1][6) :: fvem,fvnc
  double sw = 0.5 - 2 * fgSinThetaW2;
  double m1 = 0, m2 = 0;
  
  if(n == 2){
    m1 = 1;
    m2 = 2;
  }else if(n == 3){
    m1 = 3;
    m2 = 5;
  }
  
  for(int i = m1; i < m2; i++){
    fvnc(0,i) = sw * fvem(0,i) - 0.5 * fvem(1,i);
    fvnc(1,i) = sw * fvem(1,i) - 0.5 * fvem(0,i);
  }
}
//____________________________________________________________________________
double LARNCGammaXSec::DeltaPropagator(TensorOrder1&am)
{
  //D(Q)=1/(Q**2-M**2)  
  double temp = SMomentum(am);
  return 1./(temp - fgMnSq);
}
//____________________________________________________________________________
TComplex LARNCGammaXSec::Propagator(Resonance_t res, TensorOrder1&sp)
{
  //The propagator of the excited resonances 
  double smr;
  if     (res == kNoResonance) smr = LARNCGammaXSec::fgMn;
  else if(res == kP33_1232)    smr = LARNCGammaXSec::fgMdelta;
  else if(res == kP11_1440)    smr = LARNCGammaXSec::fgMP11;
  else if(res == kD13_1520)    smr = LARNCGammaXSec::fgMD13;
  else if(res == kS11_1535)    smr = LARNCGammaXSec::fgMS11;
  else LOG("LARNCGammaXSec", kERROR) << "Unsupported resonance!";
  
  double sp2 = DotProdMetric(sp, sp);

  TComplex cu(0,1);
  TComplex prop;
  if(res == kP33_1232){
    // call medium(sp2,cprop)          for the Delta, we use the deltamedium.f90
    prop = LARNCGammaXSec::cDdelta(sp);
  }else{
    double wid = LARNCGammaXSec::Width(res, sp2);
    prop = 1. / (sp2 - TMath::Power(smr, 2) + cu * wid * smr);
  }
  return prop;
}
//____________________________________________________________________________
void LARNCGammaXSec::TraceLight(TensorOrder2& c_tr_l, TensorOrder2& c_tr_l_anti)
{
  double stl = 0;
  TComplex cu(0,1);
  for(unsigned int n_alpha = 0; n_alpha < 4; n_alpha++){
    for(unsigned int n_beta = 0; n_beta < 4; n_beta++){
      c_tr_l(n_alpha,n_beta) = 0.;
      stl = 0.;
      for(unsigned int n_sigma = 0; n_sigma < 4; n_sigma++){
	for(unsigned int n_delta = 0; n_delta < 4; n_delta++){
	  stl = stl - ((double)(*LeviCivitaOrder4())(n_alpha, n_beta, n_sigma, n_delta)
		       * xk(n_sigma)
		       * xkp(n_delta)
		       * (*Metric())(n_sigma, n_sigma)
		       * (*Metric())(n_delta, n_delta)
		       );
	}
      }
      c_tr_l(n_alpha,n_beta) = 8. * (xk(n_alpha) * xkp(n_beta)
                                     + xk(n_beta) * xkp(n_alpha) 
                                     - DotProdMetric(xk, xkp) * (*Metric())(n_alpha, n_beta)
                                     + cu * stl);

      c_tr_l_anti(n_alpha,n_beta) = 8. * (xk(n_alpha) * xkp(n_beta)
                                          + xk(n_beta) * xkp(n_alpha) 
                                          - DotProdMetric(xk, xkp) * (*Metric())(n_alpha, n_beta)
                                          - cu * stl);
    }
  }

  return;
}
//____________________________________________________________________________
double LARNCGammaXSec::DeltaPi(TensorOrder1& sp){
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,dimension(0:3) :: spp

  double temp = SMomentum(sp);

  double del = 1. / (temp - LARNCGammaXSec::fgMpi2);
  return del;
}
//____________________________________________________________________________
void LARNCGammaXSec::AEM(double t1, double t2, TensorOrder1& sp, TensorOrder1& sq, TensorOrder4& caem)
{
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,dimension(0:3) :: sp,sq,spd
  // real*8,dimension(3:6) :: fcat
  // real*8,dimension(3:5) :: fcv,fcvt
  // complex*16,dimension(4,4) :: cq,cem
  // complex*16,dimension(0:3,0:3,4,4) :: CAEM
  // real*8,external :: DotProdMetric
  TensorOrder1 fcv ("fcv", 3);
  TensorOrder1 fcvt("fcvt",3);
  TensorOrder1 fcat("fcat",4);
  TensorOrder2 cem ("cem", 4);
    
  TensorOrder1 spd("spd",4);
  spd= sp + sq;
  TensorOrder2* cq = SlashDirac(sq);
  //        t1=0.d0
  //        t2=DotProdMetric(xq,xq)
  LARNCGammaXSec::FactorC(t1, t2, fcv, fcvt, fcat);             
  
  for(int n1 = 0; n1 < 4; n1++){
    for(int n2 = 0; n2 < 4; n2++){
      TComplex cs2 = (fcv(4-3) / 
                      LARNCGammaXSec::fgMnSq *
                      ((*Metric())(n1,n2) * DotProdMetric(sq,spd) - sq(n1) * spd(n2))
                      + fcv(5-3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1,n2) * DotProdMetric(sq,sp)
                                                             - sq(n1) * sp(n2)));
      
      for(int i = 0; i < 4; i++){
	for(int j = 0; j < 4; j++){
	  
	  for(int j2 = 0; j2 < 4; j2++){
	    TComplex cs1 = fcv(3-3) / LARNCGammaXSec::fgMn * ((*Metric())(n1, n2) *
                                                              (*cq)(i, j2) - sq(n1) *
                                                              (*DiracMatrix(n2))(i, j2));
	    cem(i, j2) = cs1 + cs2 * (*UnityOrder2())(i, j2);
	  }
	  caem(n1, n2, i, j) = 0.;
	    
	  for(int k2 = 0; k2 < 4; k2++){
	    caem(n1, n2, i, j) = caem(n1, n2, i, j) + cem(i, k2) * (*DiracMatrix(5))(k2, j);
	  }
	}
	  
      }
    }
  }
  return;
}
//____________________________________________________________________________
void LARNCGammaXSec::ANC(double t1, double t2, TensorOrder1& sp, TensorOrder1& sq, TensorOrder4& canc)
{
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,dimension(0:3) :: sp,sq,spd
  // real*8,dimension(3:6) :: fcat
  // real*8,dimension(3:5) :: fcv,fcvt
  // complex*16,dimension(4,4) :: cq,cnc1,cnc2
  // complex*16,dimension(0:3,0:3,4,4) :: CANC
  // real*8,external:: DotProdMetric

  TensorOrder1 fcv ("fcv",  3);
  TensorOrder1 fcvt("fcvt", 3);
  TensorOrder1 fcat("fcat", 4);
  TensorOrder2 cnc1("cnc1", 4);
  TensorOrder2 cnc2("cnc2", 4);

  TensorOrder1 spd;

  spd = sp + sq;
  TensorOrder2* cq = SlashDirac(sq);
  //       t1=0.d0
  //       t2=DotProdMetric(sq,sq)
  LARNCGammaXSec::FactorC(t1,t2,fcv,fcvt,fcat);

  for(int n1 = 0; n1 < 4; n1++){
    for(int n2 = 0; n2 < 4; n2++){
      
      TComplex cs2 = (fcvt(4-3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1,n2) * DotProdMetric(sq,spd) - sq(n1) * spd(n2)) +
                      fcvt(5-3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1,n2) * DotProdMetric(sq,sp)  - sq(n1) * sp(n2)));
      TComplex ctemp1 = (fcat(4-3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1,n2) * DotProdMetric(sq,spd) - sq(n1) * spd(n2)) +
                         fcat(5-3) * (*Metric())(n1,n2) +
                         fcat(6-3) / LARNCGammaXSec::fgMnSq * sq(n1) * sq(n2));
      
      for(int i = 0; i < 4; i++){
	for(int j = 0; j < 4; j++){
	  
	  for(int j1 = 1; j1 < 4; j1++){
	    TComplex cs1 = fcvt(3-3) / LARNCGammaXSec::fgMn * ((*Metric())(n1,n2) * (*cq)(i,j1)
                                                               - sq(n1) * (*DiracMatrix(n2))(i,j1));
	    cnc1(i,j1) = cs1 + cs2 * (*UnityOrder2())(i,j1);
	  }
	  cnc2(i,j) = 0.;
	  for(int k = 0; k < 4; k++){
	    cnc2(i,j) = cnc2(i,j) + cnc1(i,k) * (*DiracMatrix(5))(k,j);
	  }
	  
	  canc(n1, n2, i, j) = cnc2(i, j) + ctemp1 * (*UnityOrder2())(i, j);
	}
      }
    }
  }  
  return;
}
//____________________________________________________________________________
TComplex LARNCGammaXSec::cDdelta(TensorOrder1& spp)
{
  //D_delta FUNCTION
  // Whatever that is..
  
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,dimension(0:3) :: spp
  // real*8,external :: Flam

  //   spp2=SMomentum(spp)
  double spp2 = spp(0)*spp(0) - spp(1)*spp(1) - spp(2)*spp(2) - spp(3)*spp(3);
  double wd;
  
  if(spp2 - TMath::Power(LARNCGammaXSec::fgMn + LARNCGammaXSec::fgMpi2, 2.) > 0.){
    double slam = Flam(spp2,TMath::Power(LARNCGammaXSec::fgMpi2, 2), LARNCGammaXSec::fgMnSq);
    wd = 1. / (6. * TMath::Pi()) * TMath::Power(LARNCGammaXSec::fgFStar / LARNCGammaXSec::fgMpi2,2) * LARNCGammaXSec::fgMn / TMath::Power(spp2, 2.)
      * TMath::Power(TMath::Sqrt(slam) / (2.), 3.); // *Ftheta(spsqr-(fgMn+LARNCGammaXSec::fgMpi2)**2);// csw-fgMn-LARNCGammaXSec::fgMpi2
  }else{
    wd = 0.;
  }
  TComplex cu(0, 1);
  TComplex cDdelta = 1./ (spp2 - TMath::Power(LARNCGammaXSec::fgMdelta, 2) + cu * wd * LARNCGammaXSec::fgMdelta);
  // 1.d0/(csw+LARNCGammaXSec::fgMdelta)/(csw-LARNCGammaXSec::fgMdelta+cu*cwd/2.d0)
  // ***************++hills
  //     IF ((spp2-(fgMn+LARNCGammaXSec::fgMpi2)**2).GT.0) then
  //     wd=0.12*( sqrt(Flam(spp2,LARNCGammaXSec::fgMpi2**2,LARNCGammaXSec::fgMnSq)) *LARNCGammaXSec::fgMdelta/(sqrt(spp2)*sqrt(Flam(LARNCGammaXSec::fgMdelta**2,LARNCGammaXSec::fgMpi2**2,fgMn2))))**3
  //    else 
  //    wd=0.
  //    endif

  //   cDdelta=1.d0/(spp2-LARNCGammaXSec::fgMdelta**2+cu*wd*LARNCGammaXSec::fgMdelta)
  return cDdelta;
}
//____________________________________________________________________________
void LARNCGammaXSec::FactorC(double t1, double t2, TensorOrder1& fcv, TensorOrder1& fcvt, TensorOrder1& fcat)
{
  //    the factor with ~t is for NC,  without ~t is for EM;
  if(LARNCGammaXSec::fgnFF == 1){
    double sma2 = TMath::Power(1.05, 2);
    double sfca50 = 1.2;
    double sfcv30 = 1.95;
    
    fcv(2) = sfcv30;   // (1.- t1/xmv**2)**2/(1.-t1/(4.*xmv**2))    // fcv30=1.95   xmv=0.84GeV;
    fcv(3) = -LARNCGammaXSec::fgMn / LARNCGammaXSec::fgMdelta * fcv(2);
    fcv(4) = 0;
    
    fcat(2) = 0;
    fcat(4) = sfca50 / TMath::Power(1 - t2 / sma2, 2) / (1 - t2 / (3 * sma2));    // xma=1.05 GeV, fca50=1.2;
    fcat(3) = -fcat(4) / 4.;
    fcat(5) =  fcat(4) * LARNCGammaXSec::fgMnSq / (TMath::Power(LARNCGammaXSec::fgMpi2, 2) - t2);
    
    fcvt(2) = (1 - 2 * LARNCGammaXSec::fgSinThetaW2) * sfcv30 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (4 * LARNCGammaXSec::fgMv2));
    fcvt(3) = -LARNCGammaXSec::fgMn / LARNCGammaXSec::fgMdelta * fcvt(2);
    fcvt(4) = 0;
    
  }else if(LARNCGammaXSec::fgnFF == 2){
    // Hills form factor
    double sma2 = TMath::Power(1.05, 2);
    double sfca50 = 1.2;
    
    fcv(2) =  2.13;  // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2))
    fcv(3) = -1.51;  // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2))
    fcv(4) =  0.48;  // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2))
    
    fcat(2) = 0;
    fcat(4) =  sfca50 / TMath::Power(1 - t2 / sma2, 2) / (1 - t2 / (3 * sma2));
    fcat(3) = -fcat(4) / 4.;
    fcat(5) =  fcat(4) * LARNCGammaXSec::fgMnSq / (TMath::Power(LARNCGammaXSec::fgMpi2, 2) - t2);
    
    fcvt(2) =  2.13 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (4     * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    fcvt(3) = -1.51 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (4     * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    fcvt(4) =  0.48 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (0.776 * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    
  }else if(LARNCGammaXSec::fgnFF == 3){
    double fca50 = LARNCGammaXSec::fgFca5P33;
    fcv(2) =  2.13; // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2));
    fcv(3) = -1.51; // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2));
    fcv(4) =  0.48; // (1.-t1/xmv**2)**2/(1.-t1/(4.*xmv**2));
    
    fcat(2) =  0;
    fcat(4) =  fca50 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMaSq, 2);
    fcat(3) = -fcat(4) / 4.;
    fcat(5) =  fcat(4) * LARNCGammaXSec::fgMnSq / (TMath::Power(LARNCGammaXSec::fgMpi2, 2) - t2);
    
    fcvt(2) =  2.13 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (4     * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    fcvt(3) = -1.51 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (4     * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    fcvt(4) =  0.48 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2) / (1 - t2 / (0.776 * LARNCGammaXSec::fgMv2)) * (1 - 2 * LARNCGammaXSec::fgSinThetaW2);
    
  }else if(LARNCGammaXSec::fgnFF == 4){
    double fca50 = LARNCGammaXSec::fgFca5P33;
    double q2 = 0, sc3v, sc4v, sc5v, sc3a, sc4a, sc5a, sc6a;
    LARNCGammaXSec::P33_1232FormFactor(q2, sc3v, sc4v, sc5v, sc3a, sc4a, sc5a, sc6a);
    fcv(2) = sc3v;
    fcv(3) = sc4v;
    fcv(4) = sc5v;
    double s3v, s4v, s5v, s3a, s4a, s5a, s6a;
    LARNCGammaXSec::P33_1232FormFactor(t2, s3v, s4v, s5v, s3a, s4a, s5a, s6a);
    fcvt(2) = (1 - 2 * LARNCGammaXSec::fgSinThetaW2) * s3v;
    fcvt(3) = (1 - 2 * LARNCGammaXSec::fgSinThetaW2) * s4v;
    fcvt(4) = (1 - 2 * LARNCGammaXSec::fgSinThetaW2) * s5v;
    
    fcat(2) =  0;
    fcat(4) =  fca50 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMaSq, 2);
    fcat(3) = -fcat(4) / 4.;
    fcat(5) =  fcat(4) * LARNCGammaXSec::fgMnSq / (TMath::Power(LARNCGammaXSec::fgMpi2, 2) - t2);
  }else if(LARNCGammaXSec::fgnFF == 5){
    // Hills form factor
    //double smv2 = TMath::Power(0.8, 2);
    double sma2 = 1;
    double sfca50 = 1.2;
    double sfcv30 = 2.0 ;
    
    fcv(2) = sfcv30;   // (1.- t1/xmv**2)**2/(1.-t1/(4.*xmv**2))    // fcv30=1.95   xmv=0.84GeV;
    fcv(3) = -LARNCGammaXSec::fgMn / LARNCGammaXSec::fgMdelta * fcv(2);
    fcv(4) = 0.;
    
    fcat(2) = 0.;
    fcat(4) = sfca50 / TMath::Power(1 - t2 / sma2, 2); // (1.- t2/(3.*fgMaSq))  ;
    fcat(3) = -fcat(4) / 4.;
    fcat(5) =  fcat(4) * LARNCGammaXSec::fgMnSq / (TMath::Power(LARNCGammaXSec::fgMpi2, 2) - t2);
    
    fcvt(2) = (1 - 2 * LARNCGammaXSec::fgSinThetaW2) * sfcv30 / TMath::Power(1 - t2 / LARNCGammaXSec::fgMv2, 2); ///(1.-t2/(4.*xmv**2));
    fcvt(3) = -LARNCGammaXSec::fgMn / LARNCGammaXSec::fgMdelta * fcvt(2);
    fcvt(4) = 0;
  }else{
    std::cout << "fgnFF should be 1 to 4" << std::endl;
    return;
  }
}
//____________________________________________________________________________
void LARNCGammaXSec::VertexAB(double t1, double t2, TensorOrder4& ch_vertex1, TensorOrder4& ch_vertex2)
{
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // complex*16,dimension(0:3,0:3,4,4) :: ch_vertex1,ch_vertex2
  // complex*16,dimension(0:3,4,4) :: cA_nc,cA_em,cA_nc_con,cA_em_con cl_vertex
  // complex*16,dimension(4,4) :: cA_nc2,cA_em2,cA_nc_con2,cA_em_con2
  // complex*16,dimension(4,4) :: cpqm,cppqm,c_ga_n1,c_ga_n2,c_ga_con1,c_ga_con2,c_ga0,c_p0,c_q0
  double f1, f2, ft1, ft2, fta;
  LARNCGammaXSec::ElectroMagneticFormFactor(t1, fNucleon, f1,  f2);
  LARNCGammaXSec::NeutralCurrentFormFactor (t2, fNucleon, ft1, ft2, fta);
  TensorOrder3 cA_nc;
  TensorOrder3 cA_em;
  TensorOrder3 cA_nc_con;
  TensorOrder3 cA_em_con;
  TComplex cu(0,1);
  
  for(int n1 = 0; n1 < 4; n1++){
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){
	TComplex c_sigma_q = 0.;
	TComplex c_sigma_qf = 0.;
	for(int n2 = 0; n2 < 4; n2++){
	  c_sigma_q  = c_sigma_q +  (*DiracTensor(n1, n2))(i, j) * xq(n2)  * (*Metric())(n2, n2);
	  c_sigma_qf = c_sigma_qf + (*DiracTensor(n1, n2))(i, j) * xqf(n2) * (*Metric())(n2, n2);
	}
	
	double c_ga_5 = 0.;
	for(int k = 0; k < 4; k++){
	  c_ga_5 = c_ga_5 + (*DiracMatrix(n1))(i, k) * (*DiracMatrix(5))(k, j);
	}
	
	// the Fig A and Fig B
	cA_nc(n1, i, j) = ft1 * (*DiracMatrix(n1))(i, j) + ft2 * cu * c_sigma_q  / (2. * LARNCGammaXSec::fgMn) + c_ga_5 * fta; // For the vertex of hadron
	cA_em(n1, i, j) = f1  * (*DiracMatrix(n1))(i, j) - f2  * cu * c_sigma_qf / (2. * LARNCGammaXSec::fgMn);
	cA_nc_con(n1, i, j) = ft1 * (*DiracMatrix(n1))(i, j) - ft2 * cu * c_sigma_q  / (2. * LARNCGammaXSec::fgMn) + c_ga_5 * fta; // for the conj vertex           
	cA_em_con(n1, i, j) = f1  * (*DiracMatrix(n1))(i, j) + f2  * cu * c_sigma_qf / (2. * LARNCGammaXSec::fgMn);
      }
    }
  }
  TensorOrder2 cpqm ;
  TensorOrder2 cppqm;
  cpqm  = c_p  + c_q + LARNCGammaXSec::fgMn * (*UnityOrder2());
  cppqm = c_pp - c_q + LARNCGammaXSec::fgMn * (*UnityOrder2());
  
  TensorOrder1 sum1;
  TensorOrder1 sum2;
  sum1 = xp  + xq;
  sum2 = xpp - xq;
  double dpq  = LARNCGammaXSec::DeltaPropagator(sum1);
  double dppq = LARNCGammaXSec::DeltaPropagator(sum2);
  
  //  if(fNucleon==1.and.w_sqf0.lt.0.1) then
  //        dpq=0.
  //        dppq=0.
  //  endif

  for(int n_miu = 0; n_miu < 4; n_miu++){
    for(int n_alpha = 0; n_alpha < 4; n_alpha++){
      TensorOrder2* cA_em2 = Dim3to2(cA_em, n_miu);
      TensorOrder2* cA_nc2 = Dim3to2(cA_nc, n_alpha);
      TensorOrder2* cA_em_con2 = Dim3to2(cA_em_con, n_miu);
      TensorOrder2* cA_nc_con2 = Dim3to2(cA_nc_con, n_alpha);
      TensorOrder2* c_ga_n1 = Mult3Matrices((*cA_em2), cpqm,  (*cA_nc2));
      TensorOrder2* c_ga_n2 = Mult3Matrices((*cA_nc2), cppqm, (*cA_em2));
      TensorOrder2* c_ga_con1 = Mult3Matrices((*cA_nc_con2), cpqm,  (*cA_em_con2));
      TensorOrder2* c_ga_con2 = Mult3Matrices((*cA_em_con2), cppqm, (*cA_nc_con2));
      for(int i = 0; i < 4; i++){
	for(int j= 0; j < 4; j++){
	  ch_vertex1(n_miu, n_alpha, i, j) = (*c_ga_n1)(i,j)   * dpq + (*c_ga_n2)(i,j)   * dppq;
	  ch_vertex2(n_miu, n_alpha, i, j) = (*c_ga_con1)(i,j) * dpq + (*c_ga_con2)(i,j) * dppq;
	}
      }
    }
  }

}

//____________________________________________________________________________
void LARNCGammaXSec::VertexCD(double t1, double t2, TensorOrder4& ch_ver_cd, TensorOrder4& ch_ver_cdt)
{
  // Gamma_delta for the Fia C and Fig D
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // complex*16,dimension(0:3,0:3,4,4) :: ch_ver_cd,ch_ver_cdt,cAC_emt,cAC_nc,        &
  //      cAD_em,cAD_nct,cAC_em,cAD_nc
  // complex*16,dimension(4,4) :: clpq,clppq,cg_c,cg_d,cgtem_c,cgtem_d,               &
  //      cacemt,cacnc,cadnct,cadem
  
  TComplex cdeltapq  =  LARNCGammaXSec::cDdelta(xpd);  // xpd=xp+xq
  TComplex cdeltappq =  LARNCGammaXSec::cDdelta(xpdc); // xpdc=xpp-xq
  //  in_medium affection of Delta 

  // call medium(xpd,densi,cdeltapq)
  // call medium(xpdc,densi,cdeltappq)
  
  TensorOrder4 cAC_em;
  TensorOrder4 cAD_em;
  TensorOrder4 cAD_nc;
  TensorOrder4 cAC_nc;

  AEM(t1, t2, xpp, xqf, cAC_em);
  TensorOrder4 cAC_emt(cAC_em);
  cAC_emt.Conjugate();

  ANC(t1, t2, xp,   xq, cAC_nc);
  TensorOrder1 minusxq;
  minusxq = -xq;
  TensorOrder1 minusxqf;
  minusxqf = -xqf;

  ANC(t1, t2, xpp, minusxq, cAD_nc);

  TensorOrder4 cAD_nct(cAD_nc);
  cAC_emt.Conjugate();
    
  AEM(t1, t2, xp, minusxqf, cAD_em);

  for(int nmiu = 0; nmiu < 4; nmiu++){
    for(int nalpha = 0; nalpha < 4; nalpha++){
      
      TensorOrder2 cgtem_c;
      TensorOrder2 cgtem_d;
      
      
      for(int ndelta = 0; ndelta < 4; ndelta++){
	for(int nsigma = 0; nsigma < 4; nsigma++){
	  
	  // For Fig C    
	  TensorOrder2* cacemt = Dim4to2(cAC_emt,ndelta,nmiu);
	  TensorOrder2* cacnc  = Dim4to2(cAC_nc, nsigma,nalpha);
	  TensorOrder2* clpq = Lambda(ndelta, nsigma, xpd);
	  TensorOrder2* cg_c = Mult3Matrices((*cacemt), (*clpq), (*cacnc));
              
	  cgtem_c = cgtem_c  + (*cg_c) * (*Metric())(ndelta,ndelta)
	    * (*Metric())(nsigma,nsigma);

	  // For Fig D
	  TensorOrder2* cadnct = Dim4to2(cAD_nct,ndelta,nalpha);
	  TensorOrder2* cadem  = Dim4to2(cAD_em,nsigma,nmiu);
	  TensorOrder2* clppq = Lambda(ndelta, nsigma, xpdc);
	  TensorOrder2* cg_d  = Mult3Matrices((*cadnct), (*clppq), (*cadem));
		
	  cgtem_d = cgtem_d  + (*cg_d) * (*Metric())(ndelta,ndelta)
	    * (*Metric())(nsigma,nsigma);
	}
      }
      if(fNDiag == 2 && fMDelta == 1)
	cdeltappq = 0.;
      for(int i = 0; i < 4; i++){
	for(int j = 0; j < 4; j++){
	  ch_ver_cd(nmiu,nalpha,i,j) = cgtem_c(i,j) * cdeltapq + cgtem_d(i,j) * cdeltappq;
	}
      }
    }
  }
  ch_ver_cdt = ch_ver_cd;
  ch_ver_cdt.Conjugate();
  return;
}
//____________________________________________________________________________
void LARNCGammaXSec::VertexE(TensorOrder4& ch_ver_e)
{
  // use parameter
  //   implicit real*8 (a,b,d-h,o-z)
  //   implicit complex*16 (c)
  //   complex*16,dimension(0:3,0:3,4,4) :: ch_ver_e
  //   complex*16,dimension(4,4) :: csp5,csppp
  //   real*8,dimension(0:3,0:3) :: st
  //   integer,external :: Lvten
  TComplex cu(0,1);
  TensorOrder1 subtraction = xpp - xp;
  double sdel = LARNCGammaXSec::DeltaPi(subtraction);
  TComplex csfactor = - cu * (double)fNucleon * (-LARNCGammaXSec::fgGa) / (8. * TMath::Power(TMath::Pi(), 2.) * TMath::Power(LARNCGammaXSec::fgFPi, 2.)) * (0.5 - 2. * LARNCGammaXSec::fgSinThetaW2) * sdel; // 4.*sin_theta_w2   (1.-0.23122);
  TensorOrder2 csppp = c_pp - c_p;
  TensorOrder2* csp5  = MatMult(csppp, (*DiracMatrix(5)));
  TensorOrder2 st;

  for(int nmiu = 0; nmiu < 4; nmiu++){
    for(int nalpha = 0; nalpha < 4; nalpha++){
      for(int nsigma = 0; nsigma < 4; nsigma++){
	for(int ndelta = 0; ndelta < 4; ndelta++){
	  st(nmiu, nalpha) = st(nmiu, nalpha) +
	    (double)(*LeviCivitaOrder4())(ndelta, nsigma, nmiu, nalpha) *
	    xqf(ndelta) * xq(nsigma);
	}
      }
      for(int i = 1; i < 4; i++){
	for(int j = 1; j < 4; j++){
	  ch_ver_e(nmiu,nalpha,i,j) = st(nmiu, nalpha) * (*csp5)(i, j);
	}
      }
    }
  }
  
  ch_ver_e = ch_ver_e * csfactor;
  //    ch_ver_et=-ch_ver_e  

  return;
}
//____________________________________________________________________________
void LARNCGammaXSec::VertexJ12(double t1, TensorOrder4& ch_verj12, TensorOrder4& ch_verj12_t)
{
  // fNParity=1 for P11(1440), fNParity=-1 for S11(1535);
  TensorOrder2 cpqm ;
  TensorOrder2 cppqm;
  cpqm  = c_p  + c_q + LARNCGammaXSec::fgMn * (*UnityOrder2());
  cppqm = c_pp - c_q + LARNCGammaXSec::fgMn * (*UnityOrder2());
  int nexcit;
  // double fem, fvem, fvnc, fanc;
  TensorOrder2 fem ("fem",  2);
  TensorOrder2 fvem("fvem", 2);
  TensorOrder1 fanc("fanc", 2);
  TensorOrder2 fvnc("fvnc", 2);
    

  if(fNParity == 1){
    LARNCGammaXSec::P11_1440FormFactor(t1, fem, fvem, fvnc, fanc);
    nexcit = 2;
  }else if(fNParity == -1){
    LARNCGammaXSec::S11_1535FormFactor(t1, fem, fvem, fvnc, fanc);
    nexcit = 4;
  }else{
    std::cout << "fNParity should be 1 or -1. Here is Form Factor of J=1/2" << std::endl;
    return;
  }
  
  int ii = -(fNucleon - 3) / 2;       // For neutron, fNucleon = -1, i=2, for proton, fNucleon=1, i=1;
  double f1 = fem(ii,0);
  double f2 = fem(ii,1);
  double fa = 0;
  double ft1 = fvnc(ii,0);
  double ft2 = fvnc(ii,1);
  double fta = fanc(ii);
  TensorOrder1 minus;
  minus= -xqf;
  TensorOrder3 cA_em;
  TensorOrder3 cA_nc;
  Vertex12(f1,  f2,  fa,  minus,  cA_em);
  Vertex12(ft1, ft2, fta, xq,     cA_nc);
  
  TComplex cdpq  = LARNCGammaXSec::Propagator(nexcit, xpd);
  TComplex cdppq = LARNCGammaXSec::Propagator(nexcit, xpdc);
  
  for(int n_miu = 0; n_miu < 3; n_miu++){
    for(int n_alpha = 0; n_alpha < 3; n_alpha++){

       TensorOrder2* cA_em2 = Dim3to2(cA_em, n_miu);
       TensorOrder2* cA_nc2 = Dim3to2(cA_nc, n_alpha);
       TensorOrder2* c_ga_n1 = Mult3Matrices((*cA_em2), cpqm,  (*cA_nc2));
       TensorOrder2* c_ga_n2 = Mult3Matrices((*cA_nc2), cppqm, (*cA_em2));

      for(int i = 0; i < 5; i++){
	for(int j = 0; j < 5; j++){
	  ch_verj12(n_miu,n_alpha,i,j) = (*c_ga_n1)(i,j) * cdpq + (*c_ga_n2)(i,j) * cdppq;
	}
      }
    }
  }
  ch_verj12_t = ch_verj12;
  ch_verj12_t.Conjugate();
}
//____________________________________________________________________________
void LARNCGammaXSec::Vertex12(double f1, double f2, double fa, TensorOrder1& sq, TensorOrder3& ver12)
{
  TensorOrder2* csq = SlashDirac(sq);
  TensorOrder3 ver12p;
  TComplex cu(0,1);
  for(int mu = 0; mu < 3; mu++){
    for(int i = 1; i < 4; i++){
      for(int j = 1; j < 4; j++){
	double c_sigma_q = 0.;
	for(int nu = 0; nu < 3; nu++){
	  c_sigma_q = c_sigma_q + (*DiracTensor(mu, nu))(i, j) * sq(nu) * (*Metric())(nu, nu);      // sigma(mu,nu)*q(nu)
	}
	
	double c_ga_5 = 0.;
	for(int k = 1; k < 4; k++){
	  c_ga_5 = c_ga_5 + (*DiracMatrix(mu))(i, k) * (*DiracMatrix(5))(k, j);                         // ga(mu)*ga5
	}
	
	ver12p(mu, i, j) = f1 * ((*csq)(i,j) * sq(mu) - DotProdMetric(sq, sq) * (*DiracMatrix(mu))(i, j)) / (4 * LARNCGammaXSec::fgMnSq) + f2 * cu * c_sigma_q / (2. * LARNCGammaXSec::fgMn) + c_ga_5 * fa;
      }
      
      for(int j2 = 1; j2 < 4; j2++){
	if(fNParity == 1){
	  ver12(mu, i, j2) = ver12p(mu, i, j2);     // For the positive parity P11(1440)
	}else if(fNParity == -1){
	  ver12(mu, i,j2) = 0.;
          for(int k2 = 1; k2 < 4; k2++){
	    ver12(mu, i, j2) = ver12(mu, i, j2) + ver12p(mu, i, k2) * (*DiracMatrix(5))(k2, j2);        // For the negative parity S11(1535)
	  }
	}
      }
    }
  }
}
//____________________________________________________________________________
void LARNCGammaXSec::Vertex32(TensorOrder1& fcv, TensorOrder1& fca, TensorOrder1& sp, TensorOrder1& sq, TensorOrder4& cver32)
{
  //V-A(3/2) vertex32
  fcv(5) = 0;

  TensorOrder2* csq = SlashDirac(sq);
  TensorOrder1 spq;
  spq = sp + sq;
  TensorOrder2 cvector0;
  TensorOrder2 caxial0;
  TensorOrder4 cvector;
  TensorOrder4 caxial;

  //TensorOrder2& caxial0 = TensorOrder2(4);
  for(int n1 = 0; n1 < 4; n1++){
    for(int n2 = 0; n2 < 4; n2++){
      TComplex cvector1 = fcv(3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1, n2) * DotProdMetric(sq, spq) - sq(n1) * spq(n2)) + fcv(5) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1, n2) * DotProdMetric(sq, sp) - sq(n1) * sp(n2));//  +fcv(6)*(*Metric())(n1,n2);//    the vector part
      TComplex caxial1  = fca(3) / LARNCGammaXSec::fgMnSq * ((*Metric())(n1, n2) * DotProdMetric(sq, spq) - sq(n1) * spq(n2)) + fca(5) * (*Metric())(n1, n2) + fca(6) / LARNCGammaXSec::fgMnSq * sq(n1) * sq(n2); //         &    the minus axial part
      for(int i = 1; i < 5; i++){
	for(int j2 = 1; j2 < 5; j2++){
	  TComplex cvector2 = fcv(3) / LARNCGammaXSec::fgMn * ((*Metric())(n1,n2) * (*csq)(i,j2) -
                                                               sq(n1)*(*DiracMatrix(n2))(i,j2));
	  cvector0(i,j2) = cvector2 + cvector1 * (*UnityOrder2())(i,j2);
	  // caxial2=fca(3)/LARNCGammaXSec::fgMn*((*Metric())(n1,n2)*csq(i,j2)-sq(n1)*DiracMatrix()(n2,i,j2)    fca(3)=0
          // caxial0(i,j2)=caxial2 + caxial1*UnityOrder2()(i,j2) 
	  caxial0(i,j2) = caxial1 * (*UnityOrder2())(i,j2);
	}
	
	if(fNParity == 1){ // the parity of P33(1232) is positive
	  for(int j = 1; j < 5; j++){ 
	    cvector(n1,n2,i,j)=0.;
	    for(int k2 = 1; k2< 5; k2++){
	      cvector(n1,n2,i,j) = cvector(n1,n2,i,j)+ cvector0(i,k2) * (*DiracMatrix(5))(k2,j);
	    }
	    cver32(n1,n2,i,j) = cvector(n1,n2,i,j) + caxial0(i,j);
	  }
	}else if(fNParity==-1){ // the parity of D13(1520) is negative
	  for(int j = 1; j < 5; j++){
	    caxial(n1,n2,i,j)=0.;
	    for(int k2 = 1; k2 < 5; k2++){
	      caxial(n1,n2,i,j) = caxial(n1,n2,i,j) + caxial0(i,k2) * (*DiracMatrix(5))(k2,j);
	    }
	    cver32(n1,n2,i,j)= caxial(n1,n2,i,j) + cvector0(i,j);
	  }
	}else{
	  std::cout << "parity of the P33(1232) or D13(1520) is not correct" <<std::endl;
	  return;
	}
      }
    }
  }
}
//____________________________________________________________________________
void LARNCGammaXSec::VertexJ32(double t1, TensorOrder4& ch_verj32, TensorOrder4& ch_verj32_t)
{
  //For D13(1520)
  // fNParity=1 for P33(1232), fNParity=-1 for D13(1520)
  int nexcit;
  // double fem, fvem, fvnc, fanc;
  TensorOrder2 fem ("fem",  2, 5);
  TensorOrder2 fvem("fvem", 2);
  TensorOrder2 fanc("fanc", 2, 5);
  TensorOrder2 fvnc("fvnc", 2, 5);
    
  if(fNParity == 1){
    LARNCGammaXSec::P33_1232FormFactor(t1,fem,fvem,fvnc,fanc);
    nexcit=1;
  }else if(fNParity==-1){
    LARNCGammaXSec::D13_1520FormFactor(t1,fem,fvem,fvnc,fanc);
    nexcit=3;
  }else{ 
    std::cout << "fNParity should be 1 or -1. Here is Form Factor of J=3/2" << std::endl;
    return;
  }
  int ii;
  if(nexcit == 1)
    ii = 1;
  else if(nexcit==3)
    ii = -(fNucleon - 3) / 2; //        For neuturon, fNucleon=-1, i=2, for proton, fNucleon=1, i=1
  else{
    std::cout << "nexcit should be 1 or 3. Here is ver_j32" << std::endl;
    return;
  }

  TensorOrder1 fcv  ("fcv",  5);
  TensorOrder1 fcvt ("fcvt", 5);
  TensorOrder1 fcat ("fcat", 5);
  TensorOrder1 fca  ("fca",  5);

  for(int j = 0; j < 5; j++){
    fcv(j)  = fem(ii,  j);
    fcvt(j) = fvnc(ii, j);
    fcat(j) = fanc(ii, j);
    fca(j)  = 0.;
  }

  // For the direct diagram
  TensorOrder1 minus = -xqf;
  TensorOrder4 cAem;
  TensorOrder4 cAem_0;
  Vertex32(fcv, fca, xp,  minus,          cAem);
  Vertex32(fcv, fca, xpp, xqf, cAem_0);
  TensorOrder4 cAem_t = cAem_0;
  cAem_0.Conjugate();
  // For the crossed diagram
  TensorOrder2 minus2;
  minus2= -xq;
  TensorOrder4 cAnc;
  TensorOrder4 cAnc_0;
  Vertex32(fcvt, fcat, xp,  xq, cAnc);
  Vertex32(fcvt, fcat, xpp, minus,         cAnc_0);
  TensorOrder4 cAnc_t = cAnc_0;
  cAnc_t.Conjugate();
  
  TComplex cdpq  = LARNCGammaXSec::Propagator(nexcit, xpd);
  TComplex cdppq = LARNCGammaXSec::Propagator(nexcit, xpdc);

  for(int nmiu = 0; nmiu < 4; nmiu++){
    for(int nalpha = 0; nalpha < 4; nalpha++){
      TensorOrder2 cgtem_d; //  For the direct diagram
      TensorOrder2 cgtem_c; //  For the crossed diagram
      
      for(int ndelta = 0; ndelta < 4; ndelta++){
	for(int nsigma = 0; nsigma < 4; nsigma++){
	  // For the direct diagram
	  TensorOrder2* cA_emt = Dim4to2(cAem_t, ndelta, nmiu);
	  TensorOrder2* cA_nc  = Dim4to2(cAnc,   nsigma, nalpha);
          TensorOrder2* clpq = Lambda(ndelta,nsigma,xpd);
	  TensorOrder2* cg_d = Mult3Matrices((*cA_emt),(*clpq),(*cA_nc));

	  cgtem_d = cgtem_d + (*cg_d) * (*Metric())(ndelta,ndelta) * 
	    (*Metric())(nsigma,nsigma);

	  // For the crossed diagram 
	  TensorOrder2* cA_nct = Dim4to2(cAnc_t,ndelta,nalpha);
	  TensorOrder2* cA_em  = Dim4to2(cAem,nsigma,nmiu);
	  TensorOrder2* clppq = Lambda(ndelta,nsigma, xpdc);
	  TensorOrder2* cg_c = Mult3Matrices((*cA_nct),(*clppq),(*cA_em));

	  cgtem_c = cgtem_c + (*cg_c) * (*Metric())(ndelta,ndelta) *
	    (*Metric())(nsigma,nsigma);
	}
      }
      for(int i = 1; i < 5; i++){
	for(int j = 1; j < 5; j++){
	  ch_verj32(nmiu,nalpha,i,j) = cgtem_d(i,j) * cdpq + cgtem_c(i,j) * cdppq;
	}
      }
    }
  }

  ch_verj32_t = TensorOrder4(ch_verj32);
  ch_verj32_t.Conjugate();
}
//____________________________________________________________________________
void LARNCGammaXSec::AmpNum(double t1, double t2, TensorOrder2& c_lh_both)
{
  c_p  = (*SlashDirac(xp));
  c_pp = (*SlashDirac(xpp));
  c_k  = (*SlashDirac(xk));
  c_kp = (*SlashDirac(xkp));
  c_q  = (*SlashDirac(xq));
  c_qf = (*SlashDirac(xqf));

  TensorOrder2 TheIdMatrix = (*UnityOrder2());
  
  cpm  = c_p  + fgMn*TheIdMatrix;
  cppm = c_pp + fgMn*TheIdMatrix;

  TensorOrder4 ch_ver     ;
  TensorOrder4 ch_ver_t   ;
  TensorOrder4 ch_ver_ab  ;
  TensorOrder4 ch_ver_abt ;
  TensorOrder4 ch_ver_cd  ;
  TensorOrder4 ch_ver_cdt ;
  TensorOrder4 ch_ver_e   ;
  TensorOrder4 ch_ver_et  ;
  TensorOrder4 ch_ver_d13 ;
  TensorOrder4 ch_ver_d13t;
  TensorOrder4 ch_ver_p11 ;
  TensorOrder4 ch_ver_p11t;
  TensorOrder4 ch_ver_s11 ;
  TensorOrder4 ch_ver_s11t;

  if(fNDiag == 0){
    LARNCGammaXSec::VertexAB(t1, t2, ch_ver_ab, ch_ver_abt);
    LARNCGammaXSec::VertexCD(t1, t2, ch_ver_cd, ch_ver_cdt); 
    LARNCGammaXSec::VertexE(ch_ver_e);
    ch_ver_et = -ch_ver_et;
    fNParity = -1;
    LARNCGammaXSec::VertexJ32(t2, ch_ver_d13, ch_ver_d13t);
    fNParity = 1;
    LARNCGammaXSec::VertexJ12(t2, ch_ver_p11, ch_ver_p11t);
    fNParity = -1;
    LARNCGammaXSec::VertexJ12(t2, ch_ver_s11, ch_ver_s11t);

    ch_ver   = ch_ver_cd  + ch_ver_e  + ch_ver_ab  + ch_ver_s11  + ch_ver_p11  + ch_ver_d13;
    ch_ver_t = ch_ver_cdt + ch_ver_et + ch_ver_abt + ch_ver_s11t + ch_ver_p11t + ch_ver_d13t;
  }else if(fNDiag == 1){
    LARNCGammaXSec::VertexAB(t1, t2, ch_ver, ch_ver_t);
  }else if(fNDiag == 2){
    LARNCGammaXSec::VertexCD(t1, t2, ch_ver, ch_ver_t);
    fNParity = 1;
  }else if(fNDiag == 3){
    LARNCGammaXSec::VertexE(ch_ver);
    ch_ver_t = -ch_ver;
  }else if(fNDiag == 4){
    fNParity = -1;
    LARNCGammaXSec::VertexJ32(t2, ch_ver, ch_ver_t);
  }else if(fNDiag == 5){
    fNParity = 1;
    LARNCGammaXSec::VertexJ12(t2, ch_ver, ch_ver_t);
  }else if(fNDiag == 6){
    fNParity=-1;
    LARNCGammaXSec::VertexJ12(t2, ch_ver, ch_ver_t);
  }else if(fNDiag == 7){
    LARNCGammaXSec::VertexAB(t1, t2, ch_ver_ab, ch_ver_abt);
    LARNCGammaXSec::VertexCD(t1, t2, ch_ver_cd, ch_ver_cdt);
    LARNCGammaXSec::VertexE(ch_ver_e);
    ch_ver_et = -ch_ver_e;
    ch_ver    =  ch_ver_cd  + ch_ver_ab  + ch_ver_e;// + ch_ver_s11 +ch_ver_p11 + ch_ver_d13;
    ch_ver_t  =  ch_ver_cdt + ch_ver_abt + ch_ver_et;
  }else if(fNDiag == 8){
    LARNCGammaXSec::VertexAB(t1, t2, ch_ver_ab, ch_ver_abt);
    LARNCGammaXSec::VertexCD(t1, t2, ch_ver_cd, ch_ver_cdt); 
    ch_ver = ch_ver_cd  + ch_ver_ab;
    ch_ver_t = ch_ver_cdt + ch_ver_abt;
  }

  TensorOrder2 c_lh       ;
  TensorOrder2 c_lh_anti  ;
  TensorOrder2 c_tr_l     ;
  TensorOrder2 c_tr_l_anti;
  LARNCGammaXSec::TraceLight(c_tr_l, c_tr_l_anti);
  for(int n_alpha = 0; n_alpha < 4; n_alpha++){
    for(int n_beta = 0; n_beta < 4; n_beta++){
      double c_tr_h = Mult2(cpm, ch_ver_t, cppm, ch_ver, n_alpha, n_beta);
      c_lh      = c_lh      + c_tr_l     (n_alpha,n_beta) * c_tr_h * (*Metric())(n_alpha,n_alpha) * (*Metric())(n_beta,n_beta);
      c_lh_anti = c_lh_anti + c_tr_l_anti(n_alpha,n_beta) * c_tr_h * (*Metric())(n_alpha,n_alpha) * (*Metric())(n_beta,n_beta);
    }
  }
  // real part is for neutrino, image part is for antineutrino
  if(fNCheck == 520){
    if((Abs(Imaginary(c_lh)))     > TMath::Power(10., -6.) ||
       (Abs(Imaginary(c_lh_anti)) > TMath::Power(10., -6.))){
      LOG("LARNCGammaXSec", pFATAL) << "The image part of c_lh or c_lh_anti is not zero";
      LOG("LARNCGammaXSec", pFATAL) << "c_lh      = ";
      LOG("LARNCGammaXSec", pFATAL) << "c_lh_anti = ";
      exit(1);
    }
  }

  TComplex complexUnity(0,1);

  c_lh_both = Real(c_lh) + complexUnity * Real(c_lh_anti);

  if(fNCheck == 520){
    if((Real(c_lh_both)) < 0. ||
       (Imaginary(c_lh_both)) < 0.){
      LOG("LARNCGammaXSec", pFATAL) << "the LH tensor is less than zero";
      LOG("LARNCGammaXSec", pFATAL) << "c_lh_both      = ";
      exit(1);
    }
  }
  return;
  
}
double LARNCGammaXSec::SMomentum(TensorOrder1& sp){
  
  double smomentum=0.;

  for(int i = 0; i < 4; i++){
    smomentum=smomentum + sp(i)*sp(i)*(*Metric())(i,i);  
  }

  return smomentum;
}


double LARNCGammaXSec::Flam(double sx, double sy, double sz){
  return sx*sx + sy*sy + sz*sz -2.*(sx*sy + sy*sz + sx*sz);
}


TensorOrder2* LARNCGammaXSec::Dim3to2(TensorOrder3& tensor, int n){

  TensorOrder2* mat = new TensorOrder2();
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      (*mat)(i,j) = tensor(n, i, j);
  return mat;
}

TensorOrder2* LARNCGammaXSec::Dim4to2(TensorOrder4& tensor, int n1, int n2){

  TensorOrder2* mat = new TensorOrder2();
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      (*mat)(i,j) = tensor(n1, n2, i, j);
  return mat;
}


TensorOrder2* LARNCGammaXSec::Mult3Matrices(TensorOrder2& mat1,
                                            TensorOrder2& mat2,
                                            TensorOrder2& mat3){
  TensorOrder2* result;
  TensorOrder2* mult = MatMult(mat1, mat2);
  result = MatMult((*mult), mat3);
  return result;
}

TensorOrder2* LARNCGammaXSec::MatMult(TensorOrder2& mat1,
                                      TensorOrder2& mat2){

  if(mat1.GetOrderDim(1) != mat2.GetOrderDim(0))
    LOG("LARNCGammaXSec",kError)
      << "Matrices with entered orders can't be multiplied with each other.";
  
  TensorOrder2* multiply = new TensorOrder2(Form("(%s x %s)", mat1.GetName(), mat2.GetName()),
                                           mat1.GetOrderDim(0), mat2.GetOrderDim(1));
  TComplex sum (0.,0.);
  
  for (unsigned int c = 0; c < mat1.GetOrderDim(0); c++) {
    for (unsigned int d = 0; d < mat2.GetOrderDim(1); d++) {
      for (unsigned int k = 0; k < mat1.GetOrderDim(1); k++) {
        sum = sum + mat1(c,k) * mat2(k,d);
      }
      
      (*multiply)(c,d) = sum;
      sum = TComplex(0.,0.);
    }
  }

  return multiply;
  
}

TComplex LARNCGammaXSec::Mult2(TensorOrder2& c_a, TensorOrder4& c_b,
                               TensorOrder2& c_c, TensorOrder4& c_d,
                               int n_alpha,int n_beta){
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // parameter(n=4)                                   ! n is the dimension of the matrixs
  // complex*16,dimension(n,n) :: c_a,c_c,c_m1,c_m2,c_m 
  // complex*16,dimension(0:3,0:3,n,n) :: c_b,c_d

  TComplex ch(0.,0.);
  TensorOrder2 c_m1;
  TensorOrder2 c_m2;
  TensorOrder2* c_m;
  
  for(int n_miu=0; n_miu<4; n_miu++){
    for(int n_niu=0; n_niu<4; n_niu++){
      
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          c_m1(i,j)=0.;
          c_m2(i,j)=0.;
          for(int k=0; k<4; k++){
            c_m1(i,j)=c_m1(i,j) + c_a(i,k)*c_b(n_miu,n_alpha,k,j);
            c_m2(i,j)=c_m2(i,j) + c_c(i,k)*c_d(n_niu,n_beta,k,j);
          }
        }
      }
      c_m = MatMult(c_m1,c_m2);
      
      TComplex ctr(0.,0.);
      for(int i=0; i<4; i++){
        ctr=ctr + (*c_m)(i,i);
      }

      ch=ch+ctr*(*UnityOrder2(4))(n_miu,n_niu);
    }
  }
  delete c_m;
  return -ch/2.;
}

TensorOrder2* LARNCGammaXSec::Lambda(int ns, int nd, TensorOrder1& ppd){
  //   use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // complex*16,dimension(4,4) :: cppd,clam,cpmd,csbar
  // real*8,dimension(0:3) :: ppd
  
  TensorOrder2* cppd =  DiracSlach(ppd);
  TensorOrder2 cpmd = (*cppd)+xmd*unm;
  TensorOrder2 csbar;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      TComplex csdij(0.,0.);
      for(int k=0; k<4; k++){
        csdij=csdij+DM[ns](i,k)*DM[nd](k,j);
      }
      
      csbar(i,j) = metric(ns,nd)*unity2(i,j)-2./3.*ppd(ns)*ppd(nd)/xmd*xmd*unity2(i,j) 
        +(ppd(ns)*DM[nd](i,j)-ppd(nd)*DM[ns](i,j))/(3.*xmd)
        -csdij/3.;
    }
  }
  TensorOrder2* = new TensorOrder2();
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      (*clam)(i,j)=0.;
      for(int k=0; k<4; k++){
        (*clam)(i,j)=(*clam)(i,j)-cpmd(i,k)*csbar(k,j);
      }
    }
  }
  return clam;
  
}

TComplex LARNCGammaXSec::Width(Resonance_t res, double sp2){
  // use parameter
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,dimension(4,5) :: widd
  // real*8,external :: flam
  // ! xmsigma=0.475           ! For Sigma-> Pi Pi in S-wave
  // ! xmeta=0.548

  TComplex wid;
  switch(res){
  case kNoResonance:
    wid=0.;
    break;
  case kP33_1232:   
    if((sp2 - (fgMnucl+fgMpi)*(fgMnucl+fgMpi)) > 0){
      double slam = flam(sp2, fgMpi*fgMpi, fgMnSq);
      wd= 1. / (6.*TMath::Pi())*(fgFStar/fgMpi)*(fgFStar/fgMpi)*fgMnucl/(sp2*sp2) * TMath::Power(TMath::Sqrt(slam)/(2.),3.);
    }else{
      wd=0.;
    }
    double smr=xmd;
    double pwidnpi0 = 0.117;
    pwidth1(sp2, smr,fgMnucl,fgMpi,1,pwidnpi0,pwidnpi);                           //decay mode N-Pi
    wid=wd;
    // P11(1440)
    break;
  case kP11_1440:
    smr=xmp11;
    pwidnpi0=0.3*0.65;  
    pwidth1(sp2,smr,fgMnucl,fgMpi,1,pwidnpi0,pwidnpi);//                           ! decay mode N-Pi
    pwiddpi0=0.3*0.2;
    pwiddnpi0=0.117;//            ! the width of Delta 
    pwidth2(sp2,smr,xmd,fgMpi,1,pwiddpi0,fgMnucl,fgMpi,1,pwiddnpi0,pwiddpi);//      ! decay mode Delta-pi
    pwidns0=0.3*0.15;
    pwidspi0=0.6;//               ! the width of Sigma->Pi Pi
    pwidth2(sp2,smr,xmsigma,fgMnucl,0,pwidns0,fgMpi,fgMpi,0,pwidspi0,pwidns);//     ! decay mode N-Sigma
    
    wid= pwidnpi +pwiddpi + pwidns;
    break;
  case kD13_1520:   
    smr=xmd13;        
    pwidnpi0=0.115*0.6;
    pwidth1(sp2,smr,fgMnucl,fgMpi,2,pwidnpi0,pwidnpi);//                             ! n_pi
    pwiddpi00=0.115*0.15;
    pwiddpi20=0.115*0.125;     
    pwiddnpi0=0.117;
    pwidth2(sp2,smr,xmd,fgMpi,0,pwiddpi00,fgMnucl,fgMpi,1,pwiddnpi0,pwiddpi0);//    ! decay mode Delta-pi L=0
    pwidth2(sp2,smr,xmd,fgMpi,2,pwiddpi20,fgMnucl,fgMpi,1,pwiddnpi0,pwiddpi2);//    ! decay mode Delta-pi L=2
    pwidnr00=0.115*0.09;
    pwidnr20=0.115*0.035;
    pwidrpi0=0.1491;//      ! the width of Rho->pi pi P-wave 
    pwidth2(sp2,smr,xmrho,fgMnucl,0,pwidnr00,fgMpi,fgMpi,1,pwidrpi0,pwidnr0);//       ! decay mode N-Rho L=0
    pwidth2(sp2,smr,xmrho,fgMnucl,2,pwidnr20,fgMpi,fgMpi,1,pwidrpi0,pwidnr2);//       ! decay mode N-Rho L=2

    wid= pwidnpi + pwiddpi0 + pwiddpi2 + pwidnr0 + pwidnr2;
    break;
  case kS11_1535:
    smr=xms11;
    pwidnpi0=0.15*0.45;
    pwidth1(sp2,smr,fgMnucl,fgMpi,0,pwidnpi0,pwidnpi);//                             ! N-Pi
    pwidne0=0.15*0.42;
    pwidth1(sp2,smr,fgMnucl,xmeta,0,pwidne0,pwidne);//                              ! N-eta
    pwidnspi0=0.15*0.08;
    pwidns0=0.;//   ! call the width of P11(1440)
    pwidth2(sp2,smr,xmp11,fgMpi,0,pwidnspi0,fgMnucl,fgMpi,1,pwidns0,pwidnspi);//      ! N*(1440)-Pi   
    pwiddpi0=0.15*0.01;
    pwiddnpi0=0.117;
    pwidth2(sp2,smr,xmd,fgMpi,2,pwiddpi0,fgMnucl,fgMpi,1,pwiddnpi0,pwiddpi);//     ! Delta-Pi L=2
    pwidnr0=0.15*0.02;
    pwidrpi0=0.1491;
    pwidth2(sp2,smr,xmrho,fgMnucl,0,pwidnr0,fgMpi,fgMpi,1,pwidrpi0,pwidnr);//       ! decay mode N-Rho L=2
    pwidns0=0.15*0.02;
    pwidspi0=0.6;//               ! the width of Sigma->Pi Pi
    pwidth2(sp2,smr,xmsigma,fgMnucl,0,pwidns0,fgMpi,fgMpi,0,pwidspi0,pwidns);//     ! decay mode N-Sigma       
    wid=pwidnpi + pwidne  + pwiddpi + pwidnspi + pwidnr + pwidns;
    break;
  }

  return wid;
}

TComplex pwidth1(sp2,smr,sma,smb,l,pwid0){//   ! Both of the particle A and B are stable
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,external :: pcm     

  if(sp2 < (sma+smb)*(sma+smb)){
    pwid=0.;
  }else{
    sw=TMath::Sqrt(sp2);
    wpcm=pcm(sw,sma,smb);
    rpcm=pcm(smr,sma,smb);
    rho=wpcm**(2*l+1)/sw;
    rho0=rpcm**(2*l+1)/smr;
    pwid=pwid0*rho/rho0;
  }
  return pwid;
}

TComplex pwidth2(sp2,smr,sma,smb,l,pwid1,smc,smd,l2,pwid2)//   ! Both of the particle A and B are stable
  // implicit real*8 (a,b,d-h,o-z)
  // implicit complex*16 (c)
  // real*8,external :: pcm   
  
  if(sp2 < (smc+smd+smb)*(smc+smd+smb)){
    pwid=0.;
  }else{
    sw=TMath::Sqrt(sp2);
    rho_width(sw,sma,smb,l,smc,smd,l2,pwid2,rho);
    rho_width(smr,sma,smb,l,smc,smd,l2,pwid2,rho0);
    pwid=pwid1*rho/rho0 ;
  }
return pwid;
}

       SUBROUTINE rho_width(sw,sma,smb,l,smc,smd,l2,pwid2,rho)  
  implicit real*8 (a,b,d-h,o-z)
  implicit complex*16 (c)
  real*8,dimension(200) :: x,f1
  real*8,external :: pcm   
  pi=acos(-1.d0)  

  pmin=0.d0
  pmax=pcm(sw,smc+smd,smb)
  n=1
  np=20*n
  call DSG20r(pmin,pmax,n,x,np)
  do i=1,np
     f1(i)=fkernel_rho(x(i))
     ! print *,f1(i)
  end do
  ! stop
  call DRG20r(pmin,pmax,n,f1,rho)

contains

  FUNCTION fkernel_rho(pab)
    implicit real*8 (a,b,d-h,o-z)
    implicit complex*16 (c)
    ! real*8,dimension(4,5) :: swidd
    swa=sqrt(sw**2+smb**2-2.d0*sw*sqrt(pab**2+smb**2) ) 
    if(abs(sma-1.44).gt.1d-5) then
       call pwidth1(swa**2,sma,smc,smd,l2,pwid2,pwida)
       ! print *,swa,sma,smc,smd,l2,pwid2,pwida
       ! stop
    else 
       call width(2,swa**2,pwida) 
    endif

    fkernel_rho=2.d0/pi*pab**(2*l+2)/sqrt(pab**2+smb**2)*sma*pwida  &
         /( (swa**2-sma**2)**2 + sma**2*pwida**2 ) 
    !    print*,xwa,sma,xmc,xmdd,l2,pwid2,pwida,fkernel_rho,pab
    !    stop
    return
  end function fkernel_rho

END SUBROUTINE rho_width
//____________________________________________________________________________
/*
// This in theory should be useless
void LARNCGammaXSec::AmpNum2(double t1, double t2, TensorOrder2& c_lh_both_p, TensorOrder2& c_lh_both_n)
{
c_p  = SlashDirac(xp);
c_pp = SlashDirac(xpp);
c_k  = SlashDirac(xk);
c_kp = SlashDirac(xkp);
c_q  = SlashDirac(xq);
c_qf = SlashDirac(xqf);

TensorOrder4& ch_ver;
TensorOrder4& ch_ver_t;
TensorOrder4& ch_ver_ab;
TensorOrder4& ch_ver_abt;
TensorOrder4& ch_ver_cd;
TensorOrder4& ch_ver_cdt;
TensorOrder4& ch_ver_e;
TensorOrder4& ch_ver_et;
TensorOrder4& ch_ver_d13;
TensorOrder4& ch_ver_d13t;
TensorOrder4& ch_ver_p11;
TensorOrder4& ch_ver_p11t;
TensorOrder4& ch_ver_s11;
TensorOrder4& ch_ver_s11t;
  
cpm  = c_p  + fgMn * GetIdentityMatrix();
cppm = c_pp + fgMn * GetIdentityMatrix();

if(fNDiag == 7){
fNucleon=-1;
LARNCGammaXSec::VertexAB(t1,t2,ch_ver_ab_n,ch_ver_abt_n);
fNucleon=1;
LARNCGammaXSec::VertexAB(t1,t2,ch_ver_ab_p,ch_ver_abt_p);

LARNCGammaXSec::VertexCD(t1,t2,ch_ver_cd,ch_ver_cdt);
LARNCGammaXSec::VertexE(ch_ver_e);
ch_ver_et=-ch_ver_e;
     
ch_ver_p   = ch_ver_cd  + ch_ver_e  + ch_ver_ab_p;
ch_ver_t_p = ch_ver_cdt + ch_ver_et + ch_ver_abt_p;

ch_ver_n   = ch_ver_cd  - ch_ver_e  + ch_ver_ab_n;
ch_ver_t_n = ch_ver_cdt - ch_ver_et + ch_ver_abt_n;

}else if(fNDiag == 0){
fNucleon = -1;
LARNCGammaXSec::VertexAB(t1, t2, ch_ver_ab_n, ch_ver_abt_n);
fNParity = -1;
LARNCGammaXSec::VertexJ32(t2, ch_ver_d13_n, ch_ver_d13t_n);
fNParity = 1;
LARNCGammaXSec::VertexJ12(t2, ch_ver_p11_n, ch_ver_p11t_n);
fNParity = -1;
LARNCGammaXSec::VertexJ12(t2, ch_ver_s11_n, ch_ver_s11t_n);
fNucleon = 1;
LARNCGammaXSec::VertexAB(t1, t2, ch_ver_ab_p, ch_ver_abt_p);
fNParity = -1;
LARNCGammaXSec::VertexJ32(t2, ch_ver_d13_p, ch_ver_d13t_p);
fNParity = 1;
LARNCGammaXSec::VertexJ12(t2, ch_ver_p11_p, ch_ver_p11t_p);
fNParity = -1;
LARNCGammaXSec::VertexJ12(t2, ch_ver_s11_p, ch_ver_s11t_p);

LARNCGammaXSec::VertexCD(t1, t2, ch_ver_cd, ch_ver_cdt);
LARNCGammaXSec::VertexE(ch_ver_e);
ch_ver_et = -ch_ver_e;

ch_ver_p   = ch_ver_cd  + ch_ver_e  + ch_ver_ab_p  + ch_ver_d13_p  + ch_ver_p11_p  + ch_ver_s11_p;
ch_ver_t_p = ch_ver_cdt + ch_ver_et + ch_ver_abt_p + ch_ver_d13t_p + ch_ver_p11t_p + ch_ver_s11t_p;

ch_ver_n   = ch_ver_cd  - ch_ver_e  + ch_ver_ab_n  + ch_ver_d13_n  + ch_ver_p11_n  + ch_ver_s11_n;
ch_ver_t_n = ch_ver_cdt - ch_ver_et + ch_ver_abt_n + ch_ver_d13t_n + ch_ver_p11t_n + ch_ver_s11t_n;

}else{
std::cout << "FNDiag should be 0 or 7" <<std::endl;
exit(1);
}

c_lh_p = 0;
c_lh_n = 0;
c_lh_anti_p = 0;
c_lh_anti_n = 0;
LARNCGammaXSec::TraceLight(c_tr_l,c_tr_l_anti);

for(int n_alpha = 0; n_alpha < 4; n_alpha++){
for(int n_beta = 0; n_beta < 4; n_beta++){
c_mult2(cpm,ch_ver_t_p,cppm,ch_ver_p, n_alpha,n_beta,c_trh_p);
c_mult2(cpm,ch_ver_t_n,cppm,ch_ver_n, n_alpha,n_beta,c_trh_n);

c_lh_p=c_lh_p + c_tr_l(n_alpha,n_beta)*c_trh_p*g(n_alpha,n_alpha)*g(n_beta,n_beta);
c_lh_n=c_lh_n + c_tr_l(n_alpha,n_beta)*c_trh_n*g(n_alpha,n_alpha)*g(n_beta,n_beta);

c_lh_anti_p=c_lh_anti_p + c_tr_l_anti(n_alpha,n_beta)*c_trh_p*g(n_alpha,n_alpha)*g(n_beta,n_beta);
c_lh_anti_n=c_lh_anti_n + c_tr_l_anti(n_alpha,n_beta)*c_trh_n*g(n_alpha,n_alpha)*g(n_beta,n_beta);
}
}

if(fNCheck == 520){
if(abs(aimag(c_lh_p)) > TMath::Power(10, -6) ||
abs(aimag(c_lh_anti_p)) > TMath::Power(10, -6) ||
abs(aimag(c_lh_n)) > TMath::Power(10, -6) ||
abs(aimag(c_lh_anti_n)) > TMath::Power(10, -6)){
std::cout << "The image part of c_lh or c_lh_anti is not zero" << std::endl;
std::cout << "c_lh_p, c_lh_anti_p, c_lh_n, c_lh_anti_n" << std::endl;
std::cout << c_lh_p << "  " << c_lh_anti_p << "  " << c_lh_n << "  " << c_lh_anti_n << std::endl;
std::cout << c_lh_anti_p << "  " << c_lh_n <<  "  " << c_lh_anti_n << std::endl;
exit(1);
}
}

c_lh_both_p = real(c_lh_p) + cu * real(c_lh_anti_p);
c_lh_both_n = real(c_lh_n) + cu * real(c_lh_anti_n) ;

if(fNCheck == 520){
if(real(c_lh_both_p) < 0  ||
aimag(c_lh_both_p) < 0 ||
real(c_lh_both_n) < 0  ||
aimag(c_lh_both_n) < 0){
std::cout << "the LH tensor is less than zero" << std::endl;
std::cout << "c_lh_both_p, c_lh_both_n" << std::endl;
std::cout << c_lh_both_p << "  " << c_lh_both_n << std::endl;
exit(1);
}
}
  
return;
}
*/

