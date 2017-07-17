//____________________________________________________________________________
/*
   Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
   or see $GENIE/LICENSE

Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
STFC, Rutherford Appleton Laboratory - March 11, 2005

For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Conventions/RefFrame.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "BergerSehgal/BergerSehgalFMCOHPiPXSec2015.h"
#include "Utils/HadXSUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
BergerSehgalFMCOHPiPXSec2015::BergerSehgalFMCOHPiPXSec2015() :
  XSecAlgorithmI("genie::BergerSehgalFMCOHPiPXSec2015")
{

}
//____________________________________________________________________________
BergerSehgalFMCOHPiPXSec2015::BergerSehgalFMCOHPiPXSec2015(string config) :
  XSecAlgorithmI("genie::BergerSehgalFMCOHPiPXSec2015", config)
{

}
//____________________________________________________________________________
BergerSehgalFMCOHPiPXSec2015::~BergerSehgalFMCOHPiPXSec2015()
{

}
//____________________________________________________________________________
double BergerSehgalFMCOHPiPXSec2015::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // Here we are following PRD 79, 053003 (2009) by Berger and Sehgal
  // This method computes the differential cross section represented 
  // in Eq.'s 6 (CC) and 7 (NC) from that paper.

  // We have additionally modified the formulae to account for a 
  // non-infinite mass for the target nucleus

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double xsec = 0.0;

  double E      = init_state.ProbeE(kRfLab);        // nu E
  double Q2     = kinematics.Q2();
  double y      = kinematics.y();                   // inelasticity
  double t      = kinematics.t();                   // fun exists?
  assert(E > 0.);
  assert(y > 0.);
  assert(y < 1.);
  double ppistar = PionCOMAbsMomentum(interaction); // |Center of Mass Momentum|
  if (ppistar <= 0.0) { 
    LOG("BergerSehgalFMCohPi", pDEBUG) << "Pion COM momentum negative for Q2 = " << Q2 << 
      " y = " << y; 
    return 0.0; 
  }
  double front  = ExactKinematicTerm(interaction);
  if (front <= 0.0) { 
    LOG("BergerSehgalFMCohPi", pDEBUG) << "Exact kin. form = " << front << 
      " E = " << E << " Q2 = " << Q2 << " y = " << y;
    return 0.0; 
  }

  double A      = (double) init_state.Tgt().A();   // mass number
  double A2     = TMath::Power(A, 2.);
  double A_3    = TMath::Power(A, 1./3.);
  double M      = init_state.Tgt().Mass();
  double M_pi   = pionIsCharged ? kPionMass : kPi0Mass;
  double M_pi2  = M_pi * M_pi;
  double Epi    = y * E - t / (2 * M);  
  double Epi2   = Epi * Epi;
  double ma2    = fMa * fMa;
  double Ga     = ma2 / (ma2 + Q2);
  double Ga2    = Ga * Ga;
  double Ro2    = TMath::Power(fRo * units::fermi, 2.);
  double ppi2   = Epi2 - M_pi2;
  double ppi    = ppi2 > 0.0 ? sqrt(ppi2) : 0.0;
  // double fp     = 0.93 * kPionMass;  // unused  // pion decay constant

  double costheta = (t - Q2 - M_pi * M_pi) / (2 * ( (y *E - Epi) * Epi - 
       ppi * sqrt(TMath::Power(y * E - Epi, 2.) + t) ) );

  if ((costheta > 1.0) || (costheta < -1.0)) return 0.0;

  // tot. pi+N xsec
  double sTot   = 
    utils::hadxs::berger::TotalPionNucleonXSec(Epi, pionIsCharged); 
  double sTot2  = sTot * sTot;
  // inel. pi+N xsec
  double sInel  = 
    utils::hadxs::berger::InelasticPionNucleonXSec(Epi, pionIsCharged); 

  // Fabs (F_{abs}) describes the average attenuation of a pion emerging
  // from a sphere of nuclear matter with radius = R_0 A^{1/3}. it is 
  // Eq. 13 in Berger-Sehgal PRD 79, 053003
  double Fabs_input = (9.0 * A_3) / (16.0 * kPi * Ro2);
  double Fabs       = TMath::Exp( -1.0 * Fabs_input * sInel);

  // A_RS for BS version of RS, and/or Tpi>1.0
  double RS_factor = (A2 * Fabs) / (16.0 * kPi) * (sTot2); 
  double R         = fRo * A_3 * units::fermi; // nuclear radius
  double R2        = R * R;                    // 
  double b         = 0.33333 * R2;             // Eq. 12 in BS
  double expbt     = TMath::Exp( -b * t );
  double dsigEldt  = sTot2 / (16. * kPi);           // Eq. 11 in BS
  double dsigpiNdt = A2 * dsigEldt * expbt * Fabs;  // Eq. 10 in BS

  double tpilow    = 0.0;
  double siglow    = 0.0;
  double tpihigh   = 0.0;
  double sighigh   = 0.0;
  double dsigdt    = 0.0;
  double tpi       = 0.0;
  int    xsec_stat = 0;

  // c.o.m.
  tpi       = Epi - M_pi;
  xsec_stat = 
    utils::hadxs::berger::PionNucleusXSec(
        tpi, ppistar, t, A, 
        tpilow, siglow, 
        tpihigh, sighigh);
  dsigdt = siglow + (sighigh - siglow) * (tpi - tpilow) / (tpihigh - tpilow);
  dsigdt = dsigdt / (2.0 * ppistar * ppistar) * units::mb;

  double edep_dsigpiNdt = Epi > 1.0 ? dsigpiNdt : dsigdt;

  xsec = front * Ga2 * edep_dsigpiNdt;

  // Correction for finite final state lepton mass.
  // Lepton mass modification is part of Berger-Sehgal and is not optional.
  if (pionIsCharged) {
    double C = 1.;
    // First, we need to remove the leading G_{A}^2 which is required for NC.
    xsec /= Ga2;
    // Next, \cos^2 \theta_{Cabbibo} appears in the CC xsec, but not the NC.
    xsec *= fCos8c2;
    double ml    = interaction->FSPrimLepton()->Mass();
    double ml2   = TMath::Power(ml,2);
    double Q2min = ml2 * y/(1-y);
    if(Q2 > Q2min) {
      double C1 = TMath::Power(Ga - 0.5 * Q2min / (Q2 + kPionMass2), 2);
      double C2 = 0.25 * y * Q2min * (Q2 - Q2min) / 
        TMath::Power(Q2 + kPionMass2, 2);
      C = C1 + C2;
    } else {
      C = 0.;
    }
    xsec *= (2. * C); // *2 is for CC vs NC in BS 
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BergerSehgalFMCohPi", pDEBUG)
    << "\n momentum transfer .............. Q2    = " << Q2
    << "\n mass number .................... A     = " << A
    << "\n pion energy .................... Epi   = " << Epi
    << "\n propagator term ................ propg = " << propg
    << "\n Re/Im of fwd pion scat. ampl. .. Re/Im = " << fReIm
    << "\n total pi+N cross section ....... sigT  = " << sTot
    << "\n inelastic pi+N cross section ... sigI  = " << sInel
    << "\n nuclear size scale ............. Ro    = " << fRo
    << "\n pion absorption factor ......... Fabs  = " << Fabs
    << "\n t integration factor ........... tint  = " << tint;
  LOG("BergerSehgalFMCohPi", pINFO)
    << "d2xsec/dxdy[COHPi] (x= " << x << ", y="
    << y << ", E=" << E << ") = "<< xsec;
#endif

  //----- The algorithm computes d^3xsec/dQ^2dydt
  //      Check whether Jacobian tranformation is needed...

  return xsec;
}
//____________________________________________________________________________
double BergerSehgalFMCOHPiPXSec2015::ExactKinematicTerm(
    const Interaction * interaction) const
{
  // This function is a bit inefficient but is being encapsulated as 
  // such in order to possibly migrate into a general kinematics check.
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool   pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double M_pi          = pionIsCharged ? kPionMass : kPi0Mass;
  double E             = init_state.ProbeE(kRfLab);        // nu E
  double Q2            = kinematics.Q2();
  double y             = kinematics.y();                   // inelasticity
  double fp2           = (0.93 * M_pi)*(0.93 * M_pi); 

  double term = ((kGF2 * fp2) / (4.0 * kPi2)) * 
    ((E * (1.0 - y)) / sqrt(y*E * y*E + Q2)) * 
    (1.0 - Q2 / (4.0 * E*E * (1.0 - y)));
  return term;   
}
//____________________________________________________________________________
double BergerSehgalFMCOHPiPXSec2015::PionCOMAbsMomentum(
    const Interaction * interaction) const
{
  // This function is a bit inefficient but is being encapsulated as 
  // such in order to possibly migrate into a general kinematics check.
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool   pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double M_pi          = pionIsCharged ? kPionMass : kPi0Mass;
  double E             = init_state.ProbeE(kRfLab);        // nu E
  double Q2            = kinematics.Q2();
  double y             = kinematics.y();                   // inelasticity
  double MT            = init_state.Tgt().Mass(); 

  double W2      = MT * MT - Q2 + 2.0 * y * E * MT;
  double arg     = (2.0 * MT * (y * E - M_pi) - Q2 - M_pi * M_pi) * 
    (2.0 * MT * (y * E + M_pi) - Q2 - M_pi * M_pi);
  if (arg < 0) return arg;
  double ppistar = TMath::Sqrt(arg) / 2.0 / TMath::Sqrt(W2);

  return ppistar;
}
//____________________________________________________________________________
double BergerSehgalFMCOHPiPXSec2015::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BergerSehgalFMCOHPiPXSec2015::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  int nu = init_state.ProbePdg();

  if (!proc_info.IsCoherent())  return false;
  if (!proc_info.IsWeak())      return false;
  if (target.HitNucIsSet())     return false;
  if (!target.A()>1)            return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void BergerSehgalFMCOHPiPXSec2015::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalFMCOHPiPXSec2015::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalFMCOHPiPXSec2015::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fMa         = fConfig->GetDoubleDef("Ma",           gc->GetDouble("COH-Ma"));
  fRo         = fConfig->GetDoubleDef("Ro",           gc->GetDouble("COH-Ro"));
  double thc  = fConfig->GetDoubleDef("CabbiboAngle", gc->GetDouble("CabbiboAngle"));
  fCos8c2     = TMath::Power(TMath::Cos(thc), 2);

  // fRSPionXSec => Do not use the pion-nucleus cross section from Table 1 in PRD 79, 053003
  // Instead, use the Rein-Sehgal "style" pion-nucleon cross section and scale by A 
  // for all pion energies.
  fRSPionXSec = fConfig->GetBoolDef("UseRSPionXSec", gc->GetBool("COH-UseRSPionXSec"));

  //-- load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

