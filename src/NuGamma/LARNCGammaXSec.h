//____________________________________________________________________________
/*

  \class    genie::LARNCGammaXSec

  \brief    Single gamma production from resonance.
  Is a concrete implementation of the XSecAlgorithmI interface

  \ref      E. Wang, L. Alvarez-Ruso, and J. Nieves
  Phys. Rev. C 89, 015503 (2014)

  \author   Pierre Lasorak <p.j.j.lasorak \at qmul.ac.uk>
  Queen Mary University, London

  \created  November 24, 2015

  \cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LAR_NUGAMMA_XSEC_H_
#define _LAR_NUGAMMA_XSEC_H_

#include "TLorentzVector.h"
#include "Base/XSecAlgorithmI.h"
#include "NuGamma/TensorInc.h"

namespace genie {

  class LARNCGammaXSec : public XSecAlgorithmI {

  public:
    LARNCGammaXSec();
    LARNCGammaXSec(string config);
    virtual ~LARNCGammaXSec();

    //-- XSecAlgorithmI interface implementation
    double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
    double Integral        (const Interaction * i) const;
    bool   ValidProcess    (const Interaction * i) const;
    bool   ValidKinematics (const Interaction * i) const;
  
    //-- overload the Algorithm::Configure() methods to load private data
    //   members from configuration options
    void Configure(const Registry & config);
    void Configure(string config);
    void LoadConfig(void);

  private:
    
    mutable TLorentzVector *gNeutrinoInit;
    mutable TLorentzVector *gTargetInit;
    mutable TLorentzVector *gHadronSystFinal;
    mutable TLorentzVector *gNeutrinoFinal;
    mutable TLorentzVector *gRemnantNucleon;
    mutable TLorentzVector *gPhotonFinal;

    
    // Amplitude calculation
    void AmpNum(double t1, double t2, TensorUtils::TensorOrder2& c_lh_both);
    // void AmpNum2(double t1, double t2, TensorUtils::TensorOrder2& c_lh_both_p, TensorUtils::TensorOrder2& c_lh_both_n, LARContainer *cont);

    // Form factors
    void ElectroMagneticFormFactor(double t, int nucleon, double &f1,  double &f2);
    void NeutralCurrentFormFactor (double t, int nucleon, double &fv1, double &fv2, double &fva);
      
    void P33_1232FormFactor(double q2, double &c3v, double &c4v, double &c5v, double &c3a, double &c4a, double &c5a, double &c6a);
    void P11_1440FormFactor(double t1, TensorUtils::TensorOrder2& fem, TensorUtils::TensorOrder2& fvem, TensorUtils::TensorOrder2& fvnc, TensorUtils::TensorOrder1& fanc);
    void D13_1520FormFactor(double t1, TensorUtils::TensorOrder2& fem, TensorUtils::TensorOrder2& fvem, TensorUtils::TensorOrder2& fvnc, TensorUtils::TensorOrder2& fanc);
    void P33_1232FormFactor(double t1, TensorUtils::TensorOrder2& fem, TensorUtils::TensorOrder2& fvem, TensorUtils::TensorOrder2& fvnc, TensorUtils::TensorOrder2& fanc);
    void S11_1535FormFactor(double t1, TensorUtils::TensorOrder2& fem, TensorUtils::TensorOrder2& fvem, TensorUtils::TensorOrder2& fvnc, TensorUtils::TensorOrder1& fanc);
    
    void HelicityAmplitude(Resonance_t res, std::string wave, double q2, double &a12, double &a32, double &s12);
    void EMtoNCFormFactor(double n, TensorUtils::TensorOrder2&fvem, TensorUtils::TensorOrder2& fvnc);
    double DeltaPropagator(TensorUtils::TensorOrder1& am);
    TComplex Propagator(int nexcit, TensorUtils::TensorOrder1& sp);

    // The "vertex" calculation
    void VertexAB (double t1, double t2, TensorUtils::TensorOrder4& ch_vertex1, TensorUtils::TensorOrder4& ch_vertex2);
    void VertexCD (double t1, double t2, TensorUtils::TensorOrder4& ch_vertex1, TensorUtils::TensorOrder4& ch_vertex2);
    void VertexE  (TensorUtils::TensorOrder4& c_lh_both);

    void VertexJ12(double t1, TensorUtils::TensorOrder4& ch_verj12, TensorUtils::TensorOrder4& ch_verj12_t);
    void VertexJ32(double t1, TensorUtils::TensorOrder4& ch_verj32, TensorUtils::TensorOrder4& ch_verj32_t);
    void Vertex12 (double f1,  double f2, double fa, TensorUtils::TensorOrder1& sq, TensorUtils::TensorOrder3& ver12);
    void Vertex32 (TensorUtils::TensorOrder1& fcv, TensorUtils::TensorOrder1& fca, TensorUtils::TensorOrder1& sp, TensorUtils::TensorOrder1& sq, TensorUtils::TensorOrder4& cver32);

    void TraceLight(TensorUtils::TensorOrder2& c_tr_l, TensorUtils::TensorOrder2& c_tr_l_anti);
    double DeltaPi(TensorUtils::TensorOrder1& sp);
    void AEM(double t1, double t2, TensorUtils::TensorOrder1& sp, TensorUtils::TensorOrder1& sq, TensorUtils::TensorOrder4& caem);
    void ANC(double t1, double t2, TensorUtils::TensorOrder1& sp, TensorUtils::TensorOrder1& sq, TensorUtils::TensorOrder4& canc);
    TComplex cDdelta(TensorUtils::TensorOrder1& spp);
    TComplex Width(int nexcit, double sp2);
    void FactorC(double t1, double t2, TensorUtils::TensorOrder1& fcv, TensorUtils::TensorOrder1& fcvt, TensorUtils::TensorOrder1& fcat);
    TensorUtils::TensorOrder2* Lambda(int ns, int nd, TensorUtils::TensorOrder1& ppd);
    double SMomentum(TensorUtils::TensorOrder1& sp);
    double Flam(double sx, double sy, double sz);
    TensorUtils::TensorOrder2* Dim3to2(TensorUtils::TensorOrder3& tensor, int n);
    TensorUtils::TensorOrder2* Dim4to2(TensorUtils::TensorOrder4& tensor, int n1, int n2);
    
    TensorUtils::TensorOrder2* Mult3Matrices(TensorUtils::TensorOrder2& mat1,
                                             TensorUtils::TensorOrder2& mat2,
                                             TensorUtils::TensorOrder2& mat3);
    
    TensorUtils::TensorOrder2* MatMult(TensorUtils::TensorOrder2& mat1,
                                       TensorUtils::TensorOrder2& mat2);
  
    TComplex Mult2(TensorUtils::TensorOrder2& c_a, TensorUtils::TensorOrder4& c_b,
                   TensorUtils::TensorOrder2& c_c, TensorUtils::TensorOrder4& c_d,
                   int n_alpha,int n_beta);
    
    //double FMV(TensorUtils::TensorOrder1& xp1, TensorUtils::TensorOrder1& xp2);

    // For now, I want to reproduce the result from the paper, so I use their numbers. We can change that once we are sure the code is doing what we expect it to do.
    static const double fgVMassSq = 0.71;    // GeV
    static const double fgPMagMom = 2.793;   // Magnetic moment proton
    static const double fgNMagMom = -1.913;  // Magnetic moment neutron
    static const double fgAp = 0.942;
    static const double fgBp = 4.61;
    static const double fgNuclMa = 1.;       // Nucleon axial mass
    static const double fgResMa = 1.;        // Resonance axial mass
    static const double fgNDeltaMa = 0.93;   // N-Delta axial mass
    static const double fg2Ma2 = 1;          // 1.012**2              ! the heavier resonances axial mass
    static const double fgMaSq = 0.93*0.93;  // ! 1.012**2              ! For the factor of Dig CD, xma=1.05, 1.012 
    static const double fgDeltaS = 0.;        // Something related to strangeness in nucleon
    static const double fgAlpha = 1./137.;    // fine structure constant
    static const double fgMdelta = 1.232;     // delta mass
    static const double fgMP11 = 1.44;        // mass of P11(1440)
    static const double fgMD13 = 1.52;        // mass of D13(1530)
    static const double fgMS11 = 1.535;       // mass of S11(1535)
    static const double fgMnucl = 0.9382723;  // Nucleon mass
    static const double fgMpi = 0.139;        // Pion Mass
    static const double fgMpi2 = 0.138;       // The mass of pi meson (second edition), define the formfactor Fp
    static const double fgMp = 0.938272;      // Proton mass
    static const double fgMn = 0.939;         // Neutron mass
    static const double fgMnSq = 0.939*0.939; // Neutron mass sq
    static const double fgGa = -1.267;       // define the formfactor Fa, Diag E(1.26)
    static const double fgF1s_0 = 0;         // 0.53
    static const double fgFStar = 2.14;
    static const double fgF2s_0 = 0; 
    static const double fgMv2 = 0.71;        // 0.84                  0.71    For the Factor of Fig CD, xmv=0.84
    static const double fgSinThetaW2 = 0.231;
    static const double fgFPi = 0.0924;
    static const double fgFa0P11 = -0.47;
    static const double fgFa0S11 = -0.21;
    static const double fgFca5D13 = -2.14; // )!753)!2.14)         
    static const double fgFca5P33 = 1;
    static const int fgnFF = 4;      // Form factors of the helicity amplitude

    enum Nucleon{
      kProton  = 0,
      kNeutron = 1
    };

  private:
 
    double fGw;

    // Whatever was in the common block of the fortran code.
    // Includes:
    //  - some parameter for the cross section calculation
    //  - the 4vectors for all the particles
    //  - some stuff that I don't know

    // Parameters
    int fMDelta;
    int fNDiag;
    int fNParity;
    int fNCheck;
    int fNucleon;
  
    // 4vectors
    TensorUtils::TensorOrder1 xp;
    TensorUtils::TensorOrder1 xpp;
    TensorUtils::TensorOrder1 xpd;
    TensorUtils::TensorOrder1 xpdc;
    TensorUtils::TensorOrder1 xk;
    TensorUtils::TensorOrder1 xkp;
    TensorUtils::TensorOrder1 xq;
    TensorUtils::TensorOrder1 xqf;

    // Dirac slashed version of the above
    TensorUtils::TensorOrder2 c_p;
    TensorUtils::TensorOrder2 c_pp;
    TensorUtils::TensorOrder2 c_k;
    TensorUtils::TensorOrder2 c_kp;
    TensorUtils::TensorOrder2 c_q;
    TensorUtils::TensorOrder2 c_qf;
  
    TensorUtils::TensorOrder2 cpm;
    TensorUtils::TensorOrder2 cppm;
    TensorUtils::TensorOrder2 cpqm;
    TensorUtils::TensorOrder2 cppqm;
    
    TensorUtils::TensorOrder2 DM[5];
    TensorUtils::TensorOrder2 metric;
    TensorUtils::TensorOrder2 unity;

  };
} // genie namespace

#endif  // _LAR_NUGAMMA_XSEC_H_
