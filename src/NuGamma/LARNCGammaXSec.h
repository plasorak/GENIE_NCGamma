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

  enum Nucleon{
    kProton  = 0,
    kNeutron = 1
  };
    
  class LARNCGammaXSec {
    
  public:
    LARNCGammaXSec();
    ~LARNCGammaXSec();

    //-- XSecAlgorithmI interface implementation
    double DiffXSec(const Interaction * i);
    double Integral(const Interaction * i);

  private:
    
    mutable TLorentzVector *gNeutrinoInit;
    mutable TLorentzVector *gTargetInit;
    mutable TLorentzVector *gHadronSystFinal;
    mutable TLorentzVector *gNeutrinoFinal;
    mutable TLorentzVector *gRemnantNucleon;
    mutable TLorentzVector *gPhotonFinal;

    void dxsec(const double Enu, const double W, const double Qsq,
               const int nucl, const int nd,
               double& dsigmaneutrino, double& sigmaantineutrino);
    void dcross2(const int n,
                 std::complex<double>& cdc);

    void momentum(const double xkp0,     const double theta_kp, const double xqf0,
                  const double theta_qf, const double phi_qf,   const double xp0,
                  const double theta_p,  const double phi_p,
                  double&      t1,       double&      t2);
    void fgaus(const int N, int* JS, double* X,
               void (*)(const int J, double& DN, double& UP),
               std::complex<double> (*)(const int n, double* x),
               std::complex<double>& S);
    void FS(const int J, double& DN, double& UP);
    std::complex<double> ckernel(const int n, double* x);
  
    // Amplitude calculation
    void AmpNum(const double t1, const double t2, TComplex& c_lh_both);
    // void AmpNum2(double t1, double t2, TensorUtils::TensorOrder2& c_lh_both_p, TensorUtils::TensorOrder2& c_lh_both_n, LARContainer *cont);

    // Form factors
    void ElectroMagneticFormFactor(const double t, const int nucleon,
                                   double &f1,  double &f2);
    void NeutralCurrentFormFactor (const double t, const int nucleon,
                                   double &fv1, double &fv2, double &fva);
      
    void P33_1232FormFactor(const double q2,
                            double &c3v, double &c4v, double &c5v,
                            double &c3a, double &c4a, double &c5a, double &c6a);
    void P11_1440FormFactor(const double t1,
                            TensorUtils::TensorOrder2& fem,
                            TensorUtils::TensorOrder2& fvem,
                            TensorUtils::TensorOrder2& fvnc,
                            TensorUtils::TensorOrder1& fanc);
    void D13_1520FormFactor(const double t1,
                            TensorUtils::TensorOrder2& fem,
                            TensorUtils::TensorOrder2& fvem,
                            TensorUtils::TensorOrder2& fvnc,
                            TensorUtils::TensorOrder2& fanc);
    void P33_1232FormFactor(const double t1,
                            TensorUtils::TensorOrder2& fem,
                            TensorUtils::TensorOrder2& fvem,
                            TensorUtils::TensorOrder2& fvnc,
                            TensorUtils::TensorOrder2& fanc);
    void S11_1535FormFactor(const double t1,
                            TensorUtils::TensorOrder2& fem,
                            TensorUtils::TensorOrder2& fvem,
                            TensorUtils::TensorOrder2& fvnc,
                            TensorUtils::TensorOrder1& fanc);
    
    void HelicityAmplitude(const Resonance_t res, const std::string wave, const double q2,
                           double &a12, double &a32, double &s12);
    void EMtoNCFormFactor(const double n, const TensorUtils::TensorOrder2& fvem,
                          TensorUtils::TensorOrder2& fvnc);
    double DeltaPropagator(const TensorUtils::TensorOrder1& am);
    TComplex Propagator(const Resonance_t nexcit, const TensorUtils::TensorOrder1& sp);

    // The "vertex" calculation
    void VertexAB (const double t1, const double t2,
                   TensorUtils::TensorOrder4& ch_vertex1,
                   TensorUtils::TensorOrder4& ch_vertex2);
    void VertexCD (const double t1, const double t2,
                   TensorUtils::TensorOrder4& ch_vertex1,
                   TensorUtils::TensorOrder4& ch_vertex2);
    void VertexE  (TensorUtils::TensorOrder4& c_lh_both);

    void VertexJ12(const double t1,
                   TensorUtils::TensorOrder4& ch_verj12,
                   TensorUtils::TensorOrder4& ch_verj12_t);
    void VertexJ32(const double t1,
                   TensorUtils::TensorOrder4& ch_verj32,
                   TensorUtils::TensorOrder4& ch_verj32_t);
    void Vertex12 (const double f1, const double f2, const double fa,
                   const TensorUtils::TensorOrder1& sq,
                   TensorUtils::TensorOrder3& ver12);
    void Vertex32 (const TensorUtils::TensorOrder1& fcv,
                   const TensorUtils::TensorOrder1& fca,
                   const TensorUtils::TensorOrder1& sp,
                   const TensorUtils::TensorOrder1& sq,
                   TensorUtils::TensorOrder4& cver32);

    void TraceLight(TensorUtils::TensorOrder2& c_tr_l,
                    TensorUtils::TensorOrder2& c_tr_l_anti);
    double DeltaPi(const TensorUtils::TensorOrder1& sp);
    void AEM(const double t1, const double t2,
             const TensorUtils::TensorOrder1& sp, const TensorUtils::TensorOrder1& sq,
             TensorUtils::TensorOrder4& caem);
    void ANC(const double t1, const double t2,
             const TensorUtils::TensorOrder1& sp, const TensorUtils::TensorOrder1& sq,
             TensorUtils::TensorOrder4& canc);
    TComplex cDdelta(const TensorUtils::TensorOrder1& spp);
    TComplex Width(const Resonance_t nexcit, const double sp2);
    void FactorC(const double t1, const double t2,
                 TensorUtils::TensorOrder1& fcv,
                 TensorUtils::TensorOrder1& fcvt,
                 TensorUtils::TensorOrder1& fcat);
    TensorUtils::TensorOrder2* Lambda(const int ns, const int nd,
                                      const TensorUtils::TensorOrder1& ppd);
    double SMomentum(const TensorUtils::TensorOrder1& sp);
    double Flam(const double sx, const double sy, const double sz);
    TensorUtils::TensorOrder2* Dim3to2(const TensorUtils::TensorOrder3& tensor,
                                       const int n);
    TensorUtils::TensorOrder2* Dim4to2(const TensorUtils::TensorOrder4& tensor,
                                       const int n1, const int n2);
    
    TensorUtils::TensorOrder2* Mult3Matrices(const TensorUtils::TensorOrder2& mat1,
                                             const TensorUtils::TensorOrder2& mat2,
                                             const TensorUtils::TensorOrder2& mat3);
    
    TensorUtils::TensorOrder2* MatMult(const TensorUtils::TensorOrder2& mat1,
                                       const TensorUtils::TensorOrder2& mat2);
  
    TComplex Mult2(const TensorUtils::TensorOrder2& c_a,
                   const TensorUtils::TensorOrder4& c_b,
                   const TensorUtils::TensorOrder2& c_c,
                   const TensorUtils::TensorOrder4& c_d,
                   const int n_alpha, const int n_beta);
    
    double pcm(const double sw, const double sma, const double smb){
      return TMath::Sqrt(TMath::Power(sw*sw-sma*sma-smb*smb,2)-4.*sma*sma*smb*smb )/2./sw;
  };
  double fkernel_rho(const double pab, const double sw,
                     const double sma, const double smb, const int l,
                     const double smc, const double smd, const int l2,
                     const double pwid2);
  double rho_width(const double sw,  const double sma, const double smb, const int l,
                   const double smc, const double smd, const int l2,
                   const double pwid2);
    void pwidth2(const double sp2, const double smr,
                 const double sma, const double smb, const int l,  const double wid1,
                 const double smc, const double smd, const int l2, const double pwid2,
                 double& pwid);
    void pwidth1(const double sp2, const double smr,
                 const double sma, const double smb, const int l, const double pwid0,
                 double& pwid);

    //double FMV(TensorUtils::TensorOrder1& xp1, TensorUtils::TensorOrder1& xp2);

    // For now, I want to reproduce the result from the paper, so I use their numbers. We can change that once we are sure the code is doing what we expect it to do.
    static const double fgVMassSq = 0.71;    // GeV
    static const double fgPMagMom = 2.793;   // Magnetic moment proton
    static const double fgNMagMom = -1.913;  // Magnetic moment neutron
    static double gev_bar;
    static double gfermi;
    static const double fgAp = 0.942;
    static const double fgBp = 4.61;
    static const double fgNuclMa = 1.;       // Nucleon axial mass
    static const double fgResMa = 1.;        // Resonance axial mass
    static const double fgNDeltaMa = 0.93;   // N-Delta axial mass
    static const double fg2Ma2 = 1;          // 1.012**2! the heavier resonances axial mass
    static const double fgMaSq = 0.93*0.93;  // ! 1.012**2! For the factor of Dig CD, xma=1.05, 1.012 
    static const double fgDeltaS = 0.;        // Something related to strangeness in nucleon
    static const double fgAlpha = 1./137.;    // fine structure constant
    static const double fgMdelta = 1.232;     // delta mass
    static const double fgMP11 = 1.44;        // mass of P11(1440)
    static const double fgMD13 = 1.52;        // mass of D13(1530)
    static const double fgMSigma = 0.475;     // mass ofSigma
    static const double fgMEta=0.548;//                         ! the mass of eta
    static const double fgMRho=0.77549;
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


    // Integration parameters:
    static const int nprecise = 5;
    int nnstep;

    
    double W, Q2, EGamma, PhiGamma;

    double Enu;
    
    
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
