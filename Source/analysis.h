#ifndef _ANALYSIS_ZPRIME_H_
#define _ANALYSIS_ZPRIME_H_

#include "RootTuple.h"
#include "atlas_style.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <cmath> 
#include <TString.h>
#include <TH2.h>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TBrowser.h"
#include <fstream>
#include "TLegend.h"
#include <complex>
#include <iomanip>

class AnalysisZprime{
public:
  AnalysisZprime(const TString channel, const TString model, const double luminosity, const TString& inputFileName, const TString& weightsFileName, const TString& outputFileName);
  virtual ~AnalysisZprime();
  
protected:
  Long64_t TotalEvents();
  Long64_t IncrementEvent(Long64_t i);
  void SetupTreesForNewFile(const TString& s);
  void CleanUp();
  
  void SetupInputFiles();
  void SetupOutputFiles();
  void SetupWeightsFiles();
  
  void PreLoop();
  void Loop();
  void PostLoop();
  void EachEvent();
  void CreateHistograms();
  void MakeGraphs();
  void WriteHistograms();
  void GetResults();
  void CheckPerformance();

  double TotalAsymmetry(TH1D* h_A, TH1D* h_B);
  void TotalSpinAsymmetries();
  TH1D* Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B);
  void ApplyLuminosity(TH1D*);
  void AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B);
  TH1D* MakeALL();
  TH1D* MakeAL();
  vector<std::complex<double> > SolveQuadratic(double a, double b, double c);
  std::vector<TLorentzVector> ReconstructSemiLeptonic(std::vector<TLorentzVector> p, int l_Q);
  
  bool PassCuts();
  bool PassCutsMET();
  bool PassCutsMtt();
  bool PassCutsFiducial();
  bool PassCutsYtt();

  void ResetCounters();
  void InitialiseCutflow();
  void PrintCutflow();
  const void UpdateCutflow(int cut, bool passed);

  static inline void ProgressBar(unsigned int x, unsigned int n, unsigned int w);
     
private:
  AnalysisZprime();
  AnalysisZprime(const AnalysisZprime& rhs);  
  void operator = (const AnalysisZprime& rhs);  

  // Counters
  bool m_useLumi;
  unsigned int m_nReco;
  unsigned int m_nQuarksMatched;
  unsigned int m_nNeutrinoMatched;
  unsigned int m_nRealRoots;
  unsigned int m_nComplexRoots;
  
  // Parameters
  float m_pi;
  float m_GeV;
  double m_luminosity;
  double m_Wmass;
  double m_tmass;

  double m_sigma;
  vector<double> m_weights;

  // Strings
  TString m_channel;  
  TString m_model;
  TString m_inputFileName;
  TString m_weightsFileName;
  TString m_outputFileName;

  // Cuts
  enum m_cutlist{
    c_Event,
    c_Mtt,
    c_MET,
    c_Ytt,
    c_Fiducial,
    m_cuts // Keep as last entry
  };

  // Input data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;

  // Cutflow
  TH1D* h_cutflow;
  std::vector<int> m_cutflow;
  std::vector<TString> m_cutNames;

  // Regular
  TH1D* h_Mff;
  TH1D* h_ytt;
  TH1D* h_Pz_nu;
  TH1D* h_CosTheta;
  TH1D* h_CosThetaStar;
  
  // Reconstruction histograms
  TH1D* h_Pz_nu_r;
  TH1D* h_Mtt_r;
  TH1D* h_CosTheta_r;
  TH1D* h_CosThetaStar_r;
  TH1D* h_ytt_r;

  // Polarisation weighted histograms
  TH1D* h_MttLL;
  TH1D* h_MttLR;
  TH1D* h_MttRL;
  TH1D* h_MttRR;

  // Spin asymmetry historgrams
  TH1D* h_ALL;
  TH1D* h_AL;

  // Counting
  TH1D* h_imaginary;
  TH1D* h_real;
  TH1D* h_imaginary_r1;
  TH1D* h_real_r1;
  TH1D* h_imaginary_r2;
  TH1D* h_real_r2;

  // Charge asymmetry histograms
  TH1D* h_AFBstar;
  TH1D* h_AFBstarF;
  TH1D* h_AFBstarB;

  TH1D* h_AFBstar_r;
  TH1D* h_AFBstar_rF;
  TH1D* h_AFBstar_rB;

  TH1D* h_AttC;
  TH1D* h_AttCF;
  TH1D* h_AttCB;

  TH1D* h_AllC;  
  TH1D* h_AllCF;
  TH1D* h_AllCB;

  TH1D* h_AlL;
  TH1D* h_AlLF;
  TH1D* h_AlLB;

  // Final particle 4-vectors
  vector<TLorentzVector> p;
  vector<TLorentzVector> pcm;
  vector<TLorentzVector> p_r1;
  vector<TLorentzVector> p_r2;
  vector<TLorentzVector> pcm_r1;
  vector<TLorentzVector> pcm_r2;

  // Event 4-vectors
  TLorentzVector P;
  TLorentzVector Pcm;
  TLorentzVector P_r1;
  TLorentzVector P_r2;

  // Top 4-vectors
  TLorentzVector p_t;
  TLorentzVector p_tb;
  TLorentzVector p_t_r1;
  TLorentzVector p_tb_r1;
  TLorentzVector p_t_r2;
  TLorentzVector p_tb_r2;

  typedef vector<TString>::const_iterator Itr_s;
};
#endif
