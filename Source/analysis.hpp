#ifndef _ANALYSIS_ZPRIME_H_
#define _ANALYSIS_ZPRIME_H_

#include "RootTuple.hpp"
#include "atlas_style.hpp"
#include "TCanvas.h"
#include "TApplication.h"
#include <cmath>
#include <math.h>
#include <TString.h>
#include <TH2.h>
#include <TF1.h>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TBrowser.h"
#include <fstream>
#include "TLegend.h"
#include <complex>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h>
#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>

class AnalysisZprime{
public:
  AnalysisZprime(const TString, const TString, const TString, const TString, const int energy, const TString options, const int vegasIterations, const string vegasPoints, const bool add_ggG, const bool add_qqG, const int luminosity, const int btags, const bool discardComplex, const TString analysis_label);
  virtual ~AnalysisZprime();
  TString GetOutputFilename();
  void SetYttCut(const double);
  void SetXsec(const bool);
  void SetFiducial(const bool);
  void Run();

protected:
  Long64_t TotalEvents();
  Long64_t IncrementEvent(Long64_t i);
  void SetupTreesForNewFile(const TString& s);
  void CleanUp();

  void SetupInputFiles();
  void SetupOutputFiles();

  void PreLoop();
  void Loop();
  void PostLoop();
  void EachEvent();
  void CreateHistograms();
  void MakeGraphs();
  void MakeDistribution(TH1D* h, TString units);
  void Make2dDistribution(TH2D*);
  void NormalizeSliceY(TH2D*);
  void WriteHistograms();
  void CheckResults();
  void CheckPerformance();
  void CreateFilenames();
  void CheckFiles();
  void GetCrossSection(TString log);
  void GetIterationWeights(TString log);

  void GetDataDirectory();
  void GetChannelFactors();

  double TotalAsymmetry(TH1D* h_A, TH1D* h_B);
  void TotalSpinAsymmetries();
  TH1D* Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B);
  TH1D* Asymmetry2(TString name, TString title, TH1D* h_A, TH1D* h_B);
  void ApplyLuminosity(TH1D*);
  void AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B);
  TH1D* MakeALL();
  TH1D* MakeAL();
  vector<std::complex<double> > SolveQuadratic(double a, double b, double c);
  std::vector<TLorentzVector> ReconstructSemiLeptonic(std::vector<TLorentzVector> p, int l_Q);
  std::vector<TLorentzVector> ReconstructDilepton(std::vector<TLorentzVector> p);

  bool PassCuts(string type);
  bool PassCutsMET(string type);
  bool PassCutsMtt(string type);
  bool PassCutsEta(string type);
  bool PassCutsYtt(string type);
  bool PassCutsET(string type);

  void ResetCounters();
  void InitialiseCutflow();
  void PrintCutflow();
  const void UpdateCutflow(int cut, bool passed);

  static inline void ProgressBar(unsigned int x, unsigned int n, unsigned int w);

private:
  AnalysisZprime();
  AnalysisZprime(const AnalysisZprime& rhs);
  void operator = (const AnalysisZprime& rhs);

  typedef vector<TString>::const_iterator Itr_s;

  // arguments
  TString m_model;
  TString m_initial_state;
  TString m_intermediates;
  TString m_channel;
  int m_energy;
  TString m_options;
  int m_vegasIterations;
  string m_vegasPoints;
  int m_add_ggG;
  int m_add_qqG;
  int m_luminosity;
  double m_efficiency;
  const int nBtags;
  const bool m_discardComplex;
  const TString m_analysisLabel;
  bool m_xsec;
  const int m_reco;
  bool m_fid;

  // Parameters
  float m_pi;
  float m_GeV;
  double m_bmass;
  double m_Wmass;
  double m_tmass;

  // Cuts
  double m_ytt;
  double m_Emin;
  double m_Emax;
  enum m_cutlist{
    c_entries,
    c_topDecays,
    c_antitopDecays,
    c_events,
    c_realSolutions,
    c_Et,
    c_eta,
    c_MET,
    c_mtt,
    c_ytt,
    m_cuts // Keep as last entry
  };

  // others
  bool m_useLumi;
  bool m_discardEvent;

  // Counter
  unsigned int m_nReco;
  unsigned int m_nQuarksMatched;
  unsigned int m_nNeutrinoMatched;
  unsigned int m_nRealRoots;
  unsigned int m_nComplexRoots;
  bool m_R1solutionIsReal;
  bool m_R2solutionIsReal;

  double m_sigma;
  double channelFactor;
  double channelFactor_R;
  vector<double> iteration_weights;

  // Strings
  TString m_inputFileName;
  TString m_QCDfilename;
  TString m_QCDweightFile;
  TString m_weightsFilename;
  TString m_outputFilename;
  TString m_dataDirectory;

  // Input data
  vector<TString>* m_inputFiles;
  vector<TString>* m_weightFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup;

  // OutputFile
  TFile* m_outputFile;

  // Cutflow
  std::vector<int> m_cutflow;
  std::vector<TString> m_cutNames;

  // Final particle 4-vectors
  vector<TLorentzVector> p;
  vector<TLorentzVector> pcm;
  vector<TLorentzVector> ptop;
  vector<TLorentzVector> patop;
  vector<TLorentzVector> p_R1;
  vector<TLorentzVector> p_R2;
  vector<TLorentzVector> pcm_R1;
  vector<TLorentzVector> pcm_R2;
  vector<TLorentzVector> ptop_R1;
  vector<TLorentzVector> patop_R2;

  // Event 4-vectors
  TLorentzVector P;
  TLorentzVector Pcm;
  TLorentzVector Ptop;
  TLorentzVector Patop;
  TLorentzVector P_R1;
  TLorentzVector P_R2;

  // Top 4-vectors
  TLorentzVector pcm_t;
  TLorentzVector pcm_tb;
  TLorentzVector p_t_R1;
  TLorentzVector p_tb_R1;
  TLorentzVector p_t_R2;
  TLorentzVector p_tb_R2;

  // Histograms

  // naming conventions
  // p = momenta
  // P = sum of momenta in frame
  // c = cos
  // m = mass
  // A = asymmetry
  // B = backward (in costheta*)
  // F = forward (in costheta*)
  // By = backward (in delta_y)
  // Fy = forward (in delta_y)
  // Bly = backward (in delta_y for decay lepton)
  // Fly = forward (in delta_y for decay lepton)
  // C = complex only (reconstructed)
  // D = discard complex (reconstructed)
  // R = reconstructed
  // S = star [reconstructed parton centre of mass frame (now assumed in AFB, as original not useful here)]
  // LL, RR, LR, RL = helicities
  // n = normalised, e.g. mtt_F/mtt

  // Masses
  TH1D* h_mtt;
  TH1D* h_mt;
  TH1D* h_mtbar;
  TH1D* h_mtt_F;
  TH1D* h_mtt_B;
  TH1D* h_mtt_Fn;
  TH1D* h_mtt_Bn;
  TH1D* h_mtt_Fy;
  TH1D* h_mtt_By;

  TH1D* h_mtt_R;
  TH1D* h_mt_R;
  TH1D* h_mtbar_R;
  TH1D* h_mtt_FR;
  TH1D* h_mtt_BR;
  TH1D* h_mtt_FRn;
  TH1D* h_mtt_BRn;
  TH1D* h_mtt_FD;
  TH1D* h_mtt_BD;
  TH1D* h_mtt_Fl;
  TH1D* h_mtt_Bl;

  TH1D* h_mtt_LL;
  TH1D* h_mtt_LR;
  TH1D* h_mtt_RL;
  TH1D* h_mtt_RR;

  // Momenta
  TH1D* h_pv1x;
  TH1D* h_pv1x_R;
  TH1D* h_pv1y;
  TH1D* h_pv1y_R;
  TH1D* h_pv1z;
  TH1D* h_pv1z_R;
  TH1D* h_pv2x;
  TH1D* h_pv2x_R;
  TH1D* h_pv2y;
  TH1D* h_pv2y_R;
  TH1D* h_pv2z;
  TH1D* h_pv2z_R;

  // Rapidities
  TH1D* h_ytt;
  TH1D* h_ytt_R;

  // Angles
  TH1D* h_deltaPhi;
  TH1D* h_cosTheta;
  TH1D* h_cosTheta1;
  TH1D* h_cosTheta2;
  TH1D* h_cos1cos2;
  TH1D* h_cosThetaStar;
  TH1D* h_cosTheta_R;
  TH1D* h_cosThetaStar_R;

  // Counting
  TH1D* h_cutflow;

  // Forward Backward Asymmetries
  TH1D* h_AFB;
  TH1D* h_AFB_R;

  // Spin asymmetries
  TH1D* h_Ap;

  // transverse
  TH1D* h_HT;
  TH1D* h_KT;
  TH1D* h_AL1;
  TH1D* h_AL2;
  TH1D* h_AL_R;

  // 2D histograms
  TH2D* h2_mtt_cosThetaStar;
  TH2D* h2_mtt_cosThetaStar_R;

  TH2D* h2_mtt_deltaPhi;
  TH2D* h2_mtt_cosTheta1;
  TH2D* h2_mtt_cosTheta2;
  TH2D* h2_mtt_cosThetal_R;
  TH2D* h2_mtt_cos1cos2;

  TH2D* h2_HT_deltaPhi;
  TH2D* h2_KT_deltaPhi;

  //
  vector<TH1D*> h_eta;
  vector<TH1D*> h_pt;
  vector<TH1D*> h_deltaRs;
  TH1D* h_deltaRmax;
  TH1D* h_deltaRbW;
};
#endif
