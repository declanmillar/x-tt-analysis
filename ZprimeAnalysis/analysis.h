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

class AnalysisZprime{
public:
  AnalysisZprime(const TString channel, const TString model, const TString& outputFileName);
  virtual ~AnalysisZprime();
  
protected:
  // Internal for Event Looping
  Long64_t TotalEvents();
  Long64_t IncrementEvent(Long64_t i);
  void SetupTreesForNewFile(const TString& s);
  void CleanUp();
  
  void SetupInputFiles();
  
  void PreLoop();
  void Loop();
  void EachEvent();
  void PostLoop();
  void MakeGraphs();
  void TotalSpinAsymmetries();
  TH1D* Asymmetry(TH1D* h_A, TH1D* h_B);
  void ALL2to6();
  TH1D* MttALL();
  TH1D* MttAL();
  
  bool PassCuts() const;
  bool PassCuts_MET() const;
  bool PassCuts_Mtt() const;
    
  // Delta R
  // float deltaR(const float& eta1,const float& eta2,const float& phi1,const float& phi2) const;
  // float deltaPhi(const float& phi1,const float& phi2) const;    
  
  
private:
  AnalysisZprime();
  AnalysisZprime(const AnalysisZprime& rhs);  
  void operator=(const AnalysisZprime& rhs);  
  
  float m_pi;
  float m_GeV;
  TString m_channel;  
  TString m_model;
  TString m_outputFileName;

  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;
  
  // Define Histograms
  TH1D* h_pz5;
  TH1D* h_costheta5_eq;
  TH1D* h_costheta5_ee;
  TH1D* h_ct7ct5;
  TH1D* h_AFstar;
  TH1D* h_ABstar;
  TH1D* h_AFBstar;
  TH1D* h_Mtt;
  TH1D* h_MttLL;
  TH1D* h_MttLR;
  TH1D* h_MttRL;
  TH1D* h_MttRR;
  TH1D* h_MttALL;
  TH1D* h_MttAL;
  TH1D* h_MET;
  TH1D* h_HT;
  TH1D* h_Mbbll;
  TH1D* h_mll;
  TH1D* h_ETbbll;
  TH1D* h_KTbbll;
  TH1D* h_MTll;
  TH1D* h_MCTll;
  TH1D* h_MTblbl;
  TH1D* h_MCTblbl;
  typedef vector<TString>::const_iterator Itr_s;
};
#endif
