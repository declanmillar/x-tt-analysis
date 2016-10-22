#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <fstream>
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "boost/algorithm/string.hpp"
#include "boost/asio.hpp"
#include "root-tuple.hpp"
#include "atlas-style.hpp"
#include "solve-poly.hpp"


class Analysis{
    Analysis();
    Analysis(const Analysis& rhs);
    void operator = (const Analysis& rhs);

    typedef std::vector<TString>::const_iterator Itr_s;

    TString m_model;
    TString m_initial_state;
    TString m_intermediates;
    TString m_channel;
    int m_energy;
    TString m_options;
    int m_vegasIterations;
    std::string m_vegasPoints;
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
    float m_pi;
    float m_GeV;
    double m_bmass;
    double m_Wmass;
    double m_tmass;
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
    bool m_useLumi;
    bool m_discardEvent;
    unsigned int m_nReco;
    unsigned int m_nQuarksMatched;
    unsigned int m_nNeutrinoMatched;
    unsigned int m_nRealRoots;
    unsigned int m_nComplexRoots;
    bool m_R1solutionIsReal;
    bool m_R2solutionIsReal;
    double m_sigma;
    std::vector<double> iteration_weights;
    TString m_dataDirectory;
    TString m_outputFilename;
    std::vector<TString>* m_inputFiles;
    std::vector<TString>* m_weightFiles;
    RootTuple* m_ntup;
    TChain* m_chainNtup;
    TFile* m_outputFile;
    std::vector<int> m_cutflow;
    std::vector<TString> m_cutNames;
    const bool m_debug;

    // 4-vectors
    std::vector<TLorentzVector> p;
    std::vector<TLorentzVector> pcm;
    std::vector<TLorentzVector> ptop;
    std::vector<TLorentzVector> patop;
    std::vector<TLorentzVector> p_R1;
    std::vector<TLorentzVector> p_R2;
    std::vector<TLorentzVector> pcm_R1;
    std::vector<TLorentzVector> pcm_R2;
    std::vector<TLorentzVector> ptop_R1;
    std::vector<TLorentzVector> patop_R2;
    TLorentzVector P;
    TLorentzVector P_R1;
    TLorentzVector P_R2;

    // Histograms
    TH1D* h_mtt;
    TH1D* h_mt;
    TH1D* h_mtbar;
    TH1D* h_mtt_F;
    TH1D* h_mtt_B;
    TH1D* h_mtt_Fn;
    TH1D* h_mtt_Bn;
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
    TH1D* h_ytt;
    TH1D* h_ytt_R;
    TH1D* h_deltaPhi;
    TH1D* h_cosTheta;
    TH1D* h_cosTheta1;
    TH1D* h_cosTheta2;
    TH1D* h_cos1cos2;
    TH1D* h_cosThetaStar;
    TH1D* h_cosTheta_R;
    TH1D* h_cosThetaStar_R;
    TH1D* h_cutflow;
    TH1D* h_AFB;
    TH1D* h_AFB_R;
    TH1D* h_Ap;
    TH1D* h_AL1;
    TH1D* h_AL2;
    TH1D* h_AL_R;
    TH1D* h_HT;
    TH1D* h_KT;
    TH2D* h2_mtt_cosThetaStar;
    TH2D* h2_mtt_cosThetaStar_R;
    TH2D* h2_mtt_deltaPhi;
    TH2D* h2_mtt_cosTheta1;
    TH2D* h2_mtt_cosTheta2;
    TH2D* h2_mtt_cosThetal_R;
    TH2D* h2_mtt_cos1cos2;
    TH2D* h2_HT_deltaPhi;
    TH2D* h2_KT_deltaPhi;
    TH1D* h_deltaRmax;
    TH1D* h_deltaRbW;
    std::vector<TH1D*> h_eta;
    std::vector<TH1D*> h_pt;
    std::vector<TH1D*> h_deltaRs;
  protected:
    Long64_t TotalEvents();
    Long64_t IncrementEvent(Long64_t i);
    void SetupTreesForNewFile(const TString&);
    void CleanUp();
    void SetupInputFiles();
    void SetupOutputFiles();
    void PreLoop();
    void Loop();
    void PostLoop();
    void EachEvent();
    void MakeHistograms();
    void MakeDistributions();
    void MakeDistribution1D(TH1D*, TString);
    void MakeDistribution2D(TH2D*);
    void NormalizeSliceY(TH2D*);
    void WriteHistograms();
    void CheckResults();
    void CheckPerformance();
    void CreateFilenames();
    void CheckFiles();
    void GetCrossSection(TString);
    void GetIterationWeights(TString);
    void SetDataDirectory();
    void GetChannelFactors();
    void ApplyLuminosity(TH1D*);
    void AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B);
    void ResetCounters();
    void InitialiseCutflow();
    void PrintCutflow();
    const void UpdateCutflow(int cut, bool passed);  
    bool PassCuts(string type);
    bool PassCutsMET(string type);
    bool PassCutsMtt(string type);
    bool PassCutsEta(string type);
    bool PassCutsYtt(string type);
    bool PassCutsET(string type);
    TH1D* Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B);
    std::vector<TLorentzVector> ReconstructSemilepton(std::vector<TLorentzVector> p, int l_Q);
    std::vector<TLorentzVector> ReconstructDilepton(std::vector<TLorentzVector> p);
  public:
    Analysis(const TString, const TString, const TString, const TString, const int energy, const TString options, const int vegasIterations, const string vegasPoints, const bool add_ggG, const bool add_qqG, const int luminosity, const int btags, const bool discardComplex, const TString analysis_label);
    virtual ~Analysis();
    TString GetOutputFilename();
    void SetYttCut(const double);
    void SetXsec(const bool);
    void SetFiducial(const bool);
    void Run();
};
#endif
