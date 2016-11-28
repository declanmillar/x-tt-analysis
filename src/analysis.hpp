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
    TString m_process;
    TString m_options;
    bool m_add_gg;
    int m_add_qq;
    int m_energy;
    int m_luminosity;
    TString m_tag;

    bool m_xsec = false;
    int m_reco;
    bool m_fid;
    double m_ytt = 0;
    double m_Emin = -1;
    double m_Emax = -1;
    bool m_useLumi;

    int n;
    unsigned int m_nReco;
    unsigned int m_nQuarksMatched;
    unsigned int m_nNeutrinoMatched;
    unsigned int m_nRealRoots;
    unsigned int m_nComplexRoots;
    double m_sigma;
    Long64_t m_nevents;

    std::vector<double> iteration_weights;
    TString m_dataDirectory;
    TString m_outputFilename;
    TFile* m_outputFile;
    std::vector<TString>* m_inputFiles;
    std::vector<TString>* m_weightFiles;
    TChain* m_chainNtup;
    RootTuple* m_ntup;

    std::vector<int> m_cutflow;
    std::vector<TString> m_cutNames;
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

    const double m_pi = 3.14159265358979323846;
    // const double m_bmass = 4.18, m_Wmass = 80.23, m_tmass = 173.0;
    const double m_bmass = 4.18, m_Wmass = 80.4, m_tmass = 172.5;
    const int m_btags = 2;
    const double m_efficiency = 1.0;

    // Histograms
    TH1D* h_mtt;
    TH1D* h_mtt_R;

    TH1D* h_mtt_tF;
    TH1D* h_mtt_tF_R;
    TH1D* h_mtt_tB;
    TH1D* h_mtt_tB_R;

    TH1D* h_mtt_tlF;
    TH1D* h_mtt_tlF_R;
    TH1D* h_mtt_tlB;
    TH1D* h_mtt_tlB_R;

    TH1D* h_mtt_tCF;
    TH1D* h_mtt_tCF_R;
    TH1D* h_mtt_tCB;
    TH1D* h_mtt_tCB_R;

    TH1D* h_mtt_lF;
    TH1D* h_mtt_lB;

    TH1D* h_mtt_philF;
    TH1D* h_mtt_philB;

    TH1D* h_mtt_ElF;
    TH1D* h_mtt_ElB;

    TH1D* h_mtt_LL;
    TH1D* h_mtt_LR;
    TH1D* h_mtt_RL;
    TH1D* h_mtt_RR;

    // t
    TH1D* h_pxt;
    TH1D* h_pxt_R;
    TH1D* h_pyt;
    TH1D* h_pyt_R;
    TH1D* h_pzt;
    TH1D* h_pzt_R;
    TH1D* h_Et;
    TH1D* h_Et_R;
    TH1D* h_pTt;
    TH1D* h_pTt_R;
    TH1D* h_etat;
    TH1D* h_etat_R;
    TH1D* h_phit;
    TH1D* h_phit_R;
    TH1D* h_mt;
    TH1D* h_mt_R;

    // tbar
    TH1D* h_pxtbar;
    TH1D* h_pxtbar_R;
    TH1D* h_pytbar;
    TH1D* h_pytbar_R;
    TH1D* h_pztbar;
    TH1D* h_pztbar_R;
    TH1D* h_Etbar;
    TH1D* h_Etbar_R;
    TH1D* h_pTtbar;
    TH1D* h_pTtbar_R;
    TH1D* h_etatbar;
    TH1D* h_etatbar_R;
    TH1D* h_phitbar;
    TH1D* h_phitbar_R;
    TH1D* h_mtbar;
    TH1D* h_mtbar_R;

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
    TH1D* h_costheta_tt;
    TH1D* h_costheta_tt_R;
    TH1D* h_cutflow;
    TH1D* h_HT;
    TH1D* h_KT;

    // asymmetries
    TH1D* h_AtFB;
    TH1D* h_AtFB_R;
    TH1D* h_AtC;
    TH1D* h_AtC_R;
    TH1D* h_Ap;
    TH1D* h_AL1;
    TH1D* h_AL2;
    TH1D* h_AL_R;
    TH1D* h_AtlFB;
    TH1D* h_AtlFB_R;
    TH1D* h_Aphil;
    TH1D* h_AlEl;
    TH1D* h_AL;
    TH1D* h_ALL;

    // performance
    TH1D* h_pT_t_perf;
    TH1D* h_pT_tbar_perf;
    TH1D* h_costheta_tt_perf;
    TH1D* h_m_tt_perf;
    TH2D* h2_m_tt_pT_t_perf;
    TH2D* h2_m_tt_costheta_tt_perf;

    TH1D* h_deltaR_max;
    TH1D* h_deltaR_bW;
    TH1D* h_deltaR_tt;

    std::vector<TH1D*> h_eta;
    std::vector<TH1D*> h_pt;
    std::vector<TH1D*> h_deltaRs;

    TH2D* h2_mtt_cosThetaStar;
    TH2D* h2_mtt_cosThetaStar_R;
    TH2D* h2_mtt_deltaPhi;
    TH2D* h2_mtt_cosTheta1;
    TH2D* h2_mtt_cosTheta2;
    TH2D* h2_mtt_cosThetal_R;
    TH2D* h2_mtt_cos1cos2;
    TH2D* h2_HT_deltaPhi;
    TH2D* h2_KT_deltaPhi;
    
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
    void MakeDistribution1D(TH1D*, const TString&);
    void MakeDistribution2D(TH2D*, TString, TString, TString, TString);
    void NormalizeSliceY(TH2D*);
    void WriteHistograms();
    void CheckResults();
    void CheckPerformance();
    void CreateFilenames();
    void CheckFiles();
    void GetGenerationCrossSection(TString);
    void SetDataDirectory();
    void GetChannelFactors();
    void AsymmetryUncertainty(TH1D*, TH1D*, TH1D*);
    void ResetCounters();

    TH1D* MakeALL();
    TH1D* MakeAL();
    void TotalSpinAsymmetries();
    double TotalAsymmetry(TH1D* h_A, TH1D* h_B);

    void InitialiseCutflow();
    void PrintCutflow();
    void UpdateCutflow(int, bool);
    bool PassCuts(const std::vector<TLorentzVector>&, const TLorentzVector&);
    bool PassCutsMET(const std::vector<TLorentzVector>&, const TLorentzVector&);
    bool PassCutsMtt(const std::vector<TLorentzVector>&, const TLorentzVector&);
    bool PassCutsEta(const std::vector<TLorentzVector>&, const TLorentzVector&);
    bool PassCutsYtt(const std::vector<TLorentzVector>&, const TLorentzVector&);
    bool PassCutsET(const std::vector<TLorentzVector>&, const TLorentzVector&);

    TH1D* Asymmetry(const TString&, const TString&, TH1D*, TH1D*);
    std::vector<TLorentzVector> ReconstructSemilepton(const std::vector<TLorentzVector>&, const int);
    std::vector<TLorentzVector> ReconstructDilepton(const std::vector<TLorentzVector>&);
  public:
    Analysis(const TString& model, const TString& process, const TString& options, const bool gg, const bool qq, const int energy, const int luminosity, const TString& tag);
    virtual ~Analysis();
};
#endif
