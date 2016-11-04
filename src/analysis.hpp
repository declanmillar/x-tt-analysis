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
    int m_add_gg;
    int m_add_qq;
    int m_luminosity;
    double m_efficiency;
    const int m_btags;
    const TString m_tag;

    bool m_xsec;
    const int m_reco;
    bool m_fid;
    double m_ytt = 0;
    double m_Emin = -1;
    double m_Emax = -1;
    bool m_useLumi;

    unsigned int m_nReco;
    unsigned int m_nQuarksMatched;
    unsigned int m_nNeutrinoMatched;
    unsigned int m_nRealRoots;
    unsigned int m_nComplexRoots;
    double m_sigma;

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
    const double m_bmass = 4.18, m_Wmass = 80.23, m_tmass = 173.0;

    // Histograms
    TH1D* h_mtt;
    TH1D* h_mtt_F;
    TH1D* h_mtt_B;
    TH1D* h_mtt_Fn;
    TH1D* h_mtt_Bn;
    TH1D* h_mtt_R;
    TH1D* h_mtt_FR;
    TH1D* h_mtt_BR;
    TH1D* h_mtt_FRn;
    TH1D* h_mtt_BRn;
    TH1D* h_mtt_FD;
    TH1D* h_mtt_BD;
    TH1D* h_mtt_Fl;
    TH1D* h_mtt_Bl;

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
    TH1D* h_costhetatt;
    TH1D* h_costhetatt_R;
    TH1D* h_cutflow;
    TH1D* h_AFB;
    TH1D* h_AFB_R;
    TH1D* h_Ap;
    TH1D* h_AL1;
    TH1D* h_AL2;
    TH1D* h_AL_R;
    TH1D* h_HT;
    TH1D* h_KT;

    TH1D* h_deltaRmax;
    TH1D* h_deltaRbW;

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
    void AsymmetryUncertainty(TH1D*, TH1D*, TH1D*);
    void ResetCounters();

    void InitialiseCutflow();
    void PrintCutflow();
    void UpdateCutflow(int, bool);
    bool PassCuts(const string&);
    bool PassCutsMET(const string&);
    bool PassCutsMtt(const string&);
    bool PassCutsEta(const string&);
    bool PassCutsYtt(const string&);
    bool PassCutsET(const string&);

    TH1D* Asymmetry(const TString&, const TString&, TH1D*, TH1D*);
    std::vector<TLorentzVector> ReconstructSemilepton(const std::vector<TLorentzVector>&, const int);
    std::vector<TLorentzVector> ReconstructDilepton(const std::vector<TLorentzVector>&);
  public:
    Analysis(const TString&, const TString&, const TString&, const TString&, const int, const TString&, const int, const string, const bool, const bool, const int, const int, const TString&);
    virtual ~Analysis();
    TString GetOutputFilename();
    void SetEnergyRange(double, double);
    void SetYttCut(const double);
    void SetXsec(const bool);
    void SetFiducial(const bool);
    void Run();
};
#endif
