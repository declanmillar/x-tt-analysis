#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <fstream>
#include <tuple>
#include "TH2.h"
#include "TF1.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TClonesArray.h"
#include "sys/stat.h"
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/asio.hpp"
#include "atlas-style.hpp"
#include "solve-poly.hpp"
#include "neutrino-weighting.hpp"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

class Analysis {

public:
    Analysis( const std::string& model, const std::string& process, const std::string& options, const int energy, const int luminosity, const std::string& reconstruction, const std::string& tag );
    virtual ~Analysis();
    void Run();

protected:
    void SetupTreesForNewFile( const std::string& );
    void CleanUp();
    void SetupInputFiles();
    void SetupOutputFiles();

    void PreLoop();
    void Loop();
    void PostLoop();

    void EachEvent( double );
    void EveryEvent( double );
    void CleanupEvent();
    void AssignChannel();

    // histograms
    void MakeHistograms();
    void MakeDistributions();
    void MakeDistribution1D( TH1D*, const std::string& );
    void MakeDistribution2D( TH2D*, std::string, std::string, std::string, std::string );
    void NormalizeSliceY( TH2D* );
    void WriteHistograms();

    void CheckResults();
    void GetGenerationCrossSection( int );
    void GetProcessWeight( int );
    void SetDataDirectory();
    void GetChannelFactors();
    void AsymmetryUncertainty( TH1D*, TH1D*, TH1D* );
    void ResetCounters();
    void GetBranches();
    void EachFile( const std::string& );

    // cutflow
    void InitialiseCutflow();
    void PrintCutflow();
    void UpdateCutflow( int, bool );

    // event selection
    bool PassesEventSelection();
    bool ExactlyTwoLeptons();
    bool OppositeCharge();
    bool SufficientJets();
    bool SufficientBtags();
    bool SufficientHT();
    bool SufficientMET();
    bool SufficientMll(const  std::pair< TLorentzVector, TLorentzVector >& );
    bool OutsideZmassWindow(const  std::pair< TLorentzVector, TLorentzVector >& );

    Long64_t TotalEvents();
    Long64_t IncrementEvent( Long64_t i );
    double TotalAsymmetry( TH1D* h_A, TH1D* h_B );
    TH1D* Asymmetry( const std::string&, const std::string&, TH1D*, TH1D* );

    // tuple
    TClonesArray* b_Jet;
    TClonesArray* b_Electron;
    TClonesArray* b_Muon;
    TClonesArray* b_MissingET;
    TClonesArray* b_ScalarHT;


private:
    Analysis();
    Analysis( const Analysis& rhs );
    void operator = ( const Analysis& rhs );

    typedef std::vector< std::tuple< std::string, int > >::const_iterator itr_s;
    std::vector< std::tuple< std::string, int > >* m_input;
    std::vector< std::tuple< std::string, int, int, double, double, double > >* m_processes;

    std::vector< Electron* >* m_electron;
    std::vector< Muon* >* m_muon;
    std::vector< Jet* >* m_jet;

    TString m_model;
    std::string m_process;
    TString m_options;
    int m_energy;
    int m_luminosity;
    TString m_tag;
    std::string m_pdf = "CT14LL";

    bool m_xsec = false;
    bool m_fid = false;
    bool m_iso = false;
    const std::string m_reconstruction;
    double m_ytt = 0.0;
    double m_Emin = -1;
    double m_Emax = -1;
    bool m_useLumi;
    const bool m_debug = false;
    std::string m_channel;

    unsigned int m_nReco;
    unsigned int m_nQuarksMatched;
    unsigned int m_nNeutrinoMatched;
    unsigned int m_nRealRoots;
    unsigned int m_nComplexRoots;
    double m_crossSection;
    Long64_t m_nevents;

    std::vector<double> iteration_weights;
    std::string m_dataDirectory;
    std::string m_outputName;
    TFile* m_output;
    TChain* m_chain;
    ExRootTreeReader* m_tree;

    std::vector< int > m_cutflow;
    std::vector< std::string > m_cutNames;
    enum m_cutlist{
        c_events,
        c_sufficientMET,
        c_sufficientHT,
        c_twoLeptons,
        c_oppositeCharge,
        c_sufficientMll,
        c_outsideZmassWindow,
        c_sufficientJets,
        c_sufficientBtags,
        c_realSolutions,
        c_deltaR,
        m_cuts // Keep as last entry
    };

    const double m_pi = 3.14159265358979323846;
    const double m_bmass = 4.18, m_Wmass = 80.4, m_zmass = 91.19, m_tmass = 172.5;
    const int m_minimumBtags = 2;
    int m_btags;
    const double m_efficiency = 1.0;


    const bool truth = false;
    // Histograms
    TH1D* h_pt_l1;
    TH1D* h_eta_l1;
    TH1D* h_pt_l2;
    TH1D* h_eta_l2;

    TH1D* h_pt_jets;
    TH1D* h_eta_jets;
    TH1D* h_pt_bjets;
    TH1D* h_eta_bjets;
    TH1D* h_pt_qjets;
    TH1D* h_eta_qjets;

    TH1D* h_pt_alljets;
    TH1D* h_pt_allel;
    TH1D* h_pt_allmu;
    TH1D* h_eta_alljets;
    TH1D* h_eta_allel;
    TH1D* h_eta_allmu;

    TH1D* h_HT;
    TH1D* h_KT;
    TH1D* h_mvis;
    TH1D* h_HT_all;
    TH1D* h_KT_all;
    TH1D* h_mvis_all;

    TH1D* h_deltaPhi;

    TH1D* h_mW1;
    TH1D* h_mW2;

    // t
    TH1D* h_pxt;
    TH1D* h_pyt;
    TH1D* h_pzt;
    TH1D* h_Et;
    TH1D* h_pTt;
    TH1D* h_etat;
    TH1D* h_phit;
    TH1D* h_mt;

    // tbar
    TH1D* h_pxtbar;
    TH1D* h_pytbar;
    TH1D* h_pztbar;
    TH1D* h_Etbar;
    TH1D* h_pTtbar;
    TH1D* h_etatbar;
    TH1D* h_phitbar;
    TH1D* h_mtbar;

    // tt
    TH1D* h_mtt;
    TH1D* h_mtt_tF;
    TH1D* h_mtt_tB;
    TH1D* h_mtt_tlF;
    TH1D* h_mtt_tlB;
    TH1D* h_mtt_tCF;
    TH1D* h_mtt_tCB;
    TH1D* h_mtt_lF;
    TH1D* h_mtt_lB;
    TH1D* h_mtt_philF;
    TH1D* h_mtt_philB;
    TH1D* h_mtt_ElF;
    TH1D* h_mtt_ElB;
    TH1D* h_ytt;

    TH1D* h_pv1x;
    TH1D* h_pv1y;
    TH1D* h_pv1z;
    TH1D* h_pv2x;
    TH1D* h_pv2y;
    TH1D* h_pv2z;
    TH1D* h_cosTheta;
    TH1D* h_cosTheta1;
    TH1D* h_cosTheta2;
    TH1D* h_cos1cos2;
    TH1D* h_cosThetaStar;
    TH1D* h_costheta_tt;

    TH1D* h_cutflow;

    // asymmetries
    TH1D* h_AtFB;
    TH1D* h_AtC;
    TH1D* h_Ap;
    TH1D* h_AL1;
    TH1D* h_AL2;
    TH1D* h_AtlFB;
    TH1D* h_Aphil;
    TH1D* h_AlEl;
    TH1D* h_AL;
    TH1D* h_ALL;

    TH1D* h_deltaR_max;
    TH1D* h_deltaR_bW;
    TH1D* h_deltaR_tt;

    TH1D* h_nElectrons;
    TH1D* h_nMuons;
    TH1D* h_nJets;

    TH2D* h2_mtt_cosThetaStar;
    TH2D* h2_mtt_deltaPhi;
    TH2D* h2_mtt_cosTheta1;
    TH2D* h2_mtt_cosTheta2;
    TH2D* h2_mtt_cos1cos2;
    TH2D* h2_HT_deltaPhi;
    TH2D* h2_mvis_deltaPhi;
    TH2D* h2_KT_deltaPhi;
};
#endif
