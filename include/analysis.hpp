#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <fstream>
#include <tuple>
#include <float.h>
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/asio.hpp"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TVector.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "sys/stat.h"
#include "TTF.h"
#include "TCanvas.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

class Analysis {

private:
    Analysis();
    Analysis(const Analysis& rhs);
    void operator = (const Analysis& rhs);

     // args: inputfilename, proc_id
    vector<tuple<string, int>>* m_input;
    typedef vector<tuple<string, int>>::const_iterator itr_s;
    
    // args: processfilename, proc_id, n_proc, cross_section, uncertainty, weight
    vector<tuple<string, int, int, double, double, double>>* m_processes;

    vector<Electron*>* m_electrons;
    vector<GenParticle*>* m_truthElectrons;
    vector<Muon*>* m_muons;
    vector<GenParticle*>* m_truthMuons;
    vector<Jet*>* m_jets;
    vector<GenParticle*>* m_truthBquarks;
    // 
    // vector<bool>* m_electron_truth_tags;
    // vector<bool>* m_muon_truth_tags;

    GenParticle* m_hardTop;
    GenParticle* m_hardTbar;
    GenParticle* m_hardB;
    GenParticle* m_hardBbar;
    GenParticle* m_hardLepP;
    GenParticle* m_hardLepM;
    GenParticle* m_hardNu;
    GenParticle* m_hardNuBar;
    
    vector<bool>* m_lepton_truth_tags;
    vector<bool>* m_electron_truth_tags;
    vector<bool>* m_muon_truth_tags;
    vector<bool>* m_jet_truth_tags;

    string m_inputFileName;
    string m_processfilename;
    string m_model;
    string m_process;
    string m_options;
    int m_energy;
    int m_luminosity;
    string m_tag;
    string m_pdf = "CT14LL";
    bool m_useMassSlices = false;

    bool m_xSec;
    const string m_reconstruction;
    bool m_useLumi;
    const bool m_debug;
    string m_channel;

    double m_crossSection;
    Long64_t m_nevents;

    vector<double> iteration_weights;
    string m_dataDirectory;
    string m_outputFileName;
    TFile* m_outputFile;
    TChain* m_chain;
    ExRootTreeReader* m_tree;

    vector<int> m_cutflow;
    vector<string> m_cutNames;
    vector<string> m_cutTitles;
    enum m_cutlist{
        c_events,
        c_twoLeptons,
        c_oppositeCharge,
        c_sufficientMll,
        c_outsideZmassWindow,
        c_sufficientMET,
        c_sufficientJets,
        c_sufficientBtags,
        c_sufficientHT,
        c_validSolution,
        m_cuts // Keep as last entry
    };

    const double m_pi = 3.14159265358979323846;
    const double m_mass_b = 4.18, m_mass_W = 80.4, m_mass_Z = 91.19, m_mass_top = 172.5;

    const int m_minBtags;
    int m_bTags;
    
    vector<double>* m_dR_lb_truth;
    vector<double>* m_mass_ttbar_truth;

    // Histograms
    TH1D* h_cutflow;
    TProfile* h_eff_cuts_mass_ttbar_truth;
    TProfile* h_eff_cuts_pT_top_truth;
    TProfile* h_eff_cuts_pT_tbar_truth;
    TProfile* h_eff_cut_2l_mass_ttbar_truth;
    TProfile* h_eff_cut_oc_mass_ttbar_truth;
    TProfile* h_eff_cut_mll_mass_ttbar_truth;
    TProfile* h_eff_cut_mZ_mass_ttbar_truth;
    TProfile* h_eff_cut_ETmiss_mass_ttbar_truth;
    TProfile* h_eff_cut_HT_mass_ttbar_truth;
    TProfile* h_eff_cut_2j_mass_ttbar_truth;
    TProfile* h_eff_cut_2b_mass_ttbar_truth;
    TProfile* h_eff_reco_mass_ttbar_truth;
    TProfile* h_eff_reco_pT_top_truth;
    TProfile* h_eff_reco_pT_tbar_truth;

    TH1D* h_pT_l1;
    TH1D* h_eta_l1;
    TH1D* h_pT_l2;
    TH1D* h_eta_l2;

    TH1D* h_pT_jets;
    TH1D* h_eta_jets;
    TH1D* h_pT_bjets;
    TH1D* h_eta_bjets;
    TH1D* h_pT_qjets;
    TH1D* h_eta_qjets;

    TH1D* h_pT_alljets;
    TH1D* h_pT_allel;
    TH1D* h_pT_allmu;
    TH1D* h_eta_alljets;
    TH1D* h_eta_allel;
    TH1D* h_eta_allmu;

    TH1D* h_HT;
    TH1D* h_HT_truth;
    TH1D* h_HTmet_truth;
    TH1D* h_HTmet;
    TH1D* h_HTjMET;
    TH1D* h_KT_truth;
    TH1D* h_KT;
    TH1D* h_KTj;
    TH1D* h_mass_vis;
    TH1D* h_mass_bbll;
    
    TH1D* h_HT_all;
    TH1D* h_KT_all;
    TH1D* h_mass_vis_all;

    TH1D* h_deltaPhi_ll;
    TH1D* h_cosPhi;

    TH1D* h_mass_W1;
    TH1D* h_mass_W2;

    // top
    TH1D* h_pT_top;
    TH1D* h_eta_top;
    TH1D* h_phi_top;
    TH1D* h_mass_top;
    TH1D* h_E_top;
    TH1D* h_pT_top_truth;
    TH1D* h_eta_top_truth;
    TH1D* h_phi_top_truth;
    TH1D* h_mass_top_truth;

    // tbar
    TH1D* h_pT_tbar;
    TH1D* h_eta_tbar;
    TH1D* h_phi_tbar;
    TH1D* h_mass_tbar;
    TH1D* h_E_tbar;
    TH1D* h_pT_tbar_truth;
    TH1D* h_eta_tbar_truth;
    TH1D* h_phi_tbar_truth;
    TH1D* h_mass_tbar_truth;

    // ttbar
    TH1D* h_pT_ttbar;
    TH1D* h_pT_ttbar_truth;
    TH1D* h_eta_ttbar;
    TH1D* h_eta_ttbar_truth;
    TH1D* h_phi_ttbar;
    TH1D* h_phi_ttbar_truth;
    TH1D* h_mass_ttbar;
    TH1D* h_mass_ttbar_truth;
    TH1D* h_y_ttbar;
    TH1D* h_y_ttbar_truth;

    // dR between truth and reco
    TH1D* h_dR_top;
    TH1D* h_dR_tbar;
    TH1D* h_dR_ttbar;
    
    // dR between each truth lepton and its closest reco lepton
    TH1D* h_dR_l;
    
    // dR between truth leptons and b-quarks and top quarks
    TH1D* h_dR_t1t2_truth;
    TH1D* h_dR_lb_truth;
    // TH1D* h_dR_l1b1_truth;
    // TH1D* h_dR_l2b2_truth;
    TH1D* h_dR_t1l1_truth;
    TH1D* h_dR_t2l2_truth;
    TH1D* h_dR_t1b1_truth;
    TH1D* h_dR_t2b2_truth;
    
    TH2D* h2_dR_l1b1_pTl_truth;
    TH2D* h2_dR_l2b2_pTl_truth;
    TH2D* h2_dR_lb_mtt_truth;
    // TH2D* h2_dR_l1b1_mtt_truth;
    // TH2D* h2_dR_l2b2_mtt_truth;
    TH2D* h2_dR_t1l1_mtt_truth;
    TH2D* h2_dR_t2l2_mtt_truth;
    TH2D* h2_dR_t1b1_mtt_truth;
    TH2D* h2_dR_t2b2_mtt_truth;
    
    TH1D* h_ETmiss;
    TH1D* h_ETmiss_truth;
    TH2D* h2_mtt_truth_ETmiss_truth;

    // performance plots
    TH1D* h_perf_mass_top;
    TH1D* h_perf_pT_top;
    TH1D* h_perf_eta_top;
    TH1D* h_perf_phi_top;
    TH1D* h_perf_mass_tbar;
    TH1D* h_perf_pT_tbar;
    TH1D* h_perf_eta_tbar;
    TH1D* h_perf_phi_tbar;
    TH1D* h_perf_mass_ttbar;
    TH1D* h_perf_pT_ttbar;
    TH1D* h_perf_eta_ttbar;
    TH1D* h_perf_phi_ttbar;
    TH2D* h2_perf_mass_ttbar;
    TH2D* h2_perf_mass_ttbar_pTtop;
    TH2D* h2_perf_mass_ttbar_pTtbar;


    TH2D* h2_mass_top_TvR;
    TH2D* h2_pT_top_TvR;
    TH2D* h2_eta_top_TvR;
    TH2D* h2_phi_top_TvR;
    TH2D* h2_mass_tbar_TvR;
    TH2D* h2_pT_tbar_TvR;
    TH2D* h2_eta_tbar_TvR;
    TH2D* h2_phi_tbar_TvR;
    TH2D* h2_mass_ttbar_TvR;
    TH2D* h2_pT_ttbar_TvR;
    TH2D* h2_eta_ttbar_TvR;
    TH2D* h2_phi_ttbar_TvR;

    // charge asymmetries
    TH1D* h_mtt_tF;
    TH1D* h_mtt_tB;
    TH1D* h_mtt_lF;
    TH1D* h_mtt_lB;
    TH1D* h_HT_lF;
    TH1D* h_HT_lB;
    TH1D* h_KT_lF;
    TH1D* h_KT_lB;
    TH1D* h_mtt_tCF;
    TH1D* h_mtt_tCB;
    TH1D* h_mtt_lCF;
    TH1D* h_mtt_lCB;
    TH1D* h_HT_lCF;
    TH1D* h_HT_lCB;
    TH1D* h_KT_lCF;
    TH1D* h_KT_lCB;

    // TH1D* h_mtt_tlF;
    // TH1D* h_mtt_tlB;

    // top polarisation
    TH1D* h_mtt_c1F;
    TH1D* h_mtt_c1B;
    TH1D* h_mtt_c2F;
    TH1D* h_mtt_c2B;
    // TH1D* h_mtt_philF;
    // TH1D* h_mtt_philB;
    // TH1D* h_mtt_ElF;
    // TH1D* h_mtt_ElB;

    // spin correlation
    TH1D* h_mtt_c1c2F;
    TH1D* h_mtt_c1c2B;
    TH1D* h_mtt_cPhiF;
    TH1D* h_mtt_cPhiB;
    TH1D* h_mtt_DphiF;
    TH1D* h_mtt_DphiB;
    TH1D* h_HT_DphiF;
    TH1D* h_HT_DphiB;
    TH1D* h_KT_DphiF;
    TH1D* h_KT_DphiB;

    TH1D* h_cosTheta;
    TH1D* h_cosTheta1;
    TH1D* h_cosTheta2;
    TH1D* h_cos1cos2;
    TH1D* h_cosThetaStar;
    TH1D* h_cosTheta_ttbar;
    TH1D* h_cosTheta_l;
    TH1D* h_cosThetaStar_l;
    TH1D* h_deltaY_top;
    TH1D* h_deltaEta_l;
    
    TH1D* h_n_electrons;
    TH1D* h_n_muons;
    TH1D* h_n_jets;
    TH1D* h_n_bJets;
    
    TProfile* h_lepton_purity;
    TProfile* h_jet_purity;

    TH1D* h_n_truthElectrons;
    TH1D* h_n_truthMuons;
    TH1D* h_n_truthBquarks;
    
    TH1D* h_n_selElectrons;
    TH1D* h_n_selMuons;
    TH1D* h_n_selJets;
    
    TH1D* h_n_jets_no_truth_tag;
    TH1D* h_n_uniqueElectrons;
    TH1D* h_n_uniqueMuons;
    TH1D* h_n_uniqueJets;
    
    TH1D* h_n_electrons_truth_tagged;
    TH1D* h_n_selElectrons_truth_tagged;
    TH1D* h_n_uniqueElectrons_truth_tagged;
    
    TH1D* h_n_muons_truth_tagged;
    TH1D* h_n_selMuons_truth_tagged;
    TH1D* h_n_uniqueMuons_truth_tagged;
    
    TH1D* h_n_jets_truth_tagged;
    TH1D* h_n_selJets_truth_tagged;
    TH1D* h_n_uniqueJets_truth_tagged;
    TH1D* h_ntracks_truth_tagged_jets;

    TH1D* h_n_passElectrons;
    TH1D* h_n_passMuons;
    TH1D* h_n_passJets;
    TH1D* h_n_passBjets;

    TH2D* h2_mtt_cosThetaStar;
    TH2D* h2_mtt_delta_yt;
    TH2D* h2_mtt_cosThetaStar_ll;
    TH2D* h2_mtt_delta_abs_etal;
    TH2D* h2_HT_delta_abs_etal;
    TH2D* h2_KT_delta_abs_yt;
    TH2D* h2_HT_cosThetaStar_ll;
    TH2D* h2_KT_cosThetaStar_ll;
    TH2D* h2_mtt_deltaPhi;
    TH2D* h2_mtt_cosTheta1;
    TH2D* h2_mtt_cosTheta2;
    TH2D* h2_mtt_cos1cos2;
    TH2D* h2_mvis_deltaPhi;
    TH2D* h2_HT_deltaPhi;
    TH2D* h2_KT_deltaPhi;

protected:
    void SetupTreesForNewFile(const string&);
    void CleanUp();
    void SetupInputFiles();
    void SetupInputFile();
    bool SetupOutputFile();

    void PreLoop();
    void PreLoopSingle();
    void Loop();
    void PostLoop();

    void EachEvent(double);
    void EveryEvent(double);
    void CleanupEvent();
    void GetHardParticles();
    void GetTruthParticles();
    // void GetElectrons();
    // void GetMuons();
    void SelectElectrons();
    void SelectMuons();
    void SelectJets();
    void IsolateElectrons();
    void IsolateMuons();
    void OverlapRemoval();
    void RemoveJetsCloseToElectrons();
    void RemoveJetsCloseToMuons();
    void RemoveElectronsInsideJets();
    void RemoveMuonsInsideJets();
    void TruthTagLeptons();
    void TruthTagJets();
    void FillPurities(int);
    void AssignChannel();
    void FillTaggedHistograms();
    pair<TLorentzVector, TLorentzVector> GetLeptonMomenta();

    // histograms
    void MakeHistograms();
    void MakeDistributions();
    void WriteEfficiency(TH1D*, const string&, const string&);
    void MakeDistribution1D(TH1D*, const string&, bool normalise = false);
    void MakeDistribution2D(TH2D*, string, string, string, string);
    void MakeDistributionAL(TH2D*, const string&, const string&);
    void NormalizeSliceY(TH2D*);
    void AverageEachXbin(TH2D*);

    void GetGenerationCrossSection(int);
    void GetProcessWeight(int);
    void SetDataDirectory();
    void AsymmetryUncertainty(TH1D*, TH1D*, TH1D*);
    void GetBranches();
    void EachFile(const string&);

    // cutflow
    void InitialiseCutflow();
    void PrintCutflow();
    void UpdateCutflow(int, bool);
    void FillCutsEfficiencies(const vector<double>&, const int);
    void FillRecoEfficiencies(const vector<double>&, const int);

    // event selection
    bool PassesEventSelection();
    bool ExactlyTwoLeptons();
    bool OppositeCharge();
    bool SufficientJets();
    bool SufficientBtags();
    bool SufficientHT(double);
    bool SufficientMET(double);
    bool SufficientMll(const pair<TLorentzVector, TLorentzVector>&);
    bool OutsideZmassWindow(const pair<TLorentzVector, TLorentzVector>&);

    Long64_t TotalEvents();
    Long64_t IncrementEvent(Long64_t i);
    double TotalAsymmetry(TH1D* h_A, TH1D* h_B);
    void Asymmetry(const string&, const string&, const string&, TH1D*, TH1D*);

    // tuple
    TClonesArray* b_Particle;
    TClonesArray* b_Jet;
    TClonesArray* b_Electron;
    TClonesArray* b_Muon;
    TClonesArray* b_MissingET;
    TClonesArray* b_EFlowTrack;
    TClonesArray* b_EFlowPhoton;
    TClonesArray* b_EFlowNeutralHadron;
    // TClonesArray* b_ScalarHT;
    TClonesArray* b_Track;

public:
    Analysis(const string& model, const string& process, const string& options, const int energy, const int luminosity, const int minBtags, const string& reconstruction, const string& tag, const bool slice):
        m_inputFileName(""),
        m_processfilename(""),
        m_model(model),
        m_process(process),
        m_options(options),
        m_energy(energy),
        m_luminosity(luminosity),
        m_minBtags(minBtags),
        m_reconstruction(reconstruction),
        m_tag(tag),
        m_useMassSlices(slice),
        m_xSec(false),
        m_debug(false),
        m_outputFile(nullptr),
        m_input(nullptr),
        m_processes(nullptr),
        m_chain(nullptr),
        m_tree(nullptr)
    {
        this->PreLoop();
    }

    Analysis(const string& inputfilename, const string& processfilename, const int luminosity, const int minBtags, const string& reconstruction, const string& tag, const bool slice):
        m_inputFileName(inputfilename),
        m_processfilename(processfilename),
        m_model(""),
        m_process(""),
        m_options(""),
        m_energy(0),
        m_luminosity(luminosity),
        m_minBtags(minBtags),
        m_reconstruction(reconstruction),
        m_tag(tag),
        m_useMassSlices(slice),
        m_xSec(false),
        m_debug(false),
        m_outputFile(nullptr),
        m_input(nullptr),
        m_processes(nullptr),
        m_chain(nullptr),
        m_tree(nullptr)
    {   
        this->PreLoopSingle();
    }
    virtual ~Analysis();
    void Run();
};
#endif
