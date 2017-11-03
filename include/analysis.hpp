#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <fstream>
#include <tuple>
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/asio.hpp"
#include "TH2.h"
#include "TF1.h"
#include "TSystem.h"
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

    typedef vector<tuple<string, int>>::const_iterator itr_s;
    vector<tuple<string, int>>* m_input;
    vector<tuple<string, int, int, double, double, double>>* m_processes;

    vector<Electron*>* m_electrons;
    vector<GenParticle*>* m_truthElectrons;
    vector<Muon*>* m_muons;
    vector<GenParticle*>* m_truthMuons;
    vector<Jet*>* m_jets;
    vector<GenParticle*>* m_truthBquarks;

    GenParticle* m_hardTop;
    GenParticle* m_hardTbar;
    GenParticle* m_hardB;
    GenParticle* m_hardBbar;
    GenParticle* m_hardLepP;
    GenParticle* m_hardLepM;
    GenParticle* m_hardNu;
    GenParticle* m_hardNuBar;

    string m_inputfilename;
    string m_processfilename;
    string m_model;
    string m_process;
    string m_options;
    int m_energy;
    int m_luminosity;
    string m_tag;
    string m_pdf = "CT14LL";
    bool m_use_mass_slices = false;
    bool m_parallel;

    bool m_xSec = false;
    const string m_reconstruction;
    bool m_useLumi;
    const bool m_debug;
    string m_channel;
    bool m_truth;

    double m_crossSection;
    Long64_t m_nevents;

    vector<double> iteration_weights;
    string m_dataDirectory;
    string m_outputName;
    TFile* m_output;
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
        c_sufficientHT,
        c_sufficientJets,
        c_sufficientBtags,
        c_validSolution,
        // c_deltaR,
        m_cuts // Keep as last entry
    };

    const double m_pi = 3.14159265358979323846;
    const double m_mass_b = 4.18, m_mass_W = 80.4, m_mass_Z = 91.19, m_mass_top = 172.5;

    const int m_minBtags;
    int m_bTags;

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
    TH1D* h_KT;
    TH1D* h_mass_vis;
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

    TH1D* h_nTruthElectrons;
    TH1D* h_nTruthMuons;
    TH1D* h_nTruthBquarks;

    TH1D* h_nPassElectrons;
    TH1D* h_nPassMuons;
    TH1D* h_nPassJets;
    TH1D* h_nPassBjets;

    TH1D* h_nElectrons;
    TH1D* h_nMuons;
    TH1D* h_nJets;
    TH1D* h_nBjets;

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
    void SetupOutputFiles();

    void SetupInputFile();
    void SetupOutputFile();

    void PreLoop();
    void PreLoopSingle();
    void Loop();
    void PostLoop();

    void EachEvent(double);
    void EveryEvent(double);
    void CleanupEvent();
    void GetHardParticles();
    void GetTruthParticles();
    void GetElectrons();
    void GetMuons();
    void GetJets();
    void AssignChannel();
    pair<TLorentzVector, TLorentzVector> GetLeptonMomenta();

    // histograms
    void MakeHistograms();
    void MakeDistributions();
    void WriteEfficiency(TH1D*, const string&);
    void MakeDistribution1D(TH1D*, const string&, bool normalise = false);
    void MakeDistribution2D(TH2D*, string, string, string, string);
    void MakeDistributionAL(TH2D*, const string&, const string&);
    void NormalizeSliceY(TH2D*);

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
    bool SufficientHT();
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
    TClonesArray* b_ScalarHT;

public:
    Analysis(const string& model, const string& process, const string& options, const int energy, const int luminosity, const int minBtags, const string& reconstruction, const string& tag, const bool slice):
        m_inputfilename(""),
        m_processfilename(""),
        m_model(model),
        m_process(process),
        m_options(options),
        m_energy(energy),
        m_luminosity(luminosity),
        m_minBtags(minBtags),
        m_reconstruction(reconstruction),
        m_tag(tag),
        m_use_mass_slices(slice),
        m_truth(true),
        m_debug(true),
        m_output(nullptr),
        m_input(nullptr),
        m_processes(nullptr),
        m_chain(nullptr),
        m_tree(nullptr),
        m_parallel(false)
    {
        cout << "parallel: " << m_parallel << "\n";
        this->PreLoop();
    }

    Analysis(const string& inputfilename, const string& processfilename, const int luminosity, const int minBtags, const string& reconstruction, const string& tag, const bool slice):
        m_inputfilename(inputfilename),
        m_processfilename(processfilename),
        m_model(""),
        m_process(""),
        m_options(""),
        m_energy(0),
        m_luminosity(luminosity),
        m_minBtags(minBtags),
        m_reconstruction(reconstruction),
        m_tag(tag),
        m_use_mass_slices(slice),
        m_truth(true),
        m_debug(true),
        m_output(nullptr),
        m_input(nullptr),
        m_processes(nullptr),
        m_chain(nullptr),
        m_tree(nullptr),
        m_parallel(true)
    {
        cout << "parallel: " << m_parallel << "\n";
        this->PreLoopSingle();
    }
    virtual ~Analysis();
    void Run();
};
#endif
