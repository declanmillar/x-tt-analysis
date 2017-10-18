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
        c_realSolutions,
        // c_deltaR,
        m_cuts // Keep as last entry
    };

    const double m_pi = 3.14159265358979323846;
    const double m_bmass = 4.18, m_Wmass = 80.4, m_zmass = 91.19, m_tmass = 172.5;

    const int m_minBtags;
    int m_bTags;

    // Histograms
    TH1D* h_cutflow;

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

    TH1D* h_deltaPhi_ll;
    TH1D* h_cosPhi;

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

    // TH1D* h_pxt_truth;
    // TH1D* h_pyt_truth;
    // TH1D* h_pzt_truth;
    // TH1D* h_Et_truth;
    TH1D* h_pTt_truth;
    TH1D* h_etat_truth;
    TH1D* h_phit_truth;
    TH1D* h_mt_truth;
    // TH1D* h_pxtbar_truth;
    // TH1D* h_pytbar_truth;
    // TH1D* h_pztbar_truth;
    // TH1D* h_Etbar_truth;
    TH1D* h_pTtbar_truth;
    TH1D* h_etatbar_truth;
    TH1D* h_phitbar_truth;
    TH1D* h_mtbar_truth;

    // tt
    TH1D* h_mtt;
    TH1D* h_mttTruth;

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

    TH1D* h_ytt;

    TH1D* h_cosTheta;
    TH1D* h_cosTheta1;
    TH1D* h_cosTheta2;
    TH1D* h_cos1cos2;
    TH1D* h_cosThetaStar;
    TH1D* h_cosTheta_tt;
    TH1D* h_cosTheta_l;
    TH1D* h_cosThetaStar_l;
    TH1D* h_deltaR_max;
    TH1D* h_deltaR_bW;
    TH1D* h_deltaR_tt;
    TH1D* h_delta_abs_yt;
    TH1D* h_delta_abs_etal;

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
    void MakeDistribution1D(TH1D*, const string&);
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
        m_debug(false),
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
        m_debug(false),
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
