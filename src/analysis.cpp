#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"
#include "lester_mt2_bisect.h"
#include "neutrino-weighter.hpp"
#include "kinematic-reconstructer.hpp"
#include "highest-pt.hpp"
#include "match-bjets-to-leps.hpp"
#include "get-parameter.hpp"
#include <exception>
#include <regex>

Analysis::Analysis(const string& model, const string& process, const string& options, const int energy, const int luminosity, const int minimumBtags, const string& reconstruction, const string& tag, bool slice):
    m_model(model),
    m_process(process),
    m_options(options),
    m_energy(energy),
    m_luminosity(luminosity),
    m_minimumBtags(minimumBtags),
    m_reconstruction(reconstruction),
    m_tag(tag),
    m_use_mass_slices(slice),
    m_debug(false),
    m_output(nullptr),
    m_input(nullptr),
    m_processes(nullptr),
    m_chain(nullptr),
    m_tree(nullptr)
{
    this->PreLoop();
}


void Analysis::Run() {
    this->Loop();
    this->PostLoop();
}


void Analysis::EachEvent(double weight) {
    UpdateCutflow(c_events, true);

    if (m_debug) cout << "Starting EachEvent ...\n";

    MissingET* missingET = (MissingET*) b_MissingET->At(0);
    TLorentzVector p_miss;
    double ETmiss = missingET->MET;
    p_miss.SetPtEtaPhiM(ETmiss, missingET->Eta, missingET->Phi, 0.0);

    if (!this->SufficientMET()) return;
    if (!this->SufficientHT()) return;

    m_electron = new vector<Electron*>;
    for (int i = 0; i < b_Electron->GetEntries(); i++) {
        bool passed = false;
        Electron* electron = (Electron*) b_Electron->At(i);
        if (electron->PT > 25.0 and electron->Eta < 2.47) passed = true;
        if (passed) m_electron->push_back(electron);
    }

    m_muon = new vector<Muon*>;
    for (int i = 0; i < b_Muon->GetEntries(); i++) {
        bool passed = false;
        Muon* muon = (Muon*) b_Muon->At(i);
        if (muon->PT > 25.0 and muon->Eta < 2.5) passed = true;
        if (passed) m_muon->push_back(muon);
    }

    if (!this->ExactlyTwoLeptons()) return;
    this->AssignChannel();
    if (!this->OppositeCharge()) return;

    pair< TLorentzVector, TLorentzVector> p_l;
    if (m_channel == "ee") {
        Electron *electron1 = m_electron->at(0);
        Electron *electron2 = m_electron->at(1);

        double charge = electron1->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
    }
    else if (m_channel == "mumu") {
        Muon *muon1 = (Muon*) m_muon->at(0);
        Muon *muon2 = (Muon*) m_muon->at(1);

        double charge = muon1->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
    }
    else if (m_channel == "emu") {
        Electron *electron = (Electron*) m_electron->at(0);
        Muon *muon = (Muon*) m_muon->at(0);

        double charge = electron->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        }
    }
    else cout << "Error: invalid channel\n";
    if (m_debug) p_l.first.Print();
    if (m_debug) p_l.second.Print();

    if (!this->SufficientMll(p_l)) return;
    if (!this->OutsideZmassWindow(p_l)) return;

    m_jet = new vector<Jet*>;
    for (int i = 0; i < b_Jet->GetEntries(); i++) {
        bool passed = false;
        Jet *jet = (Jet*) b_Jet->At(i);
        if (jet->PT > 25.0 and jet->Eta < 2.5) passed = true;
        if (passed) m_jet->push_back(jet);
    }

    if (!this->SufficientJets()) return;
    if (!this->SufficientBtags()) return;

    vector<TLorentzVector> p_j, p_b, p_q;
    for (int i = 0; i < m_jet->size(); i++) {
        Jet *jet = (Jet*) m_jet->at(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        // if (jet->BTag & (1 << i))
        if (jet->BTag > 0) {
            if (m_debug) cout << "b-jet pT = " << p.Pt() << "\n";
            p_b.push_back(p);
            h_pt_bjets->Fill(p.Pt(), weight);
            h_eta_bjets->Fill(p.Eta(), weight);
        }
        else {
            p_q.push_back(p);
            h_pt_qjets->Fill(p.Pt(), weight);
            h_eta_qjets->Fill(p.Eta(), weight);
        }
        p_j.push_back(p);
        h_pt_jets->Fill(p.Pt(), weight);
        h_eta_jets->Fill(p.Eta(), weight);
    }

    // ScalarHT *scalarHT = (ScalarHT*) b_ScalarHT->At(0);
    // double HT = scalarHT->HT;
    double HT = p_l.first.Pt() + p_l.second.Pt() + ETmiss;
    for (auto& p : p_b) HT += p.Pt();
    for (auto& p : p_q) HT += p.Pt();

    TLorentzVector pvis = p_l.first + p_l.second;
    for (auto& p : p_j) pvis += p;
    double mvis = pvis.M();
    double pTvis = pvis.Pt();
    double KT = sqrt(mvis * mvis + pTvis * pTvis) + ETmiss;

    h_HT->Fill(HT / 1000, weight);
    h_mvis->Fill(mvis / 1000, weight);
    h_KT->Fill(KT / 1000, weight);

    if (m_debug) {
        cout << "p_l+   = "; p_l.first.Print();
        cout << "p_l-   = "; p_l.second.Print();
        for (int i = 0; i < p_b.size(); i++) {
            cout << "p_b" << i << "   = ";
            p_b.at(i).Print();
        }
        for (int i = 0; i < p_q.size(); i++) {
            cout << "p_q" << i << "   = ";
            p_q.at(i).Print();
        }
        cout << "p_miss = "; p_miss.Print();
    }

    pair< TLorentzVector, TLorentzVector> p_b_hypo;
    if (m_btags >= 2) p_b_hypo = TwoHighestPt(p_b);
    else if (m_btags == 0) p_b_hypo = TwoHighestPt(p_q);
    else if (m_btags == 1) {
        p_b_hypo.first = p_b.at(0);
        p_b_hypo.second = HighestPt(p_q);
    }

    TLorentzVector p_top, p_tbar, p_ttbar, p_b1, p_b2, p_v1, p_v2;

    if (m_debug) cout << "Reconstruction method: " << m_reconstruction << "\n";

    if (m_reconstruction == "KIN") {
        if (m_debug) cout << "Setting up kinematic reconstruction ...\n";
        vector<TLorentzVector> p_b_hiPt = { p_b_hypo.first, p_b_hypo.second };
        KinematicReconstructer KIN = KinematicReconstructer(m_bmass, m_Wmass, m_tmass);
        bool isSolution = KIN.Reconstruct(p_l, p_b_hiPt, p_miss);
        p_top = KIN.GetTop();
        p_tbar = KIN.GetTbar();
        p_ttbar = KIN.GetTtbar();
        p_b1 = KIN.GetB();
        p_b2 = KIN.GetBbar();
        p_v1 = KIN.GetNu();
        p_v2 = KIN.GetNubar();
        if (isSolution) this->UpdateCutflow(c_realSolutions, true);
        else {
            this->UpdateCutflow(c_realSolutions, false);
            return;
        }
        if (m_debug) cout << "Finished kinematic reconstruction\n";
    }
    else if (m_reconstruction == "NuW") {
        auto p_b_match = MatchBjetsToLeps(p_l, p_b_hypo);

        NeutrinoWeighter nuW = NeutrinoWeighter(1, p_l.first.Pt() + p_l.first.Phi()); // 2nd argument is random seed same for specific event
        double weight_max = nuW.Reconstruct(p_l.first, p_l.second, p_b_match.first, p_b_match.second, p_miss.Px(), p_miss.Py(), p_miss.Phi());
        if (weight_max > 0.0) {
            p_top = nuW.GetTop();
            p_tbar = nuW.GetTbar();
            p_ttbar = nuW.GetTtbar();
            p_b1 = nuW.GetB();
            p_b2 = nuW.GetBbar();
            p_v1 = nuW.GetNu();
            p_v2 = nuW.GetNubar();
        }
    }
    else cout << "Error: Reconstruction = {KIN, NuW}\n";

    TLorentzVector p_W1 = p_l.first + p_v1;
    TLorentzVector p_W2 = p_l.second + p_v2;

    auto mtt = p_ttbar.M();
    auto ytt = p_ttbar.Rapidity();
    auto costheta_tt = cos(p_top.Angle(p_tbar.Vect()));
    auto delta_abs_yt = abs(p_top.Rapidity()) - abs(p_tbar.Rapidity());

    auto pttbar_top = p_top, pttbar_tbar = p_tbar, pttbar_ttbar = p_ttbar;
    auto pttbar_b1 = p_b1, pttbar_b2 = p_b2, pttbar_v1 = p_v1, pttbar_v2 = p_v2;
    auto pttbar_l = p_l;

    TVector3 v_ttbar = - (p_ttbar).BoostVector();
    pttbar_top.Boost(v_ttbar);
    pttbar_tbar.Boost(v_ttbar);
    pttbar_ttbar.Boost(v_ttbar);
    pttbar_b1.Boost(v_ttbar);
    pttbar_b2.Boost(v_ttbar);
    pttbar_l.first.Boost(v_ttbar);
    pttbar_l.second.Boost(v_ttbar);
    pttbar_v1.Boost(v_ttbar);
    pttbar_v2.Boost(v_ttbar);

    auto ptop_top = p_top, ptop_tbar = p_tbar, ptop_ttbar = p_ttbar;
    auto ptop_b1 = p_b1, ptop_b2 = p_b2, ptop_v1 = p_v1, ptop_v2 = p_v2;
    auto ptop_l = p_l;

    TVector3 v_top = - (p_top).BoostVector();
    ptop_top.Boost(v_top);
    ptop_tbar.Boost(v_top);
    ptop_ttbar.Boost(v_top);
    ptop_b1.Boost(v_top);
    ptop_b2.Boost(v_top);
    ptop_l.first.Boost(v_top);
    ptop_l.second.Boost(v_top);
    ptop_v1.Boost(v_top);
    ptop_v2.Boost(v_top);

    auto ptbar_top = p_top, ptbar_tbar = p_tbar, ptbar_ttbar = p_ttbar;
    auto ptbar_b1 = p_b1, ptbar_b2 = p_b2, ptbar_v1 = p_v1, ptbar_v2 = p_v2;
    auto ptbar_l = p_l;

    TVector3 v_tbar = - (p_tbar).BoostVector();
    ptbar_top.Boost(v_tbar);
    ptbar_tbar.Boost(v_tbar);
    ptbar_ttbar.Boost(v_tbar);
    ptbar_b1.Boost(v_tbar);
    ptbar_b2.Boost(v_tbar);
    ptbar_l.first.Boost(v_tbar);
    ptbar_l.second.Boost(v_tbar);
    ptbar_v1.Boost(v_tbar);
    ptbar_v2.Boost(v_tbar);

    double costheta = pttbar_top.CosTheta();
    double costhetastar = int(ytt / abs(ytt)) * costheta;

    double deltaPhi = p_l.first.DeltaPhi(p_l.second) / m_pi;


    double costheta_tl1 = cos(ptop_l.first.Angle(pttbar_top.Vect()));
    double costheta_tl2 = cos(ptbar_l.second.Angle(pttbar_tbar.Vect()));

    if (m_debug) cout << "Filling histograms ...\n";

    mtt = mtt / 1000;

    h_pt_l1->Fill(p_l.first.Pt(), weight);
    h_pt_l2->Fill(p_l.second.Pt(), weight);
    h_eta_l1->Fill(p_l.first.Eta(), weight);
    h_eta_l2->Fill(p_l.second.Eta(), weight);

    h_Et->Fill(p_top.E(), weight);
    h_pTt->Fill(p_top.Pt(), weight);
    h_etat->Fill(p_top.Eta(), weight);
    h_phit->Fill(p_top.Phi() / m_pi, weight);
    h_mt->Fill(p_top.M(), weight);

    h_Etbar->Fill(p_tbar.E(), weight);
    h_pTtbar->Fill(p_tbar.Pt(), weight);
    h_etatbar->Fill(p_tbar.Eta(), weight);
    h_phitbar->Fill(p_tbar.Phi() / m_pi, weight);
    h_mtbar->Fill(p_tbar.M(), weight);

    h_mtt->Fill(mtt, weight);
    h_ytt->Fill(ytt, weight);

    h_mW1->Fill(p_W1.M(), weight);
    h_mW2->Fill(p_W2.M(), weight);

    h_deltaPhi->Fill(deltaPhi, weight);

    h2_HT_deltaPhi->Fill(HT / 1000, deltaPhi, weight);
    h2_mvis_deltaPhi->Fill(mvis / 1000, deltaPhi, weight);
    h2_KT_deltaPhi->Fill(KT / 1000, deltaPhi, weight);

    h_costheta_tt->Fill(costheta_tt, weight);
    h_cosTheta->Fill(costheta, weight);
    h_cosThetaStar->Fill(costhetastar, weight);

    if (costhetastar > 0) h_mtt_tF->Fill(mtt, weight);
    if (costhetastar < 0) h_mtt_tB->Fill(mtt, weight);

    if (delta_abs_yt > 0) h_mtt_tCF->Fill(mtt, weight);
    if (delta_abs_yt < 0) h_mtt_tCB->Fill(mtt, weight);

    double cos1cos2 = costheta_tl1 * costheta_tl2;

    h_cosTheta1->Fill(costheta_tl1, weight);
    h_cosTheta2->Fill(costheta_tl2, weight);
    h_cos1cos2->Fill(cos1cos2, weight);

    if (costheta_tl1 > 0) h_mtt_tlF->Fill(mtt, weight);
    if (costheta_tl2 < 0) h_mtt_tlB->Fill(mtt, weight);

    if (costheta_tl1 > 0) h_mtt_lF->Fill(mtt, weight);
    if (costheta_tl1 < 0) h_mtt_lB->Fill(mtt, weight);

    h2_mtt_cosThetaStar->Fill(mtt, costhetastar, weight);
    h2_mtt_deltaPhi->Fill(mtt, deltaPhi, weight);
    h2_mtt_cosTheta1->Fill(mtt, costheta_tl1, weight);
    h2_mtt_cosTheta2->Fill(mtt, costheta_tl2, weight);
    h2_mtt_cos1cos2->Fill(mtt, cos1cos2, weight);

    if (m_debug) cout << "Finished filling histograms\n";
    this->CleanupEvent();
}

void Analysis::CleanupEvent() {
    if (m_debug) cout << "Cleaning up event\n";
    delete m_electron;
    delete m_muon;
    delete m_jet;
    if (m_debug) cout << "Finished cleaning up event\n";
}

void Analysis::EveryEvent(double weight) {
    if (m_debug) cout << "Starting EveryEvent ...\n";

    // runs for every event with no event selection or cuts
    h_nElectrons->Fill(b_Electron->GetEntries(), weight);
    h_nMuons->Fill(b_Muon->GetEntries(), weight);
    h_nJets->Fill(b_Jet->GetEntries(), weight);

    if (m_debug) cout << "Fetching all jets ...\n";
    vector<TLorentzVector> p_j;
    for (int i = 0; i < b_Jet->GetEntries(); i++) {
        Jet *jet = (Jet*) b_Jet->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        h_pt_alljets->Fill(p.Pt(), weight);
        h_eta_alljets->Fill(p.Eta(), weight);
        p_j.push_back(p);
    }

    if (m_debug) cout << "Fetching all electrons ...\n";
    vector<TLorentzVector> p_el;
    for (int i = 0; i < b_Electron->GetEntries(); i++) {
        Electron *electron = (Electron*) b_Electron->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
        h_pt_allel->Fill(p.Pt(), weight);
        h_eta_allel->Fill(p.Eta(), weight);
        p_el.push_back(p);
    }

    if (m_debug) cout << "Fetching all muons ...\n";
    vector<TLorentzVector> p_mu;
    for (int i = 0; i < b_Muon->GetEntries(); i++) {
        Muon *muon = (Muon*) b_Muon->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        h_pt_allmu->Fill(p.Pt(), weight);
        h_eta_allmu->Fill(p.Eta(), weight);
        p_mu.push_back(p);
    }

    if (m_debug) cout << "Fetching missing ET ...\n";
    MissingET *missingET = (MissingET*) b_MissingET->At(0);
    TLorentzVector p_miss;
    double ETmiss = missingET->MET;

    if (m_debug) cout << "Calculating HT ...\n";
    double HT = 0;
    for (auto p : p_j) HT += p.Pt();
    for (auto p : p_el) HT += p.Pt();
    for (auto p : p_mu) HT += p.Pt();
    HT += ETmiss;
    h_HT_all->Fill(HT / 1000, weight);

    if (m_debug) cout << "Calculating pvis ...\n";
    TLorentzVector pvis(0, 0, 0, 0);
    for (auto p : p_j) pvis += p;
    for (auto p : p_el) pvis += p;
    for (auto p : p_mu) pvis += p;

    if (m_debug) cout << "Calculating mvis ...\n";
    double mvis = pvis.M();
    h_mvis_all->Fill(mvis / 1000, weight);

    if (m_debug) cout << "Calculating KT ...\n";
    double pTvis = pvis.Pt();
    double KT = sqrt(mvis * mvis + pTvis * pTvis) + ETmiss;
    h_KT_all->Fill(KT / 1000, weight);

    if (m_debug) cout << "Finished EveryEvent ...\n";
}


void Analysis::SetupInputFiles() {
    m_input = new vector< tuple<string, int> >;
    m_processes = new vector< tuple<string, int, int, double, double, double> >;
    string filename;

    string E = to_string(m_energy);

    vector<string> initials = { "gg", "qq", "dd", "uu" };

    int proc_id = 0;
    for (auto initial : initials) {

        // check initial state has been specified for analysis
        size_t pos = m_process.find("-tt");
        string final_state = m_process.substr(pos);
        if (!boost::contains(m_process, initial)) continue;

        // get model name (if only QCD processes - set to Standard Model)
        string model = "";
        if (initial == "gg" or initial == "qq") model = "SM";
        else model = m_model;

        // infer intermediate EW sector bosons based on model and intial states
        string intermediates = "";
        if (initial == "uu" or initial == "dd") {
            intermediates = intermediates + "-AZ";
            if (m_model != "SM") intermediates = intermediates + "X";
            // intermediates = "-X";
        }

        // get options
        string options = "";
        options = m_options;

        // combine for file name
        filename = initial + intermediates + final_state + "_" + model + "_" + E + "TeV" + "_" + m_pdf + options;



        cout << "Adding:         " << filename << "*_pythia_delphes.root ...\n";

        // loop over all matching files (e.g. *_01.root and *_02.root)
        boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
        int nfiles = 0;

        if (m_use_mass_slices) cout << "I AM USING MASS SLICES\n";

        int end = 1;
        if (m_use_mass_slices) end =  m_energy;

        for (int j = 0; j < end; j++) {
            string range = "";
            if (m_use_mass_slices) range = "_" + to_string(j) + "-" + to_string(j + 1);

            int nfiles_per_slice = 0;
            for (boost::filesystem::directory_iterator i(m_dataDirectory); i != end_itr; ++i) {

                if (!boost::filesystem::is_regular_file(i->status())) continue;
                if (m_debug) cout << "is file: " << i->path().filename().string() << "\n";

                if (!boost::contains(i->path().filename().string(), filename)) continue;
                if (m_debug) cout << "contains " << filename << ": " << i->path().filename().string() << "\n";

                if (boost::contains(i->path().filename().string(), "KIN")) continue;
                if (boost::contains(i->path().filename().string(), "NuW")) continue;
                if (m_debug) cout << "not output: " << i->path().filename().string() << "\n";

                if (!boost::contains(i->path().filename().string(), "_pythia_delphes")) continue;
                if (m_debug) cout << "has _pythia_delphes suffix: " << i->path().filename().string() << "\n";

                if (!boost::contains(i->path().filename().string(), range)) continue;

                regex reg(filename + "_[0-9]-[0-9]_[0-9]+_pythia_delphes");
                if (!m_use_mass_slices and regex_search( i->path().filename().string(), reg)) continue;

                if (i->path().extension() == ".root") {
                    // cout << "ends .root: " << i->path().filename().string() << "\n";
                    nfiles++;
                    nfiles_per_slice++;
                    tuple< string, int > input = make_tuple(m_dataDirectory + i->path().filename().string(), proc_id);
                    m_input->push_back(input);
                    if (nfiles < 10) cout << "Input " << nfiles << ":        " << get<0>(input);
                    else if (nfiles < 100) cout << "Input " << nfiles << ":       " << get<0>(input);
                    else cout << "Input " << nfiles << ":      " << get<0>(input);
                    cout << ", process: " << get<1>(input) << "\n";
                }
            }
            if (m_use_mass_slices and nfiles_per_slice == 0) {
                cout << "No files in energy range " << range << " [TeV]\n";
                continue;
            }
            string proc_filename = m_dataDirectory + initial + intermediates + "-tt-bbllvv" + "_" + model + "_" + E + "TeV" + "_" + m_pdf + options + range + ".txt";
            cout << "Adding process: " << proc_filename << " ...\n";
            tuple< string, int, int, double, double, double > process = make_tuple(proc_filename, proc_id, nfiles, -999, -999, -999);
            m_processes->push_back(process);
            proc_id++;
        }
    }

    // check some input files have been specified
    if (m_input->size() < 1) {
        cout << "Error: no input files specified.";
        exit(false);
    }

    // check all input files exist
    for (auto input : *m_input) {
        struct stat buffer;
        bool exists = stat((get<0>(input)).c_str(), &buffer) == 0;
        if (exists == false) {
            cout << "Error: no " << get<0>(input) << "\n";
            exit(exists);
        }
    }

    for (auto process : *m_processes) {
        struct stat buffer;
        bool exists = stat((get<0>(process)).c_str(), &buffer) == 0;
        if (exists == false) {
            cout << "Error: no " << get<0>(process) << "\n";
            exit(exists);
        }
    }
}


void Analysis::SetupOutputFiles() {
    string E = to_string(m_energy) + "TeV";
    string L = to_string(m_luminosity) + "fb-1";

    m_outputName = m_dataDirectory + m_process + "_" + m_model + "_" + E + "_" + m_pdf + m_options;
    m_outputName += "_pythia_delphes";
    if (m_use_mass_slices) m_outputName += "_sliced";
    m_outputName += "_" + m_reconstruction + m_tag;
    if (m_luminosity > 0) m_outputName += L;
    m_outputName += ".root";
    m_output = new TFile(m_outputName.c_str(), "RECREATE");
}

void Analysis::PostLoop() {
    this->MakeDistributions();
    this->PrintCutflow();
    m_output->Close();
    cout << "\nOutput\n";
    cout << m_outputName << "\n";
}


TH1D* Analysis::Asymmetry(const string& name, const string& title, TH1D* h1, TH1D* h2) {
    TH1D* h_numerator = (TH1D*) h1->Clone(name.data());
    TH1D* h_denominator = (TH1D*) h1->Clone();
    h_numerator->SetTitle(title.data());
    h_numerator->Add(h2, -1);
    h_denominator->Add(h2, 1);
    h_numerator->Divide(h_denominator);
    delete h_denominator;
    if (m_luminosity > 0) this->AsymmetryUncertainty(h_numerator, h1, h2);
    return h_numerator;
}


void Analysis::AsymmetryUncertainty(TH1D* hA, TH1D* h1, TH1D* h2) {
    double A, dA, N, N1, N2;
    for (int i = 1; i < hA->GetNbinsX() + 1; i++) {
        A = hA->GetBinContent(i);
        N1 = h1->GetBinContent(i);
        N2 = h2->GetBinContent(i);
        N = N1 + N2;
        if (N > 0) dA = sqrt((1.0 - A * A) / N);
        else dA = 0;
        hA->SetBinError(i, dA);
    }
}


void Analysis::MakeHistograms() {
    double binWidth = 0.1;
    double Emin = 0.05;
    double Emax = 12.95;
    double nbins = (Emax - Emin) / binWidth;
     cout << "Plotting range: " << Emin << " - " << Emax << " [TeV]\n";

    h_pt_l1 = new TH1D("pT_l1", "p^{l^{+}}_{T}", 40, 0.0, 1000.0);
    h_pt_l1->Sumw2();
    h_pt_l2 = new TH1D("pT_l2", "p^{l^{-}}_{T}", 40, 0.0, 1000.0);
    h_pt_l2->Sumw2();
    h_eta_l1 = new TH1D("eta_l1", "#eta_{l^{+}}", 60, -3.0, 3.0);
    h_eta_l1->Sumw2();
    h_eta_l2 = new TH1D("eta_l2", "#eta_{l^{-}}", 60, -3.0, 3.0);
    h_eta_l2->Sumw2();

    h_pt_jets = new TH1D("pT_jets", "p_{T}^{jets}", 40, 0.0, 5000.0);
    h_pt_jets->Sumw2();
    h_eta_jets = new TH1D("eta_jets", "#eta_{jets}", 60, -3.0, 3.0);
    h_eta_jets->Sumw2();
    h_pt_bjets = new TH1D("pT_bjets", "p_{T}^{b-jets}", 40, 0.0, 5000.0);
    h_pt_bjets->Sumw2();
    h_eta_bjets = new TH1D("eta_bjets", "#eta_{b-jets}", 60, -3.0, 3.0);
    h_eta_bjets->Sumw2();
    h_pt_qjets = new TH1D("pT_qjets", "p_{T}^{q-jets}", 40, 0.0, 5000.0);
    h_pt_qjets->Sumw2();
    h_eta_qjets = new TH1D("eta_qjets", "#eta_{q-jets}", 60, -3.0, 3.0);
    h_eta_qjets->Sumw2();

    h_mtt = new TH1D("m_tt", "m_{tt}", nbins, Emin, Emax);
    h_mtt->Sumw2();
    h_ytt = new TH1D("y_tt", "y_{tt}", 50, -2.5, 2.5);
    h_ytt->Sumw2();

    h_mW1 = new TH1D("mW1", "m_{W^{+}}", 150, 0.0, 150.0);
    h_mW1->Sumw2();
    h_mW2 = new TH1D("mW2", "m_{W^{-}}", 150, 0.0, 150.0);
    h_mW2->Sumw2();

    h_Et = new TH1D("E_t", "E_{t}", 100, 0.0, 5000.0);
    h_Et->Sumw2();
    h_pTt = new TH1D("pT_t", "p_{T}^{t}", 100, 0.0, 5000.0);
    h_pTt->Sumw2();
    h_etat = new TH1D("eta_t", "#eta_{t}", nbins, 0.0, 10.0);
    h_etat->Sumw2();
    h_phit = new TH1D("phi_t", "#phi_{t}", nbins, -1.0, 1.0);
    h_phit->Sumw2();
    h_mt = new TH1D("m_t", "m_{t}", 40, 100.0, 300.0);
    h_mt->Sumw2();

    h_Etbar = new TH1D("E_tbar", "E_{#bar{t}}", 100, 0.0, 5000.0);
    h_Etbar->Sumw2();
    h_pTtbar = new TH1D("pT_tbar", "p_{T}^{#bar{t}}", 100, 0.0, 5000.0);
    h_pTtbar->Sumw2();
    h_etatbar = new TH1D("eta_tbar", "#eta_{#bar{t}}", nbins, 0.0, 10.0);
    h_etatbar->Sumw2();
    h_phitbar = new TH1D("phi_tbar", "#phi_{#bar{t}}", nbins, -1.0, 1.0);
    h_phitbar->Sumw2();
    h_mtbar = new TH1D("m_tbar", "m_{#bar{t}}", 40, 100.0, 300.0);
    h_mtbar->Sumw2();

    // AtFB
    h_mtt_tF = new TH1D("mtt_tF", "m_{tt}^{tF}", nbins, Emin, Emax);
    h_mtt_tF->Sumw2();
    h_mtt_tB = new TH1D("mtt_tB", "m_{tt}^{tB}", nbins, Emin, Emax);
    h_mtt_tB->Sumw2();

    // AtC
    h_mtt_tCF = new TH1D("mtt_tCF", "m_{tt}^{tCF}", nbins, Emin, Emax);
    h_mtt_tCF->Sumw2();
    h_mtt_tCB = new TH1D("mtt_tCB", "m_{tt}^{tCB}", nbins, Emin, Emax);
    h_mtt_tCB->Sumw2();

    h2_mtt_cosThetaStar = new TH2D("mtt_costhetastar", "m_{tt} cos#theta^{*}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosThetaStar->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_cosThetaStar->GetYaxis()->SetTitle("cos#theta^*");
    h2_mtt_cosThetaStar->Sumw2();

    h_costheta_tt = new TH1D("costheta_tt", "cos#theta_{t,#bar{t}}", 20, -1.0, 1.0);
    h_costheta_tt->Sumw2();

    // AtlFB
    h_mtt_tlF = new TH1D("mtt_tlF", "m_{tt}^{tlF}", nbins, Emin, Emax);
    h_mtt_tlF->Sumw2();
    h_mtt_tlB = new TH1D("mtt_tlB", "m_{tt}^{tlB}", nbins, Emin, Emax);
    h_mtt_tlB->Sumw2();

    // Aphil
    h_mtt_philF = new TH1D("mtt_philF", "m_{tt}^{philF}", nbins, Emin, Emax);
    h_mtt_philF->Sumw2();
    h_mtt_philB = new TH1D("mtt_philB", "m_{tt}^{philB}", nbins, Emin, Emax);
    h_mtt_philB->Sumw2();

    // AlEl
    h_mtt_ElF = new TH1D("mtt_ElF", "m_{tt}^{ElF}", nbins, Emin, Emax);
    h_mtt_ElF->Sumw2();
    h_mtt_ElB = new TH1D("mtt_ElB", "m_{tt}^{ElB}", nbins, Emin, Emax);
    h_mtt_ElB->Sumw2();

    h_cosTheta = new TH1D("costheta", "cos#theta", nbins, -1.0, 1.0);
    h_cosTheta->Sumw2();
    h_cosThetaStar = new TH1D("costheta_star", "cos#theta^{*}", nbins, -1.0, 1.0);
    h_cosThetaStar->Sumw2();

    h_HT = new TH1D("HT", "H_{T}", 60, 0.0, 6.0);
    h_HT->Sumw2();
    h_HT_all = new TH1D("HT_all", "H^{all}_{T}", 50, 0.0, 6.0);
    h_HT_all->Sumw2();

    h_KT = new TH1D("KT", "K_{T}", 60, 0.0, 6.0);
    h_KT->Sumw2();
    h_KT_all = new TH1D("KT_all", "K^{all}_{T}", 50, 0.0, 6.0);
    h_KT_all->Sumw2();

    h_mvis = new TH1D("mvis", "m_{vis}", 60, 0.0, 6.0);
    h_mvis->Sumw2();
    h_mvis_all = new TH1D("mvis_all", "m^{all}_{vis}", 60, 0.0, 6.0);
    h_mvis_all->Sumw2();

    h_pt_alljets = new TH1D("pT_alljets", "p_{T}^{jet}", 1000, 0.0, 5000.0);
    h_pt_alljets->Sumw2();
    h_pt_allel = new TH1D("pT_allel", "p_{T}^{electron}", 1000, 0.0, 5000.0);
    h_pt_allel->Sumw2();
    h_pt_allmu = new TH1D("pT_allmu", "p_{T}^{muon}", 1000, 0.0, 5000.0);
    h_pt_allmu->Sumw2();
    h_eta_alljets = new TH1D("eta_alljets", "#eta^{jet}", 500, 0.0, 5.0);
    h_eta_alljets->Sumw2();
    h_eta_allel = new TH1D("eta_allel", "#eta^{electron}", 250, 0.0, 2.5);
    h_eta_allel->Sumw2();
    h_eta_allmu = new TH1D("eta_allmu", "#eta^{muon}", 250, 0.0, 2.5);
    h_eta_allmu->Sumw2();

    h_deltaPhi = new TH1D("delta_phi", "#Delta#phi", 10, 0.0, 1.0);
    h_deltaPhi->Sumw2();

    h_cosTheta1 = new TH1D("costheta_tl1", "cos#theta_{t,l+}", 10, -1.0, 1.0);
    h_cosTheta1->Sumw2();
    h_cosTheta2 = new TH1D("costheta_tl2", "cos#theta_{t,l-}", 10, -1.0, 1.0);
    h_cosTheta2->Sumw2();
    h_cos1cos2 = new TH1D("cos1cos2", "cos#theta_{t,l+}cos#theta_{t,l-}", 20, -1.0, 1.0);
    h_cos1cos2->Sumw2();

    h_mtt_lF = new TH1D("mtt_lF", "m_{tt}^{F,l}", nbins, Emin, Emax);
    h_mtt_lF->Sumw2();
    h_mtt_lB = new TH1D("mtt_lB", "m_{tt}^{B,l}", nbins, Emin, Emax);
    h_mtt_lB->Sumw2();

    h_nElectrons = new TH1D("n_electrons", "n_{electrons}", 10, 0.0, 10.0);
    h_nElectrons->Sumw2();
    h_nMuons = new TH1D("n_muons", "n_{muons}", 10, 0.0, 10.0);
    h_nMuons->Sumw2();
    h_nJets = new TH1D("n_jets", "n_{jets}", 10, 0.0, 10.0);
    h_nJets->Sumw2();

    h2_mtt_deltaPhi = new TH2D("mtt_deltaphi", "m_{tt} #Delta#phi_{l}", nbins, Emin, Emax, 10, 0.0, 1.0);
    h2_mtt_deltaPhi->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l}");
    h2_mtt_deltaPhi->Sumw2();

    h2_mtt_cosTheta1 = new TH2D("mtt_costheta_tl1", "m_{tt} cos#theta_{t,l+}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta1->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h2_mtt_cosTheta1->GetYaxis()->SetTitle("cos#theta_{l+}");
    h2_mtt_cosTheta1->Sumw2();

    h2_mtt_cosTheta2 = new TH2D("mtt_costheta_tl2", "m_{tt} cos#theta_{t,l-}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta2->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h2_mtt_cosTheta2->GetYaxis()->SetTitle("cos#theta_{l-}");
    h2_mtt_cosTheta2->Sumw2();

    h2_mtt_cos1cos2 = new TH2D("mtt_cos1cos2", "m_{tt} cos#theta_{l+}cos#theta_{l-}", nbins, Emin, Emax, 20, -1.0, 1.0);
    h2_mtt_cos1cos2->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_cos1cos2->GetYaxis()->SetTitle("cos#theta_{l+}cos#theta_{l-}");
    h2_mtt_cos1cos2->Sumw2();

    h2_HT_deltaPhi = new TH2D("HT_deltaphi", "H_{T} #Delta#phi_{l}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_HT_deltaPhi->GetXaxis()->SetTitle("H_{T} [TeV]");
    h2_HT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_HT_deltaPhi->Sumw2();

    h2_mvis_deltaPhi = new TH2D("mvis_deltaphi", "m_{vis} #Delta#phi_{l}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_mvis_deltaPhi->GetXaxis()->SetTitle("m_{vis} [TeV]");
    h2_mvis_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_mvis_deltaPhi->Sumw2();

    h2_KT_deltaPhi = new TH2D("KT_deltaphi", "K_{T} #Delta#phi_{l}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_KT_deltaPhi->GetXaxis()->SetTitle("K_{T} [TeV]");
    h2_KT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_KT_deltaPhi->Sumw2();

    h_deltaR_tt = new TH1D("deltaR_tt", "#Delta R(t,#bar{t})", 100, 0.0, 5.0);
    h_deltaR_tt->Sumw2();

    h_deltaR_bW = new TH1D("deltaR_bW", "#Delta R(b,W)", 100, 0.0, 5.0);
    h_deltaR_bW->Sumw2();
    h_deltaR_max = new TH1D("deltaR_max", "#Delta R_{max}", 100, 0.0, 5.0);
    h_deltaR_max->Sumw2();
}


void Analysis::MakeDistributions() {
    if (m_debug) cout << "Making distributions ...\n";
    this->MakeDistribution1D(h_mtt, "TeV");
    this->MakeDistribution1D(h_ytt, "");

    this->MakeDistribution1D(h_mW1, "TeV");
    this->MakeDistribution1D(h_mW2, "TeV");

    this->MakeDistribution1D(h_Et, "GeV");
    this->MakeDistribution1D(h_pTt, "GeV");
    this->MakeDistribution1D(h_etat, "");
    this->MakeDistribution1D(h_phit, "");
    this->MakeDistribution1D(h_mt, "GeV");

    this->MakeDistribution1D(h_Etbar, "GeV");
    this->MakeDistribution1D(h_pTtbar, "GeV");
    this->MakeDistribution1D(h_etatbar, "");
    this->MakeDistribution1D(h_phitbar, "");
    this->MakeDistribution1D(h_mtbar, "GeV");

    this->MakeDistribution1D(h_pt_alljets, "GeV");
    this->MakeDistribution1D(h_pt_allel, "GeV");
    this->MakeDistribution1D(h_pt_allmu, "GeV");
    this->MakeDistribution1D(h_eta_alljets, "");
    this->MakeDistribution1D(h_eta_allel, "");
    this->MakeDistribution1D(h_eta_allmu, "");

    this->MakeDistribution1D(h_cosTheta, "");
    this->MakeDistribution1D(h_cosThetaStar, "");
    this->MakeDistribution1D(h_costheta_tt, "");

    this->MakeDistribution1D(h_mtt_tF, "TeV");
    this->MakeDistribution1D(h_mtt_tB, "TeV");

    h_AtFB = this->Asymmetry("AtFB", "A^{t}_{FB^{*}}", h_mtt_tF, h_mtt_tB);
    h_AtFB->GetYaxis()->SetTitle(h_AtFB->GetTitle());
    h_AtFB->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AtFB->Write();

    this->MakeDistribution1D(h_mtt_tCF, "TeV");
    this->MakeDistribution1D(h_mtt_tCB, "TeV");

    h_AtC = this->Asymmetry("AtC", "A^{t}_{C}", h_mtt_tCF, h_mtt_tCB);
    h_AtC->GetYaxis()->SetTitle(h_AtC->GetTitle());
    h_AtC->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AtC->Write();

    this->MakeDistribution1D(h_pt_l1, "TeV");
    this->MakeDistribution1D(h_eta_l1, "");
    this->MakeDistribution1D(h_pt_l2, "TeV");
    this->MakeDistribution1D(h_eta_l2, "");
    this->MakeDistribution1D(h_pt_jets, "TeV");
    this->MakeDistribution1D(h_eta_jets, "");
    this->MakeDistribution1D(h_pt_bjets, "TeV");
    this->MakeDistribution1D(h_eta_bjets, "");
    this->MakeDistribution1D(h_pt_qjets, "TeV");
    this->MakeDistribution1D(h_eta_qjets, "");

    this->MakeDistribution1D(h_HT, "TeV");
    this->MakeDistribution1D(h_KT, "TeV");
    this->MakeDistribution1D(h_mvis, "TeV");

    this->MakeDistribution1D(h_HT_all, "TeV");
    this->MakeDistribution1D(h_KT_all, "TeV");
    this->MakeDistribution1D(h_mvis_all, "TeV");

    this->MakeDistribution1D(h_nElectrons, "");
    this->MakeDistribution1D(h_nMuons, "");
    this->MakeDistribution1D(h_nJets, "");

    this->MakeDistribution1D(h_cosTheta1, "");
    this->MakeDistribution1D(h_cosTheta2, "");
    this->MakeDistribution1D(h_cos1cos2, "");
    this->MakeDistribution1D(h_deltaPhi, "rad / #pi");

    this->MakeDistribution1D(h_mtt_tlF, "TeV");
    this->MakeDistribution1D(h_mtt_tlB, "TeV");

    h_AtlFB = this->Asymmetry("AtlFB", "A^{tl}_{FB^*}", h_mtt_tlF, h_mtt_tlB);
    h_AtlFB->GetYaxis()->SetTitle(h_AtlFB->GetTitle());
    h_AtlFB->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AtlFB->GetYaxis()->SetTitleOffset(0.9);
    h_AtlFB->GetXaxis()->SetTitleOffset(0.95);
    h_AtlFB->Write();

    this->MakeDistribution1D(h_mtt_lF, "TeV");
    this->MakeDistribution1D(h_mtt_lB, "TeV");

    h_Ap = this->Asymmetry("Ap", "A_{P}", h_mtt_lF, h_mtt_lB);
    h_Ap->GetYaxis()->SetTitle(h_Ap->GetTitle());
    h_Ap->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_Ap->GetYaxis()->SetTitleOffset(0.9);
    h_Ap->GetXaxis()->SetTitleOffset(0.95);
    h_Ap->Write();

    this->MakeDistribution2D(h2_HT_deltaPhi, "H_{T}", "GeV", "#Delta#phi", "");
    this->MakeDistribution2D(h2_mvis_deltaPhi, "m_{vis}", "GeV", "#Delta#phi", "");
    this->MakeDistribution2D(h2_KT_deltaPhi, "K_{T}", "GeV", "#Delta#phi", "");

    this->MakeDistribution2D(h2_mtt_deltaPhi, "m_{tt}", "GeV", "cos#theta cos#theta", "");
    this->MakeDistribution2D(h2_mtt_cosThetaStar, "m_{tt}", "GeV", "cos#theta^{*}", "");
    this->MakeDistribution2D(h2_mtt_cosTheta1, "m_{tt}", "GeV", "cos#theta_{l^{+}}", "");
    this->MakeDistribution2D(h2_mtt_cosTheta2, "m_{tt}", "GeV", "cos#theta_{l^{-}}", "");

    this->MakeDistributionAL(h2_mtt_cosTheta1, "AL1");
    this->MakeDistributionAL(h2_mtt_cosTheta2, "AL2");
}


void Analysis::MakeDistribution1D(TH1D* h, const string& units) {
    string ytitle, yunits, xunits;
    if (m_xsec) {
        if (m_luminosity > 0) {
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                h->SetBinError(i, sqrt(h->GetBinContent(i)));
                // cout << "N  = " << "" << h->GetBinContent(i) << "\n";
                // cout << "dN = " << "" << h->GetBinError(i) << "\n";
            }
            ytitle = "Expected events";
        }
        else {
            ytitle = "d#sigma / d" + (string) h->GetTitle();
            if (units != "") {
                yunits = " [fb/" + units + "]";
                xunits = " [" + units + "]";
            }
            else {
                yunits = "";
                xunits = "";
            }
        }
    }
    else {
        ytitle = "Generated events";
    }
    if (units != "") xunits = " [" + units + "]";
    else xunits = "";
    h->GetYaxis()->SetTitle((ytitle + yunits).data());
    h->GetXaxis()->SetTitle((h->GetTitle() + xunits).data());
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleOffset(0.95);
    m_output->cd();
    m_output->cd("/");
    h->Write();
}



void Analysis::MakeDistribution2D(TH2D* h, string xtitle, string xunits, string ytitle, string yunits) {
    string ztitle, zunits;
    if (m_xsec) {
        if (m_luminosity > 0) {
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                for (int j = 1; j < h->GetNbinsY() + 1; j++) {
                    h->SetBinError(i, j, sqrt(h->GetBinContent(i, j)));
                }
            }
            ztitle = "Expected events";
        }
        ztitle = "d#sigma / d";
        if (xunits != "" and yunits != "") {
            zunits = " [fb/" + xunits + "/" + yunits + "]";
            xunits = " [" + xunits + "]";
            yunits = " [" + yunits + "]";
        }
        else if (xunits != "" and yunits == "") {
            zunits = " [fb/" + xunits + "]";
            xunits = " [" + xunits + "]";
        }
        else if (xunits == "" and yunits != "") {
            zunits = " [fb/" + yunits + "]";
            yunits = " [" + yunits + "]";
        }
        else {
            zunits = "";
        }
        h->GetZaxis()->SetTitle((ztitle + h->GetTitle() + zunits).data());
    }
    else {
        if (m_luminosity > 0) {
            ztitle = "Expected events";
        }
        else {
            ztitle = "Generated events";
        }
        if (xunits != "") xunits = " [" + xunits + "]";
        if (yunits != "") yunits = " [" + yunits + "]";
        h->GetZaxis()->SetTitle((ztitle).data());
    }
    h->GetYaxis()->SetTitle((ytitle + yunits).data());
    h->GetXaxis()->SetTitle((xtitle + xunits).data());
    m_output->cd();
    m_output->cd("/");
    h->Write();
}

void Analysis::MakeDistributionAL(TH2D* h, const string& name) {
    TF1* func = new TF1("func1", "[0]*x + [1]", -1, 1);
    TObjArray slices;
    this->NormalizeSliceY(h2_mtt_cosTheta2);
    h2_mtt_cosTheta2->FitSlicesY(func, 0, -1, 0, "QRN", &slices);
    for (auto slice : slices) slice->Write();
    TH1D* h_AL = (TH1D*) slices[0]->Clone(name.data());
    slices.Clear();
    h_AL->Scale(2 / h2_mtt_cosTheta2->GetYaxis()->GetBinWidth(1));
    h_AL->SetTitle("A_{L}");
    h_AL->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AL->GetYaxis()->SetTitle("A_{L}");
    m_output->cd();
    m_output->cd("/");
    h_AL->Write();
}


void Analysis::NormalizeSliceY(TH2D* h) {
    double integral = 1;
    int k;
    for (int i = 1; i < h->GetNbinsX() + 1; i++) {
        integral = h->Integral(i, i, 1, h->GetNbinsY());
        for (int j = 1; j < h->GetNbinsY() + 1; j++) {
            k = h->GetBin(i, j);
            h->SetBinContent(k, h->GetBinContent(k) / integral);
            h->SetBinError(k, h->GetBinError(k) / integral);
        }
    }
}


bool Analysis::ExactlyTwoLeptons() {
    bool twoLeptons;
    if (m_electron->size() + m_muon->size() == 2) twoLeptons = true;
    else twoLeptons = false;
    this->UpdateCutflow(c_twoLeptons, twoLeptons);
    return twoLeptons;
}


void Analysis::AssignChannel() {
    if (m_electron->size() == 2) m_channel = "ee";
    else if (m_muon->size() == 2) m_channel = "mumu";
    else if (m_electron->size() == 1 and m_muon->size() == 1) m_channel = "emu";
    else cout << "Error: can't assign channel\n";
    if (m_debug) cout << "Channel assigned: " << m_channel << "\n";
}


bool Analysis::OppositeCharge() {
    double charge1, charge2;
    if (m_channel == "ee") {
        Electron *electron1 = (Electron*) m_electron->at(0);
        Electron *electron2 = (Electron*) m_electron->at(1);

        charge1 = electron1->Charge;
        charge2 = electron2->Charge;
    }
    else if (m_channel == "mumu") {
        Muon *muon1 = (Muon*) m_muon->at(0);
        Muon *muon2 = (Muon*) m_muon->at(1);

        charge1 = muon1->Charge;
        charge2 = muon2->Charge;
    }
    else if (m_channel == "emu") {
        Muon *muon = (Muon*) m_muon->at(0);
        Electron *electron = (Electron*) m_electron->at(0);

        charge1 = muon->Charge;
        charge2 = electron->Charge;
    }
    else cout << "Error: invalid channel\n";

    bool oppositeCharge;

    if (charge1 == charge2) oppositeCharge = false;
    else oppositeCharge = true;
    this->UpdateCutflow(c_oppositeCharge, oppositeCharge);

    return oppositeCharge;
}


bool Analysis::SufficientMll(const pair<TLorentzVector, TLorentzVector>& p_l) {
    // suppress hadronic background, e.g. $J/Psi$
    if (m_debug) cout << "cutting on mll...\n";
    bool sufficientMll;
    double mll = (p_l.first + p_l.second).M();
    if (mll > 15.0) sufficientMll = true;
    else sufficientMll = false;
    this->UpdateCutflow (c_sufficientMll, sufficientMll);
    if (m_debug) cout << "cut on mll\n";
    return sufficientMll;
}


bool Analysis::OutsideZmassWindow(const pair<TLorentzVector, TLorentzVector>& p_l) {
    // suppress DY+jets background
    if (m_debug) cout << "cutting on |mll - mZ| > 10 ...\n";
    bool outsideZmassWindow;
    double mll = (p_l.first + p_l.second).M();
    if ((abs(mll - m_zmass)) > 10.0) outsideZmassWindow = true;
    else outsideZmassWindow = false;
    if (m_channel == "emu") outsideZmassWindow = true;
    this->UpdateCutflow (c_outsideZmassWindow, outsideZmassWindow);
    if (m_debug) cout << "cut on |mll - mZ| > 10 ...\n";
    return outsideZmassWindow;
}


bool Analysis::SufficientBtags() {
    bool sufficientBtags;
    m_btags = 0;
    for (int i = 0; i < m_jet->size(); i++) {
        Jet *jet = (Jet*) m_jet->at(i);
        if (jet->BTag > 0) m_btags++;
    }
    if (m_btags >= m_minimumBtags) sufficientBtags = true;
    else sufficientBtags = false;
    this->UpdateCutflow(c_sufficientBtags, sufficientBtags);
    return sufficientBtags;
}


bool Analysis::SufficientMET() {
    // account for the neutrinos
    // further reduce DY background (no cut for e-mu)
    bool sufficientMET;
    MissingET* missingET = (MissingET*) b_MissingET->At(0);
    double ETmiss = missingET->MET;
    if (ETmiss > 60) sufficientMET = true;
    else sufficientMET = false;
    if (m_channel == "emu") sufficientMET = true;
    this->UpdateCutflow(c_sufficientMET, sufficientMET);
    return sufficientMET;
}


bool Analysis::SufficientHT() {
    // suppress Z/gamma* ( \tau^+ \tau^-) + jets$
    bool sufficientHT;
    ScalarHT *scalarHT = (ScalarHT*) b_ScalarHT->At(0);
    double HT = scalarHT->HT;
    if (HT > 130) sufficientHT = true;
    else sufficientHT = false;
    if (m_channel == "ee" or m_channel == "mumu") sufficientHT = true;
    this->UpdateCutflow(c_sufficientHT, sufficientHT);
    return sufficientHT;
}


bool Analysis::SufficientJets() {
    bool sufficientJets;
    if (m_jet->size() >= 2) sufficientJets = true;
    else sufficientJets = false;
    this->UpdateCutflow(c_sufficientJets, sufficientJets);
    return sufficientJets;
}


void Analysis::PreLoop() {
    this->SetDataDirectory();
    this->SetupInputFiles();
    this->SetupOutputFiles();
    this->InitialiseCutflow();
    this->MakeHistograms();
    cout << "\n";
}


void Analysis::SetDataDirectory() {
    // sets directory based on hostname

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    string Hostname(hostname);

    if (Hostname == "Lorkhan")
        m_dataDirectory = "/Users/declan/Data/zprime/";
    else if ((Hostname.find("lxplus") != string::npos) or (Hostname.find("cern") != string::npos))
        m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime/";
    else if ((Hostname.find("cyan") != string::npos) or (Hostname.find("blue") != string::npos) or (Hostname.find("green") != string::npos))
        m_dataDirectory = "/scratch/dam1g09/zprime/";
    else
        cout << "Hostname " << Hostname << " not recognised.\n";
}


void Analysis::GetBranches() {
    cout << "Fetching branches ...\n";
    b_Jet = m_tree->UseBranch("Jet");
    b_Electron = m_tree->UseBranch("Electron");
    b_Muon = m_tree->UseBranch("Muon");
    b_MissingET = m_tree->UseBranch("MissingET");
    b_ScalarHT = m_tree->UseBranch("ScalarHT");
}


void Analysis::GetGenerationCrossSection(int proc_id) {
    if (m_debug) cout << "Getting generation cross section ...\n";
    string proc_filename = get<0>(m_processes->at(proc_id));
    cout << "Process ID = " << proc_id << "\n";
    cout << "File = " << proc_filename << "\n";

    ifstream proc_file;
    proc_file.open(proc_filename);
    if (!proc_file.is_open()) cout << "Error: Unable to open file.";
    if (m_debug) cout << "m_processes size = " << m_processes->size() << "\n";
    get<3>(m_processes->at(proc_id)) = get_parameter(&proc_file);
    get<4>(m_processes->at(proc_id)) = get_parameter(&proc_file);
    proc_file.close();
    cout << "Generation cross section = " << get<3>(m_processes->at(proc_id)) << " +/- " << get<4>(m_processes->at(proc_id)) << " [fb]\n";
}


void Analysis::GetProcessWeight(int proc_id) {
    int nproc = get<2>(m_processes->at(proc_id));
    if (m_debug) cout << "nproc = " << nproc << "\n";
    int nevents = m_nevents;
    cout << "Number of events = " << nevents << "\n";
    get<5>(m_processes->at(proc_id)) = 1.0 * get<3>(m_processes->at(proc_id)) / (nevents * nproc);
    if (m_luminosity > 0) get<5>(m_processes->at(proc_id)) *= m_luminosity;
    cout << "Event weight = "<< get<5>(m_processes->at(proc_id)) << "\n";
}


void Analysis::Loop() {
    for (itr_s it = m_input->begin(); it != m_input->end(); ++it) {
        int i = it - m_input->begin();
        cout << "Input " << i + 1 << "\n";
        cout << get<0>(*it) << "\n";
        this->EachFile(get<0>(*it));
        m_nevents = this->TotalEvents();
        int proc_id = get<1>(m_input->at(i));
        this->GetGenerationCrossSection(proc_id);
        this->GetProcessWeight(proc_id);
        for (Long64_t jevent = 0; jevent < m_nevents; ++jevent) {
            Long64_t ievent = this->IncrementEvent(jevent);
            if (ievent < 0) break;
            this->EveryEvent(get<5>(m_processes->at(proc_id)));
            this->EachEvent(get<5>(m_processes->at(proc_id)));
            if (!m_debug) ProgressBar(jevent, m_nevents - 1, 50);
        }
        this->CleanUp();
        cout << "\n";
    }
}


void Analysis::EachFile (const string& filename) {
    this->SetupTreesForNewFile(filename);
    this->GetBranches();
}


Analysis::~Analysis() {
    delete m_input;
    delete m_output;
}


Long64_t Analysis::TotalEvents() {
    if (m_tree != 0) return m_tree->GetEntries();
    return -999;
}


Long64_t Analysis::IncrementEvent(Long64_t i) {
    Long64_t ev(-1);
    if (m_tree != 0) ev = m_tree->ReadEntry(i);
    return ev;
}


void Analysis::SetupTreesForNewFile(const string& s) {
    string treeToUse = "Delphes";
    m_chain = new TChain (treeToUse.data(), "");
    string TStringNtuple = s + "/" + treeToUse;
    m_chain->Add(TStringNtuple.data(), 0);
    m_tree = new ExRootTreeReader(m_chain);
    Long64_t numberOfEntries = m_tree->GetEntries();
    // cout << "Number of entries = " << numberOfEntries << "\n";
}


void Analysis::CleanUp() {
    delete m_chain;
    delete m_tree;
}


double Analysis::TotalAsymmetry(TH1D* h_A, TH1D* h_B) {
    double A = h_A->Integral("width");
    double B = h_B->Integral("width");
    double Atot = (A - B) / (A + B);
    return Atot;
}


void Analysis::InitialiseCutflow() {
    m_cutflow = vector<int >(m_cuts, -999);
    m_cutNames = vector<string >(
    m_cuts,                            "no name               ");
    m_cutNames[c_sufficientMET]      = "Sufficient MET        ";
    m_cutNames[c_events]             = "Events                ";
    m_cutNames[c_twoLeptons]         = "Two leptons           ";
    m_cutNames[c_oppositeCharge]     = "Opposite Charge       ";
    m_cutNames[c_sufficientMll]      = "Sufficient mll        ";
    m_cutNames[c_outsideZmassWindow] = "Outside Z mass window ";
    m_cutNames[c_sufficientJets]     = "Sufficient jets       ";
    m_cutNames[c_sufficientBtags]    = "Sufficient b-tags     ";
    m_cutNames[c_sufficientHT]       = "Sufficient HT         ";
    m_cutNames[c_realSolutions]      = "Has real roots        ";
    m_cutNames[c_deltaR]             = "deltaR                ";

    h_cutflow = new TH1D("cutflow", "cutflow", m_cuts, 0.0, m_cuts);
}


void Analysis::UpdateCutflow(const int cut, const bool passed) {
    if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
    if (passed) m_cutflow[cut] += 1;
}


void Analysis::PrintCutflow() {
    cout << "\nCutflow:\n";
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut].c_str());

        cout << m_cutNames[cut] << " " << m_cutflow[cut] << "\n";
    }
    h_cutflow->Write();
}
