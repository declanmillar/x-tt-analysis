#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"
#include "lester_mt2_bisect.h"
#include "neutrino-weighter.hpp"
#include "kinematic-reconstructer.hpp"
#include "two-highest.hpp"
#include "match-bjets-to-leps.hpp"

Analysis::Analysis(const TString& model, const TString& process, const TString& options, const int energy, const int luminosity, const std::string& reconstruction, const TString tag):
    m_model(model),
    m_process(process),
    m_options(options),
    m_energy(energy),
    m_luminosity(luminosity),
    m_reconstruction(reconstruction),
    m_tag(tag),
    m_outputFile(nullptr),
    m_inputFiles(nullptr),
    m_chain(nullptr),
    m_tree(nullptr)
{
    this->PreLoop();
}


void Analysis::Run()
{
    this->Loop();
    this->PostLoop();
}


void Analysis::EachEvent()
{
    UpdateCutflow(c_events, true);

    double weight = 1.0;

    if (!this->SufficientBtags()) return;
    if (!this->TwoLeptons()) return;
    if (!this->OppositeCharge()) return;

    TLorentzVector pvis(0,0,0,0);
    std::vector<TLorentzVector> p_b, p_q, p_j;

    for (int i = 0; i < b_Jet->GetEntries(); i++) {
        Jet *jet = (Jet*) b_Jet->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        if (jet->BTag > 0) {
        // if (jet->BTag & (1 << i)) {
            if (m_debug) std::cout << "b-jet pT = " <<p.Pt() << "\n";
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
        pvis += p;
    }

    std::pair<TLorentzVector, TLorentzVector> p_l;
    if (m_channel == "electron") {
        Electron *electron1 = (Electron*) b_Electron->At(0);
        Electron *electron2 = (Electron*) b_Electron->At(1);

        double charge1 = electron1->Charge;
        double charge2 = electron2->Charge;

        if (charge1 > 0) {
            p_l.first.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
    }
    else if (m_channel == "muon") {
        Muon *muon1 = (Muon*) b_Muon->At(0);
        Muon *muon2 = (Muon*) b_Muon->At(1);

        double charge1 = muon1->Charge;
        double charge2 = muon2->Charge;

        if (charge1 > 0) {
            p_l.first.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
    }
    if (m_debug) p_l.first.Print();
    if (m_debug) p_l.second.Print();

    MissingET *missingET = (MissingET*) b_MissingET->At(0);
    TLorentzVector p_miss;
    double ETmiss = missingET->MET;
    p_miss.SetPtEtaPhiM(ETmiss, missingET->Eta, missingET->Phi, 0.0);

    // ScalarHT *scalarHT = (ScalarHT*) b_ScalarHT->At(0);
    // double HT = scalarHT->HT;
    double HT = p_l.first.Pt() + p_l.second.Pt() + ETmiss;
    for (auto& p : p_b) HT += p.Pt();
    for (auto& p : p_q) HT += p.Pt();

    pvis = pvis + p_l.first + p_l.second;
    for (auto& p : p_b) pvis += p;
    for (auto& p : p_q) pvis += p;
    double mvis = pvis.M();
    double pTvis = pvis.Pt();
    double KT = sqrt(mvis * mvis + pTvis * pTvis) + ETmiss;

    h_HT->Fill(HT / 1000, weight);
    h_mvis->Fill(mvis / 1000, weight);
    h_KT->Fill(KT / 1000, weight);

    if (m_debug) {
        std::cout << "p_l+   = "; p_l.first.Print();
        std::cout << "p_l-   = "; p_l.second.Print();
        for (int i = 0; i < p_b.size(); i++) {
            std::cout << "p_b" << i << "   = ";
            p_b.at(i).Print();
        }
        for (int i = 0; i < p_q.size(); i++) {
            std::cout << "p_q" << i << "   = ";
            p_q.at(i).Print();
        }
        std::cout << "p_miss = "; p_miss.Print();
    }

    std::vector<TLorentzVector> p_bv;
    TLorentzVector p_b1, p_b2, p_v1, p_v2, ttbar, p_top, p_tbar, p_ttbar, p_W1, p_W2;

    std::pair<TLorentzVector, TLorentzVector> p_b_hi = TwoHighestPt(p_b);

    if (m_reconstruction == "KIN") {
        std::vector<TLorentzVector> p_b_hiPt = {p_b_hi.first, p_b_hi.second};
        KinematicReconstructer KIN = KinematicReconstructer(m_bmass, m_Wmass, m_tmass);
        bool isSolution = KIN.Reconstruct(p_l, p_b_hiPt, p_q, p_miss);
        p_top   = KIN.GetTop();
        p_tbar  = KIN.GetTbar();
        p_ttbar = KIN.GetTtbar();
        p_b1    = KIN.GetB();
        p_b2    = KIN.GetBbar();
        p_v1    = KIN.GetNu();
        p_v2    = KIN.GetNubar();
        if (isSolution) {
            this->UpdateCutflow(c_realSolutions, true);
        }
        else {
            this->UpdateCutflow(c_realSolutions, false);
            return;
        }
    }
    else if (m_reconstruction == "NuW") {
        auto p_b_match = MatchBjetsToLeps(p_l, p_b_hi);

        NeutrinoWeighter nuW = NeutrinoWeighter(1, p_l.first.Pt() + p_l.first.Phi()); // 2nd argument is random seed same for specific event
        double weight_max  = nuW.Reconstruct(p_l.first, p_l.second, p_b_match.first, p_b_match.second, p_miss.Px(), p_miss.Py(), p_miss.Phi());
        if (weight_max > 0.0){
            p_top   = nuW.GetTop();
            p_tbar  = nuW.GetTbar();
            p_ttbar = nuW.GetTtbar();
            p_b1    = nuW.GetB();
            p_b2    = nuW.GetBbar();
            p_v1    = nuW.GetNu();
            p_v2    = nuW.GetNubar();
        }

        TLorentzVector p_W1 = p_l.first + p_v1;
        TLorentzVector p_W2 = p_l.second + p_v2;
    }
    else std::cout << "reconstruction = {KIN, NuW}\n";

    double mtt = p_ttbar.M();
    double ytt = p_ttbar.Rapidity();
    double costheta_tt = cos(p_top.Angle(p_tbar.Vect()));
    double delta_abs_yt = std::abs(p_top.Rapidity()) - std::abs(p_tbar.Rapidity());

    // double costheta = pcm_t.CosTheta();
    // double costhetastar = int(ytt / std::abs(ytt)) * costheta;

    // fill histograms
    h_pt_l1->Fill(p_l.first.Pt(), weight);
    h_pt_l2->Fill(p_l.second.Pt(), weight);
    h_eta_l1->Fill(p_l.first.Eta(), weight);
    h_eta_l2->Fill(p_l.second.Eta(), weight);
    h_HT->Fill(HT, weight);

    h_pxt->Fill(p_top.Px(), weight);
    h_pyt->Fill(p_top.Py(), weight);
    h_pzt->Fill(p_top.Pz(), weight);
    h_Et->Fill(p_top.E(), weight);
    h_pTt->Fill(p_top.Pt(), weight);
    h_etat->Fill(p_top.Eta(), weight);
    h_phit->Fill(p_top.Phi() / m_pi, weight);
    h_mt->Fill(p_top.M(), weight);

    h_pxtbar->Fill(p_tbar.Px(), weight);
    h_pytbar->Fill(p_tbar.Py(), weight);
    h_pztbar->Fill(p_tbar.Pz(), weight);
    h_Etbar->Fill(p_tbar.E(), weight);
    h_pTtbar->Fill(p_tbar.Pt(), weight);
    h_etatbar->Fill(p_tbar.Eta(), weight);
    h_phitbar->Fill(p_tbar.Phi() / m_pi, weight);
    h_mtbar->Fill(p_tbar.M(), weight);

    h_mtt->Fill(mtt / 1000, weight);
    h_ytt->Fill(ytt, weight);

    h_mW1->Fill(p_W1.M(), weight);
    h_mW2->Fill(p_W2.M(), weight);
    //     h_costheta_tt->Fill(costheta_tt, weight);
    //     h_cosTheta->Fill(costheta, weight);
    //     h_cosThetaStar->Fill(costhetastar, weight);

    //     if (costhetastar > 0) h_mtt_tF->Fill(mtt, weight);
    //     if (costhetastar < 0) h_mtt_tB->Fill(mtt, weight);

    //     if (delta_abs_yt > 0) h_mtt_tCF->Fill(mtt, weight);
    //     if (delta_abs_yt < 0) h_mtt_tCB->Fill(mtt, weight);

    //     if (n == 6) {

    //         std::vector<TLorentzVector> ptop = p;
    //         v = -1 * (p[0] + p[2] + p[3]).BoostVector();
    //         for (auto& l : ptop) l.Boost(v);

    //         std::vector<TLorentzVector> patop = p;
    //         v = -1 * (p[1] + p[4] + p[5]).BoostVector();
    //         for (auto& l : patop) l.Boost(v);


    //         double costheta_tl1 = cos(ptop[2].Angle(pcm_t.Vect()));
    //         double costheta_tl2 = cos(patop[4].Angle(pcm_tbar.Vect()));
    //         double cos1cos2 = costheta_tl1 * costheta_tl2;
    double deltaPhi = p_l.first.DeltaPhi(p_l.second) / m_pi;
    h_deltaPhi->Fill(deltaPhi, weight);
    //         h_cosTheta1->Fill(costheta_tl1, weight);
    //         h_cosTheta2->Fill(costheta_tl2, weight);
    //         h_cos1cos2->Fill(cos1cos2, weight);
    //         h_deltaR_bW->Fill(deltaR_bW, weight);
    //         h_deltaR_max->Fill(*deltaR_max, weight);

    //         for (int i = 0; i < (int) deltaRs.size(); i++)
    //             h_deltaRs[i]->Fill(deltaRs[i], weight);

    //         for (int i = 0; i < (int) p.size(); i++) {
    //             h_eta[i]->Fill(p[i].Eta(), weight);
    //             h_pt[i]->Fill(p[i].Pt(), weight);
    //         }

    //         if (costheta_tl1 > 0) h_mtt_tlF->Fill(mtt, weight);
    //         if (costheta_tl2 < 0) h_mtt_tlB->Fill(mtt, weight);

    //         if (costheta_tl1 > 0) h_mtt_lF->Fill(mtt, weight);
    //         if (costheta_tl1 < 0) h_mtt_lB->Fill(mtt, weight);

    //         if (phil < m_pi and phil > phi0) h_mtt_philF->Fill(mtt, weight);
    //         if (phil < phi0) h_mtt_philB->Fill(mtt, weight);

    //         if (El > E0) h_mtt_ElF->Fill(mtt, weight);
    //         if (El < E0) h_mtt_ElB->Fill(mtt, weight);

    //         h2_mtt_cosThetaStar->Fill(mtt, costhetastar, weight);
    //         h2_mtt_deltaPhi->Fill(mtt, deltaPhi, weight);
    //         h2_mtt_cosTheta1->Fill(mtt, costheta_tl1, weight);
    //         h2_mtt_cosTheta2->Fill(mtt, costheta_tl2, weight);
    //         h2_mtt_cos1cos2->Fill(mtt, cos1cos2, weight);
    h2_HT_deltaPhi->Fill(HT, deltaPhi, weight);
    h2_mvis_deltaPhi->Fill(mvis, deltaPhi, weight);
    h2_KT_deltaPhi->Fill(KT, deltaPhi, weight);
}


void Analysis::SetupInputFiles()
{
    m_inputFiles = new std::vector<std::string>;
    std::string filename;

    std::string E = std::to_string(m_energy);

    std::vector<std::string> initials = {"gg", "qq", "dd", "uu"};

    int nfiles = 0;
    for (auto initial : initials) {
        std::size_t pos = m_process.find("-tt");
        std::string final = m_process.substr(pos);
        if (boost::contains(m_process, initial)) {
            std::string model = "";
            if (initial == "gg" || initial == "qq") model = "SM";
            else model = m_model;
            std::string options = "";

            options = m_options;
            std::string intermediates = "";
            if (initial == "uu" || initial == "dd") {
                intermediates = intermediates + "-AZ";
                if (m_model != "SM") intermediates = intermediates + "X";
                intermediates = "-X";
            }
            filename = initial + intermediates + final + "_" + model + "_" + E + "TeV" + "_" + m_pdf + options;

            // std::cout << "filename = " << filename << "\n";
            // loop over all matching files (e.g. *.01.root and *.02.root)
            boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
            for (boost::filesystem::directory_iterator i(m_dataDirectory + "/"); i != end_itr; ++i) {
                if (!boost::filesystem::is_regular_file(i->status())) continue;
                std::cout << "file: " << i->path().filename().string() << "\n";
                if (!boost::contains(i->path().filename().string(), filename)) continue;
                std::cout << "file: " << i->path().filename().string() << "\n";
                if (boost::contains(i->path().filename().string(), "KIN")) continue;
                if (boost::contains(i->path().filename().string(), "NuW")) continue;
                if (!boost::contains(i->path().filename().string(), "_pythia_delphes")) continue;
                if (i->path().extension() == ".root") {
                    nfiles++;
                    m_inputFiles->push_back(m_dataDirectory + "/" + i->path().filename().string());
                    std::cout << "Input " << nfiles << ":        " << i->path().filename().string() << "\n";
                }
            }
        }
    }

    // check some input files have been specified
    if (m_inputFiles->size() < 1) {
        std::cout << "Error: no input files specified.";
        exit(false);
    }

    // check all input files exist
    for (auto inputFile : *m_inputFiles) {
        struct stat buffer;
        bool exists = stat(inputFile.c_str(), &buffer) == 0;
        if (exists == false) {
            std::cout << "Error: no " << inputFile << "\n";
            exit(exists);
        }
    }
}


void Analysis::SetupOutputFiles()
{
    std::string E = std::to_string(m_energy);

    m_outputFilename = m_dataDirectory + "/" + m_process + "_" + m_model + "_" + E + "TeV" + "_" + m_pdf + m_options;
    m_outputFilename += "_pythia_delphes";
    m_outputFilename = m_outputFilename + "." + m_reconstruction + m_tag;
    if (m_reconstruction == 2 && m_btags != 2) m_outputFilename += "_b" + std::to_string(m_btags);
    std::string ytt = std::to_string(m_ytt);
    if (m_ytt > 0) m_outputFilename += "_y" + ytt.erase(ytt.find_last_not_of('0') + 1, std::string::npos);
    if (m_Emin >= 0 || m_Emax >= 0) m_outputFilename += "_E" + std::to_string(m_Emin) + "-" + std::to_string(m_Emax);
    std::string eff = std::to_string(m_efficiency);
    if (m_efficiency < 1.0) m_outputFilename += "_e" + eff.erase(eff.find_last_not_of('0') + 1, std::string::npos);
    if (m_luminosity > 0) m_outputFilename += "_L" + std::to_string(m_luminosity);
    if (m_fid) m_outputFilename += "_fid";
    if (m_iso) m_outputFilename += "_iso";
    m_outputFilename += ".root";
    m_outputFile = new TFile(m_outputFilename.c_str(), "RECREATE");
}

void Analysis::PostLoop()
{
    std::cout << "Results\n";
    this->CheckResults();
    if (m_reconstruction == 2) this->CheckPerformance();
    this->MakeDistributions();
    // this->WriteHistograms();
    this->PrintCutflow();
    std::cout << "\n";
    std::cout << "Output\n";
    std::cout << m_outputFilename.c_str() << "\n";
}


void Analysis::CheckResults()
{
    double sigma = h_mtt->Integral("width");
    std::cout << "Analysis cross section: " << sigma << " [fb]\n";
}


void Analysis::CheckPerformance()
{
    double quarkRecoRatio = m_nQuarksMatched / (double) m_nReco;
    double neutrinoRecoRatio = m_nNeutrinoMatched / (double) m_nReco;
    std::cout << "-- Performance --\n";
    printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio * 100);
    printf("pzNu assignment: %.1f%% correct\n", neutrinoRecoRatio * 100);
    std::cout << "\n";
}


TH1D* Analysis::Asymmetry(const TString& name, const TString& title, TH1D* h1, TH1D* h2)
{
    TH1D* h_numerator = (TH1D*) h1->Clone(name);
    TH1D* h_denominator = (TH1D*) h1->Clone();
    h_numerator->SetTitle(title);
    h_numerator->Add(h2, -1);
    h_denominator->Add(h2, 1);
    h_numerator->Divide(h_denominator);
    delete h_denominator;
    if (m_useLumi) this->AsymmetryUncertainty(h_numerator, h1, h2);
    return h_numerator;
}


void Analysis::AsymmetryUncertainty(TH1D* hA, TH1D* h1, TH1D* h2)
{
    double A, dA, N, N1, N2;
    for (int i = 1; i < hA->GetNbinsX() + 1; i++) {
        A = hA->GetBinContent(i);
        N1 = h1->GetBinContent(i);
        N2 = h2->GetBinContent(i);
        N = N1 + N2;
        if (N > 0) dA = sqrt( (1.0 - A * A) / N);
        else dA = 0;
        hA->SetBinError(i, dA);
    }
}


void Analysis::MakeHistograms()
{
    double binWidth = 0.1;
    double Emin = 0.05;
    double Emax = 12.95;
    double nbins = (Emax - Emin) / binWidth;
    std:: cout << "Range:          " << Emin << " -- " << Emax << " [TeV]\n";

    h_pt_l1 = new TH1D("pT_l1", "p^{l^{+}}_{T}", 40, 0, 1000);
    h_pt_l1->Sumw2();
    h_pt_l2 = new TH1D("pT_l2", "p^{l^{-}}_{T}", 40, 0, 1000);
    h_pt_l2->Sumw2();
    h_eta_l1 = new TH1D("eta_l1", "#eta_{l^{+}}", 60, -3, 3);
    h_eta_l1->Sumw2();
    h_eta_l2 = new TH1D("eta_l2", "#eta_{l^{-}}", 60, -3, 3);
    h_eta_l2->Sumw2();

    h_pt_jets = new TH1D("pT_jets", "p_{T}^{jets}", 40, 0, 5000);
    h_pt_jets->Sumw2();
    h_eta_jets = new TH1D("eta_jets", "#eta_{jets}", 60, -3, 3);
    h_eta_jets->Sumw2();
    h_pt_bjets = new TH1D("pT_bjets", "p_{T}^{b-jets}", 40, 0, 5000);
    h_pt_bjets->Sumw2();
    h_eta_bjets = new TH1D("eta_bjets", "#eta_{b-jets}", 60, -3, 3);
    h_eta_bjets->Sumw2();
    h_pt_qjets = new TH1D("pT_qjets", "p_{T}^{q-jets}", 40, 0, 5000);
    h_pt_qjets->Sumw2();
    h_eta_qjets = new TH1D("eta_qjets", "#eta_{q-jets}", 60, -3, 3);
    h_eta_qjets->Sumw2();

    h_mtt = new TH1D("m_tt", "m_{tt}", nbins, Emin, Emax);
    h_mtt->Sumw2();
    h_ytt = new TH1D("y_tt", "y_{tt}", 50, -2.5, 2.5);
    h_ytt->Sumw2();

    h_mW1 = new TH1D("mW1", "m_{W^{+}}", 150, 0, 150);
    h_mW1->Sumw2();
    h_mW2 = new TH1D("mW2", "m_{W^{-}}", 150, 0, 150);
    h_mW2->Sumw2();

    h_pxt = new TH1D("px_t", "p_{x}_{t}", 100, 0, 5000);
    h_pxt->Sumw2();
    h_pyt = new TH1D("py_t", "p_{y}_{t}", 100, 0, 5000);
    h_pyt->Sumw2();
    h_pzt = new TH1D("pz_t", "p_{z}^{t}", 100, 0, 5000);
    h_pzt->Sumw2();
    h_Et = new TH1D("E_t", "E_{t}", 100, 0, 5000);
    h_Et->Sumw2();
    h_pTt = new TH1D("pT_t", "p_{T}^{t}", 100, 0, 5000);
    h_pTt->Sumw2();
    h_etat = new TH1D("eta_t", "#eta_{t}", nbins, 0, 10);
    h_etat->Sumw2();
    h_phit = new TH1D("phi_t", "#phi_{t}", nbins, -1, 1);
    h_phit->Sumw2();
    h_mt = new TH1D("m_t", "m_{t}", 40, 100, 300);
    h_mt->Sumw2();

    h_pxtbar = new TH1D("px_tbar", "p_{x}^{#bar{t}}", 100, 0, 5000);
    h_pxtbar->Sumw2();
    h_pytbar = new TH1D("py_tbar", "p_{y}^{#bar{t}}", 100, 0, 5000);
    h_pytbar->Sumw2();
    h_pztbar = new TH1D("pz_tbar", "p_{z}^{#bar{t}}", 100, 0, 5000);
    h_pztbar->Sumw2();
    h_Etbar = new TH1D("E_tbar", "E_{#bar{t}}", 100, 0, 5000);
    h_Etbar->Sumw2();
    h_pTtbar = new TH1D("pT_tbar", "p_{T}^{#bar{t}}", 100, 0, 5000);
    h_pTtbar->Sumw2();
    h_etatbar = new TH1D("eta_tbar", "#eta_{#bar{t}}", nbins, 0, 10);
    h_etatbar->Sumw2();
    h_phitbar = new TH1D("phi_tbar", "#phi_{#bar{t}}", nbins, -1, 1);
    h_phitbar->Sumw2();
    h_mtbar = new TH1D("m_tbar", "m_{#bar{t}}", 40, 100, 300);
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

    h_HT = new TH1D("HT", "H_{T}", 50, 0, 6);
    h_HT->Sumw2();

    h_KT = new TH1D("KT", "K_{T}", 50, 0, 6);
    h_KT->Sumw2();

    h_mvis = new TH1D("mvis", "m_{vis}", 40, 0, 4);
    h_mvis->Sumw2();

    h_deltaPhi = new TH1D("delta_phi", "#Delta#phi", 10, 0, 1);
    h_deltaPhi->Sumw2();
    h_pv1x = new TH1D("pv1x", "p_{x}^{#nu_{1}}", nbins, -500.0, 500.0);
    h_pv1x->Sumw2();
    h_pv1y = new TH1D("pv1y", "p_{y}^{#nu_{1}}", nbins, -500.0, 500.0);
    h_pv1y->Sumw2();
    h_pv1z = new TH1D("pv1z", "p_{z}^{#nu_{1}}", nbins, -500.0, 500.0);
    h_pv1z->Sumw2();
    h_pv2x = new TH1D("pv2x", "p_{x}^{#nu_{2}}", nbins, -500.0, 500.0);
    h_pv2x->Sumw2();
    h_pv2y = new TH1D("pv2y", "p_{y}^{#nu_{2}}", nbins, -500.0, 500.0);
    h_pv2y->Sumw2();
    h_pv2z = new TH1D("pv2z", "p_{z}^{#nu_{2}}", nbins, -500.0, 500.0);
    h_pv2z->Sumw2();

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

    h2_mtt_deltaPhi = new TH2D("mtt_deltaphi", "m_{tt} #Delta#phi_{l}", nbins, Emin, Emax, 10, 0, 1);
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

    h2_HT_deltaPhi = new TH2D("HT_deltaphi", "H_{T} #Delta#phi_{l}", 40, 0, 4, 10, 0, 1);
    h2_HT_deltaPhi->GetXaxis()->SetTitle("H_{T} [TeV]");
    h2_HT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_HT_deltaPhi->Sumw2();

    h2_mvis_deltaPhi = new TH2D("mvis_deltaphi", "m_{vis} #Delta#phi_{l}", 40, 0, 4, 10, 0, 1);
    h2_mvis_deltaPhi->GetXaxis()->SetTitle("m_{vis} [TeV]");
    h2_mvis_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_mvis_deltaPhi->Sumw2();

    h2_KT_deltaPhi = new TH2D("KT_deltaphi", "K_{T} #Delta#phi_{l}", 40, 0, 4, 10, 0, 1);
    h2_KT_deltaPhi->GetXaxis()->SetTitle("K_{T} [TeV]");
    h2_KT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_KT_deltaPhi->Sumw2();

    h_deltaR_tt = new TH1D("deltaR_tt", "#Delta R(t,#bar{t})", 100, 0, 5);
    h_deltaR_tt->Sumw2();

    h_deltaR_bW = new TH1D("deltaR_bW", "#Delta R(b,W)", 100, 0, 5);
    h_deltaR_bW->Sumw2();
    h_deltaR_max = new TH1D("deltaR_max", "#Delta R_{max}", 100, 0, 5);
    h_deltaR_max->Sumw2();
}


void Analysis::MakeDistributions()
{
    std::cout << "Making distributions ...\n";
    this->MakeDistribution1D(h_mtt, "TeV");
    this->MakeDistribution1D(h_ytt, "");

    this->MakeDistribution1D(h_mW1, "TeV");
    this->MakeDistribution1D(h_mW2, "TeV");

    this->MakeDistribution1D(h_pxt, "GeV");
    this->MakeDistribution1D(h_pyt, "GeV");
    this->MakeDistribution1D(h_pzt, "GeV");
    this->MakeDistribution1D(h_Et, "GeV");
    this->MakeDistribution1D(h_pTt, "GeV");
    this->MakeDistribution1D(h_etat, "");
    this->MakeDistribution1D(h_phit, "");
    this->MakeDistribution1D(h_mt, "GeV");

    this->MakeDistribution1D(h_pxtbar, "GeV");
    this->MakeDistribution1D(h_pytbar, "GeV");
    this->MakeDistribution1D(h_pztbar, "GeV");
    this->MakeDistribution1D(h_Etbar, "GeV");
    this->MakeDistribution1D(h_pTtbar, "GeV");
    this->MakeDistribution1D(h_etatbar, "");
    this->MakeDistribution1D(h_phitbar, "");
    this->MakeDistribution1D(h_mtbar, "GeV");

    // this->MakeDistribution1D(h_cosTheta, "");
    // this->MakeDistribution1D(h_cosThetaStar, "");
    // this->MakeDistribution1D(h_costheta_tt, "");

    // this->MakeDistribution1D(h_mtt_tF, "TeV");
    // this->MakeDistribution1D(h_mtt_tB, "TeV");

    // h_AtFB = this->Asymmetry("AtFB", "A^{t}_{FB^{*}}", h_mtt_tF, h_mtt_tB);
    // h_AtFB->GetYaxis()->SetTitle(h_AtFB->GetTitle());
    // h_AtFB->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_AtFB->Write();

    // this->MakeDistribution1D(h_mtt_tCF, "TeV");
    // this->MakeDistribution1D(h_mtt_tCB, "TeV");

    // h_AtC = this->Asymmetry("AtC", "A^{t}_{C}", h_mtt_tCF, h_mtt_tCB);
    // h_AtC->GetYaxis()->SetTitle(h_AtC->GetTitle());
    // h_AtC->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_AtC->Write();

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

    // this->MakeDistribution1D(h_cosTheta1, "");
    // this->MakeDistribution1D(h_cosTheta2, "");
    // this->MakeDistribution1D(h_cos1cos2, "");
    this->MakeDistribution1D(h_deltaPhi, "rad / #pi");

    // this->MakeDistribution1D(h_mtt_tlF, "TeV");
    // this->MakeDistribution1D(h_mtt_tlB, "TeV");

    // h_AtlFB = this->Asymmetry("AtlFB", "A^{tl}_{FB^*}", h_mtt_tlF, h_mtt_tlB);
    // h_AtlFB->GetYaxis()->SetTitle(h_AtlFB->GetTitle());
    // h_AtlFB->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_AtlFB->GetYaxis()->SetTitleOffset(0.9);
    // h_AtlFB->GetXaxis()->SetTitleOffset(0.95);
    // h_AtlFB->Write();

    // this->MakeDistribution1D(h_mtt_lF, "TeV");
    // this->MakeDistribution1D(h_mtt_lB, "TeV");

    // h_Ap = this->Asymmetry("Ap", "A_{P}", h_mtt_lF, h_mtt_lB);
    // h_Ap->GetYaxis()->SetTitle(h_Ap->GetTitle());
    // h_Ap->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_Ap->GetYaxis()->SetTitleOffset(0.9);
    // h_Ap->GetXaxis()->SetTitleOffset(0.95);
    // h_Ap->Write();

    // h_AlEl = this->Asymmetry("AlEl", "A^{l}_{E_l}", h_mtt_ElF, h_mtt_ElB);
    // h_AlEl->GetYaxis()->SetTitle(h_AlEl->GetTitle());
    // h_AlEl->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_AlEl->GetYaxis()->SetTitleOffset(0.9);
    // h_AlEl->GetXaxis()->SetTitleOffset(0.95);
    // h_AlEl->Write();

    // h_Aphil = this->Asymmetry("Aphil", "A^{l}_{#phi}", h_mtt_philF, h_mtt_philB);
    // h_Aphil->GetYaxis()->SetTitle(h_Aphil->GetTitle());
    // h_Aphil->GetXaxis()->SetTitle("m_{tt} [TeV]");
    // h_Aphil->GetYaxis()->SetTitleOffset(0.9);
    // h_Aphil->GetXaxis()->SetTitleOffset(0.95);
    // h_Aphil->Write();

    // this->MakeDistribution1D(h_pv1x, "GeV");
    // this->MakeDistribution1D(h_pv1y, "GeV");
    // this->MakeDistribution1D(h_pv1z, "GeV");
    // this->MakeDistribution1D(h_pv2x, "GeV");
    // this->MakeDistribution1D(h_pv2y, "GeV");
    // this->MakeDistribution1D(h_pv2z, "GeV");

    // for (auto& h : h_pt) {
    //     h->GetYaxis()->SetTitle("d#sigma / d p_{T}");
    //     h->GetXaxis()->SetTitle("p_{T}");
    // }

    // for (auto& h : h_eta) {
    //     h->GetYaxis()->SetTitle("d#sigma / d #eta");
    //     h->GetXaxis()->SetTitle("#eta");
    // }

    // for (auto& h_deltaR : h_deltaRs) {
    //     h_deltaR->GetYaxis()->SetTitle("d#sigma / d #Delta R");
    //     h_deltaR->GetXaxis()->SetTitle("#Delta R");
    // }

    // this->MakeDistribution1D(h_deltaR_bW, "");
    // this->MakeDistribution1D(h_deltaR_max, "");
    // this->MakeDistribution1D(h_deltaR_tt, "");

    this->MakeDistribution2D(h2_HT_deltaPhi, "H_{T}", "GeV", "#Delta#phi", "");
    this->MakeDistribution2D(h2_mvis_deltaPhi, "m_{vis}", "GeV", "#Delta#phi", "");
    this->MakeDistribution2D(h2_KT_deltaPhi, "K_{T}", "GeV", "#Delta#phi", "");

    //     this->MakeDistribution2D(h2_mtt_deltaPhi, "m_{tt}", "GeV", "cos#theta cos#theta", "");
    //     this->MakeDistribution2D(h2_mtt_cosThetaStar, "m_{tt}", "GeV", "cos#theta^{*}", "");
    //     this->MakeDistribution2D(h2_mtt_cosTheta1, "m_{tt}", "GeV", "cos#theta_{l^{+}}", "");
    //     this->MakeDistribution2D(h2_mtt_cosTheta2, "m_{tt}", "GeV", "cos#theta_{l^{-}}", "");
}


void Analysis::MakeDistribution1D(TH1D* h, const TString& units)
{
    TString ytitle, yunits, xunits;
    if (m_xsec) {
        if (m_useLumi) {
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                h->SetBinError(i, sqrt(h->GetBinContent(i)));
                // std::cout << "N  = " << "" << h->GetBinContent(i) << "\n";
                // std::cout << "dN = " << "" << h->GetBinError(i) << "\n";
            }
            ytitle = "Expected events";
        }
        else {
            ytitle = "d#sigma / d" + (TString) h->GetTitle();
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
    h->GetYaxis()->SetTitle(ytitle + yunits);
    h->GetXaxis()->SetTitle(h->GetTitle() + xunits);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleOffset(0.95);
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::MakeDistribution2D(TH2D* h, TString xtitle,  TString xunits, TString ytitle, TString yunits) {
    TString ztitle, zunits;
    if (m_xsec) {
        // h->Scale(m_sigma / m_nevents, "width");
        if (m_useLumi) {
            // std::cout << "xtitle  = " << xtitle << ", ytitle = " << ytitle << "\n";
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                for (int j = 1; j < h->GetNbinsY() + 1; j++) {
                    h->SetBinError(i, j, sqrt(h->GetBinContent(i, j)));
                    // std::cout << "N  = " << "" << h->GetBinContent(i, j) << ", dN = " << "" << h->GetBinError(i, j) << "\n";
                }
            }
            // std::cout << "\n";
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
        h->GetZaxis()->SetTitle(ztitle + h->GetTitle() + zunits);
    }
    else {
        if (m_useLumi) {
            ztitle = "Expected events";
        }
        else {
            ztitle = "Generated events";
        }
        if (xunits != "") xunits = " [" + xunits + "]";
        if (yunits != "") yunits = " [" + yunits + "]";
        h->GetZaxis()->SetTitle(ztitle);
    }
    h->GetYaxis()->SetTitle(ytitle + yunits);
    h->GetXaxis()->SetTitle(xtitle + xunits);
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::WriteHistograms()
{
    std::cout << "Writing histograms ...\n";
    m_outputFile->cd();
    m_outputFile->cd("/");

    TF1 *func = new TF1("func1", "[0]*x + [1]", -1, 1);
    TObjArray slices1;
    this->NormalizeSliceY(h2_mtt_cosTheta1);
    h2_mtt_cosTheta1->FitSlicesY(func, 0, -1, 0, "QRN", &slices1);
    for (auto slice : slices1) slice->Write();
    TH1D* h_AL1 = (TH1D*) slices1[0]->Clone("AL1");
    slices1.Clear();
    h_AL1->Scale(2 / h2_mtt_cosTheta1->GetYaxis()->GetBinWidth(1));
    h_AL1->SetTitle("A_{L}");
    h_AL1->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AL1->GetYaxis()->SetTitle("A_{L}");
    h_AL1->Write();

    TObjArray slices2;
    this->NormalizeSliceY(h2_mtt_cosTheta2);
    h2_mtt_cosTheta2->FitSlicesY(func, 0, -1, 0, "QRN", &slices2);
    for (auto slice : slices2) slice->Write();
    TH1D* h_AL2 = (TH1D*) slices2[0]->Clone("AL2");
    slices2.Clear();
    h_AL2->Scale(2 / h2_mtt_cosTheta2->GetYaxis()->GetBinWidth(1));
    h_AL2->SetTitle("A_{L}");
    h_AL2->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AL2->GetYaxis()->SetTitle("A_{L}");
    h_AL2->Write();
}

void Analysis::NormalizeSliceY(TH2D* h)
{
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

bool Analysis::PassFiducialCuts(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    if (!m_fid) return true;
    if (this->PassCutsET(p, P)) {
        if (this->PassCutsEta(p, P)) {
            if (!m_iso) return true;
            if (this->PassCutsDeltaR(p, P)) {
                return true;
            }
        }
    }
    return false;
}


bool Analysis::PassCutsMET(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    bool pass = true;
    this->UpdateCutflow(c_MET, pass);
    return pass;
}


bool Analysis::TwoLeptons()
{
    bool twoLeptons;
    if (m_channel == "electron") {
        if (b_Electron->GetEntries() == 2) twoLeptons = true;
        else twoLeptons = false;
    }
    else if (m_channel == "muon") {
        if (b_Muon->GetEntries() == 2) twoLeptons = true;
        else twoLeptons = false;
    }
    this->UpdateCutflow(c_twoLeptons, twoLeptons);
    return twoLeptons;
}


bool Analysis::OppositeCharge()
{
    double charge1, charge2;
    if (m_channel == "electron") {
        Electron *electron1 = (Electron*) b_Electron->At(0);
        Electron *electron2 = (Electron*) b_Electron->At(1);

        charge1 = electron1->Charge;
        charge2 = electron2->Charge;
    }
    else if (m_channel == "muon") {
        Muon *muon1 = (Muon*) b_Muon->At(0);
        Muon *muon2 = (Muon*) b_Muon->At(1);

        charge1 = muon1->Charge;
        charge2 = muon2->Charge;
    }

    bool oppositeCharge;

    if (charge1 == charge2) oppositeCharge = false;
    else oppositeCharge = true;
    this->UpdateCutflow(c_oppositeCharge, oppositeCharge);

    return oppositeCharge;
}


bool Analysis::SufficientBtags()
{
    bool sufficientBtags;
    int nBtags = 0;
    for (int i = 0; i < b_Jet->GetEntries(); i++) {
      Jet *jet = (Jet*) b_Jet->At(i);
      if (jet->BTag > 0) nBtags++;
    }
    if (nBtags >= m_btags) sufficientBtags = true;
    else sufficientBtags = false;
    this->UpdateCutflow(c_sufficientBtags, sufficientBtags);
    return sufficientBtags;
}


bool Analysis::PassCutsMtt(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    if (m_Emin < 0 and m_Emax < 0) return true;

    double mtt = P.M() / 1000;

    if (mtt > m_Emin) {
        if (mtt < m_Emax) {
            UpdateCutflow(c_mtt, true);
            return true;
        }
    }
    UpdateCutflow(c_mtt, false);
    return false;
}


bool Analysis::PassCutsDeltaR(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    std::vector<double> deltaRs;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            if (p[i].DeltaR(p[j]) < 0.4) {
                UpdateCutflow(c_deltaR, false);
                return false;
            }
    UpdateCutflow(c_deltaR, true);
    return true;
}

bool Analysis::PassCutsEta(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    for (unsigned int i = 0; i < p.size(); i++) {
        // bool outsideCrack = abs(p[i].PseudoRapidity()) <= 1.37 || abs(p[i].PseudoRapidity()) >= 1.52;
        bool central = std::abs(p[i].PseudoRapidity()) <= 2.5;
        // bool passesEtaCuts = outsideCrack && central;
        if (central == false) {
            UpdateCutflow(c_eta, false);
            return false;
        }
    }
    UpdateCutflow(c_eta, true);
    return true;
}

bool Analysis::PassCutsET(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    for (unsigned int i = 0; i < p.size(); i++) {
        if(p[i].Et() <= 25) {
            UpdateCutflow(c_Et, false);
            return false;
        }
        else continue;
    }
    UpdateCutflow(c_Et, true);
    return true;
}

bool Analysis::PassCutsYtt(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    if (m_ytt < 0) return true;

    double ytt = std::abs(P.Rapidity());

    if (ytt > m_ytt) {
        UpdateCutflow(c_ytt, true);
        return true;

    }
    UpdateCutflow(c_ytt, false);
    return false;
}


void Analysis::PreLoop()
{
    this->SetDataDirectory();
    this->SetupInputFiles();
    this->SetupOutputFiles();
    this->ResetCounters();
    this->InitialiseCutflow();
    this->MakeHistograms();
    std::cout << "\n";
}


void Analysis::SetDataDirectory()
{
    // sets directory based on hostname

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    std::string Hostname(hostname);

    if (Hostname == "Sunder")
        m_dataDirectory = "/Users/declan/Data/zprime";
    else if ((Hostname.find("lxplus") != std::string::npos) || (Hostname.find("cern") != std::string::npos))
        m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime";
    else if ((Hostname.find("cyan") != std::string::npos) || (Hostname.find("blue") != std::string::npos) || (Hostname.find("green") != std::string::npos))
        m_dataDirectory = "/scratch/dam1g09/zprime";
    else
        std::cout << "Hostname " << Hostname.c_str() << " not recognised.\n";
}


void Analysis::ResetCounters()
{
    if (m_luminosity >= 0) m_useLumi = true;
    else m_useLumi = false;
    m_nQuarksMatched = 0;
    m_nNeutrinoMatched = 0;
    m_nReco = 0;
    m_nRealRoots = 0;
    m_nComplexRoots = 0;
}

void Analysis::GetBranches()
{
    std::cout << "Fetching branches ...\n";
    b_Jet = m_tree->UseBranch("Jet");
    b_Electron = m_tree->UseBranch("Electron");
    b_Muon = m_tree->UseBranch("Muon");
    b_MissingET = m_tree->UseBranch("MissingET");
    b_ScalarHT = m_tree->UseBranch("ScalarHT");
}

void Analysis::GetGenerationCrossSection(TString filename)
{
    std::string Filename = filename.Data();

    std::cout << "Generation Cross section = " << 1000 << " [fb]\n";

    // m_sigma = 1000 * m_crossSection;
}

void Analysis::Loop()
{
    for (itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        std::cout << "Input " << i - m_inputFiles->begin() + 1 << "\n";
        std::cout << (*i) << "\n";
        this->EachFile((*i));
        m_nevents = this->TotalEvents();
        std::cout << "Events: " << m_nevents << "\n";
        for (Long64_t jevent = 0; jevent < m_nevents; ++jevent) {
            Long64_t ievent = this->IncrementEvent(jevent);
            if (ievent < 0) break;
            // std::cout << "event = " << jevent << "\n";
            this->EachEvent();
            if (!m_debug) ProgressBar(jevent, m_nevents - 1, 50);
        }
        this->CleanUp();
        std::cout << "\n";
    }
}

void Analysis::EachFile(TString filename)
{
    this->SetupTreesForNewFile(filename);
    this->GetBranches();
}


Analysis::~Analysis()
{
    delete m_inputFiles;
}


Long64_t Analysis::TotalEvents()
{
    if (m_tree != 0) {return m_tree->GetEntries();}
    return -999;
}


Long64_t Analysis::IncrementEvent(Long64_t i)
{
    Long64_t ev(-1);
    if (m_tree != 0) {ev = m_tree->ReadEntry(i);}
    return ev;
}


void Analysis::SetupTreesForNewFile(const TString& s)
{
    TString treeToUse = "Delphes";
    m_chain = new TChain(treeToUse,"");
    TString TStringNtuple = s + "/" + treeToUse;
    m_chain->Add(TStringNtuple,0);
    m_tree = new ExRootTreeReader(m_chain);
    Long64_t numberOfEntries = m_tree->GetEntries();
    // std::cout << "Number of entries = " << numberOfEntries << "\n";
}


void Analysis::CleanUp()
{
    delete m_chain;
    delete m_tree;
}


double Analysis::TotalAsymmetry(TH1D* h_A, TH1D* h_B) {
    double A = h_A->Integral("width");
    double B = h_B->Integral("width");
    double Atot = (A - B) / (A + B);
    return Atot;
}


void Analysis::InitialiseCutflow()
{
    m_cutflow = std::vector<int>(m_cuts, -999);
    m_cutNames = std::vector<std::string>(
    m_cuts,                         "no name          ");
    m_cutNames[c_events]          = "Events           ";
    m_cutNames[c_twoLeptons]      = "Two leptons      ";
    m_cutNames[c_oppositeCharge]  = "Opposite Charge  ";
    m_cutNames[c_sufficientBtags] = "Sufficient b-tags";
    m_cutNames[c_MET]             = "MET              ";
    m_cutNames[c_realSolutions]   = "Has real roots   ";
    m_cutNames[c_mtt]             = "mtt              ";
    m_cutNames[c_ytt]             = "ytt              ";
    m_cutNames[c_eta]             = "eta              ";
    m_cutNames[c_Et]              = "ET               ";
    m_cutNames[c_deltaR]          = "deltaR           ";

    h_cutflow = new TH1D("cutflow", "cutflow", m_cuts, 0, m_cuts);
}


void Analysis::UpdateCutflow(const int cut, const bool passed)
{
    if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
    if (passed) m_cutflow[cut] += 1;
}


void Analysis::PrintCutflow()
{
    std::cout << "Cutflow: \n";
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut].c_str());

        std::cout << m_cutNames[cut].c_str() << " " << m_cutflow[cut] << "\n";
    }
    h_cutflow->Write();
    m_outputFile->Close();
    delete m_outputFile;
}
