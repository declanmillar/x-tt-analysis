#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"
// #include "delphes-branches.hpp"

Analysis::Analysis(const TString& model, const TString& process, const TString& options, const int energy, const int luminosity, const int reconstruction, const TString tag):
    m_model(model),
    m_process(process),
    m_options(options),
    m_energy(energy),
    m_luminosity(luminosity),
    m_reconstruction(reconstruction),
    m_tag(tag),
    m_outputFile(nullptr),
    m_inputFiles(nullptr),
    m_weightFiles(nullptr),
    m_chain(nullptr),
    m_tree(nullptr)
{
    this->PreLoop();
    this->Loop();
    this->PostLoop();
}


void Analysis::EachEvent()
{
    UpdateCutflow(c_events, true);

    double weight = 1;

    if (!this->TwoElectrons()) return;
    if (!this->OppositeCharge()) return;
    if (!this->SufficientBtags()) return;

    TLorentzVector pvis(0,0,0,0);
    std::vector<TLorentzVector> p_b, p_q, p_j;

    for (int i = 0; i < b_Jet->GetEntries(); i++) {
        Jet *jet = (Jet*) b_Jet->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        if (jet->BTag > 0) {
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

    Electron *electron1 = (Electron*) b_Electron->At(0);
    Electron *electron2 = (Electron*) b_Electron->At(1);

    double charge1 = electron1->Charge;
    double charge2 = electron2->Charge;

    std::pair<TLorentzVector, TLorentzVector> p_l;
    if (charge1 > 0) {
        p_l.first.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0);
        p_l.second.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0);
    }
    else {
        p_l.second.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0);
        p_l.first.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0);
    }

    MissingET *missingET = (MissingET*) b_MissingET->At(0);
    TLorentzVector p_miss;
    double ETmiss = missingET->MET;
    p_miss.SetPtEtaPhiM(ETmiss, missingET->Eta, missingET->Phi, 0);

    ScalarHT *scalarHT = (ScalarHT*) b_ScalarHT->At(0);
    double HT = scalarHT->HT / 1000;

    pvis = pvis + p_l.first + p_l.second;
    for (auto& p : p_b) pvis += p;
    for (auto& p : p_q) pvis += p;
    double mvis = pvis.M();
    double pTvis = pvis.Pt();
    double KT = sqrt(mvis * mvis + pTvis * pTvis) + ETmiss;
    mvis = mvis / 1000;
    KT = KT / 1000;

    h_mvis->Fill(mvis, weight);
    h_KT->Fill(KT, weight);

    if (verbose) {
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
    TLorentzVector b1, b2, v1, v2;

    p_bv = this->ReconstructDilepton(p_l, p_b, p_q, p_miss);

    if (verbose) std::cout << "number of solutions: " << p_bv.size() << "\n";
    if (p_bv.size() == 0) {
        this->UpdateCutflow(c_realSolutions, false);
        return;
    }
    else
        this->UpdateCutflow(c_realSolutions, true);

    b1 = p_bv[0];
    b2 = p_bv[1];
    v1 = p_bv[2];
    v2 = p_bv[3];

    // NeutrinoWeighting* neutrinoWeighting = new NeutrinoWeighting;
    // neutrinoWeighting->apply(p_l[0], p_b[0], p_l[1], p_b[1], p_miss);

    TLorentzVector p_W1 = p_l.first + v1;
    TLorentzVector p_W2 = p_l.second + v2;

    TLorentzVector p_t = b1 + p_W1;
    TLorentzVector p_tbar = b2 + p_W2;

    TLorentzVector P = p_t + p_tbar;
    double mtt = P.M();
    double ytt = P.Rapidity();
    double costheta_tt = cos(p_t.Angle(p_tbar.Vect()));
    double delta_abs_yt = std::abs(p_t.Rapidity()) - std::abs(p_tbar.Rapidity());

    // double costheta = pcm_t.CosTheta();
    // double costhetastar = int(ytt / std::abs(ytt)) * costheta;

    // fill histograms
    h_pt_l1->Fill(p_l.first.Pt(), weight);
    h_pt_l2->Fill(p_l.second.Pt(), weight);
    h_eta_l1->Fill(p_l.first.Eta(), weight);
    h_eta_l2->Fill(p_l.second.Eta(), weight);
    h_HT->Fill(HT, weight);

    h_pxt->Fill(p_t.Px(), weight);
    h_pyt->Fill(p_t.Py(), weight);
    h_pzt->Fill(p_t.Pz(), weight);
    h_Et->Fill(p_t.E(), weight);
    h_pTt->Fill(p_t.Pt(), weight);
    h_etat->Fill(p_t.Eta(), weight);
    h_phit->Fill(p_t.Phi() / m_pi, weight);
    h_mt->Fill(p_t.M(), weight);

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

    //     // if (n < 6) {
    //         // h_mtt_LL->Fill(mtt, m_tree->weightLL());
    //         // h_mtt_LR->Fill(mtt, m_tree->weightLR());
    //         // h_mtt_RL->Fill(mtt, m_tree->weightRL());
    //         // h_mtt_RR->Fill(mtt, m_tree->weightRR());
    //     // }

    //     if (n == 6) {

    //         std::vector<TLorentzVector> ptop = p;
    //         v = -1 * (p[0] + p[2] + p[3]).BoostVector();
    //         for (auto& l : ptop) l.Boost(v);

    //         std::vector<TLorentzVector> patop = p;
    //         v = -1 * (p[1] + p[4] + p[5]).BoostVector();
    //         for (auto& l : patop) l.Boost(v);

    //         double HT = 0;
    //         for (auto& l : p) HT += l.Pt();
    //         HT = HT / 1000;
    //         TLorentzVector pvis = p[0] + p[1] + p[2] + p[4];
    //         double mvis = pvis.M();
    //         double pTvis = pvis.Pt();
    //         double KT = sqrt(mvis * mvis + pTvis * pTvis) + (p[3] + p[5]).Pt();
    //         KT = KT / 1000;

    //         std::vector<double> deltaRs;
    //         for (int i = 0; i < 6; i++)
    //             for (int j = i + 1; j < 6; j++)
    //                 deltaRs.push_back(p[i].DeltaR(p[j]));

    //         double costheta_tl1 = cos(ptop[2].Angle(pcm_t.Vect()));
    //         double costheta_tl2 = cos(patop[4].Angle(pcm_tbar.Vect()));
    //         double cos1cos2 = costheta_tl1 * costheta_tl2;
    double deltaPhi = p_l.first.DeltaPhi(p_l.second) / m_pi;

    //         TLorentzVector p_W = p[2] + p[3];
    //         double deltaR_bW = p_W.DeltaR(p[0]);
    //         std::vector<double> deltaR_ts;
    //         deltaR_ts.push_back(p_t.DeltaR(p[0]));
    //         deltaR_ts.push_back(p_t.DeltaR(p[2]));
    //         deltaR_ts.push_back(p_t.DeltaR(p[3]));
    //         auto deltaR_max = max_element(std::begin(deltaR_ts), std::end(deltaR_ts));

    //         double phil = p[2].Phi();
    //         double El = p[2].E();

    //         double phi0 = 0.7, E0 = 80;

    //         h_KT->Fill(KT, weight);
    //         h_pv1x->Fill(p[3].Px(), weight);
    //         h_pv1y->Fill(p[3].Py(), weight);
    //         h_pv1z->Fill(p[3].Pz(), weight);
    //         h_pv2x->Fill(p[5].Px(), weight);
    //         h_pv2y->Fill(p[5].Py(), weight);
    //         h_pv2z->Fill(p[5].Pz(), weight);
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
    m_inputFiles = new std::vector<TString>;
    m_weightFiles = new std::vector<TString>;
    TString filename;

    std::string E = std::to_string(m_energy);

    std::vector<std::string> initials = {"gg", "qq", "dd", "uu"};

    for (auto initial : initials) {
        std::size_t pos = m_process.find("-tt");
        std::string final = m_process.substr(pos);
        if (boost::contains(m_process, initial)) {
            std::string model = "";
            if (initial == "gg" || initial == "qq") model = "SM";
            else model = m_model;
            std::string options = "";
            if (initial == "gg") options = ".3-5.20x2M.weighted";
            else if (initial == "qq") options = ".3-5.20x2M.weighted";
            options = m_options;
            std::string intermediates = "";
            if (initial == "uu" || initial == "dd") {
                intermediates = intermediates + "-AZ";
                if (m_model != "SM") intermediates = intermediates + "X";
            }
            filename = m_dataDirectory + "/" + initial + intermediates + final + "." + model + "." + E + "TeV" + "." + m_pdf + options;
            filename += ".pythia.delphes";
            m_inputFiles->push_back(filename + ".root");
            m_weightFiles->push_back(filename + ".log");
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
        bool exists = stat(inputFile.Data(), &buffer) == 0;
        if (exists == false) {
            std::cout << "Error: no " << inputFile.Data() << std::endl;
            exit(exists);
        }
    }
}


void Analysis::SetupOutputFiles()
{
    std::string E = std::to_string(m_energy);

    m_outputFilename = m_dataDirectory + "/" + m_process + "." + m_model + "." + E + "TeV" + "." + m_pdf + m_options;
    m_outputFilename += ".pythia.delphes";
    m_outputFilename = m_outputFilename + ".r" + std::to_string(m_reconstruction) + m_tag;
    if (m_reconstruction == 2 && m_btags != 2) m_outputFilename += ".b" + std::to_string(m_btags);
    std::string ytt = std::to_string(m_ytt);
    if (m_ytt > 0) m_outputFilename += ".y" + ytt.erase(ytt.find_last_not_of('0') + 1, std::string::npos);
    if (m_Emin >= 0 || m_Emax >= 0) m_outputFilename += ".E" + std::to_string(m_Emin) + "-" + std::to_string(m_Emax);
    std::string eff = std::to_string(m_efficiency);
    if (m_efficiency < 1.0) m_outputFilename += ".e" + eff.erase(eff.find_last_not_of('0') + 1, std::string::npos);
    if (m_luminosity > 0) m_outputFilename += ".L" + std::to_string(m_luminosity);
    if (m_fid) m_outputFilename += ".fid";
    if (m_iso) m_outputFilename += ".iso";
    m_outputFilename += ".root";
    m_outputFile = new TFile(m_outputFilename, "RECREATE");
}

void Analysis::PostLoop()
{
    this->CheckResults();
    if (m_reconstruction == 2) this->CheckPerformance();
    this->MakeDistributions();
    // this->WriteHistograms();
    this->PrintCutflow();
    std::cout << std::endl << "Output: " << m_outputFilename.Data() << std::endl;
}


void Analysis::CheckResults()
{
    double sigma = h_mtt->Integral("width");
    std::cout << "Analysis cross section: " << sigma << " [fb]" << std::endl;
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
    std:: cout << "Range:          " << Emin << " -- " << Emax << " [TeV]" << std::endl;

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

    h_mW1 = new TH1D("mW1", "m_{W^{+}}", 100, 0, 5000);
    h_mW1->Sumw2();
    h_mW2 = new TH1D("mW2", "m_{W^{-}}", 100, 0, 5000);
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
    std::cout << "Making distributions ..." << std::endl;
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
                // std::cout << "N  = " << "" << h->GetBinContent(i) << std::endl;
                // std::cout << "dN = " << "" << h->GetBinError(i) << std::endl;
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
            // std::cout << "xtitle  = " << xtitle << ", ytitle = " << ytitle << std::endl;
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                for (int j = 1; j < h->GetNbinsY() + 1; j++) {
                    h->SetBinError(i, j, sqrt(h->GetBinContent(i, j)));
                    // std::cout << "N  = " << "" << h->GetBinContent(i, j) << ", dN = " << "" << h->GetBinError(i, j) << std::endl;
                }
            }
            // std::cout << std::endl;
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


void Analysis::WriteHistograms()
{
    std::cout << "Writing histograms ..." << std::endl;
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


bool Analysis::PassCuts(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    // if (this->PassCutsET(p, P)) {
    //     if (this->PassCutsEta(p,P)) {
    //         if (this->PassCutsMET(p, P)) {
            if (this->PassCutsYtt(p, P)) {
                if (this->PassCutsMtt(p, P)) {
                    return true;
                }
            }
    //         }
    //     }
    // }
    return false;
}

bool Analysis::PassCutsMET(const std::vector<TLorentzVector>& p, const TLorentzVector& P)
{
    bool pass = true;
    this->UpdateCutflow(c_MET, pass);
    return pass;
}

bool Analysis::TwoElectrons()
{
    bool twoElectrons;
    if (b_Electron->GetEntries() == 2) twoElectrons = true;
    else twoElectrons = false;
    this->UpdateCutflow(c_twoElectrons, twoElectrons);
    return twoElectrons;
}

bool Analysis::OppositeCharge()
{
    Electron *electron1 = (Electron*) b_Electron->At(0);
    Electron *electron2 = (Electron*) b_Electron->At(1);

    double charge1 = electron1->Charge;
    double charge2 = electron2->Charge;

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
    std::cout << std::endl;
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
    std::cout << "Fetching branches ..." << std::endl;
    b_Jet = m_tree->UseBranch("Jet");
    b_Electron = m_tree->UseBranch("Electron");
    b_MissingET = m_tree->UseBranch("MissingET");
    b_ScalarHT = m_tree->UseBranch("ScalarHT");
}


void Analysis::Loop()
{
    for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        std::cout << "Input : " << i - m_inputFiles->begin() + 1;
        std::cout << (*i) << std::endl;
        this->EachFile((*i));
        m_nevents = this->TotalEvents();
        std::cout << "Events: " << m_nevents << std::endl;
        for (Long64_t jevent = 0; jevent < m_nevents; ++jevent) {
            Long64_t ievent = this->IncrementEvent(jevent);
            if (ievent < 0) break;
            this->EachEvent();
            if (!verbose) ProgressBar(jevent, m_nevents - 1, 50);
        }
        this->CleanUp();
        std::cout << std::endl;
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
    // std::cout << "Number of entries = " << numberOfEntries << std::endl;
}


void Analysis::CleanUp()
{
    delete m_chain;
    delete m_tree;
}


std::vector<TLorentzVector> Analysis::ReconstructSemilepton(const std::vector<TLorentzVector>& p, int Q_l)
{
    // Returns a std::vector of 4-momenta for all 6 particles in the final state with matching of b-quarks to each top
    // and matching of
    // Takes a std::vector of true final-state particle momenta as the argument and the charge of the final
    // state lepton: if +, t decayed leptonically; if -, t~ decayed leptonically.
    // As going from bbllnn->bblnqq/bbqqln requires only a simple re-weighting for parton truth,
    // it saves on storage space and processing time to store all events as bbllnn. However,
    // when we reconstruct the neutrino, we must account for the fact either the top, or the anti-top
    // may decay hadronically. This means there are two distinguishable final states:
    // Q_l = +1 : pp -> b b~ l+ nu q q'
    // Q_l = -1 : pp -> b b~ q q' l- nu
    // Note that the order here is important, as the order of indices in the std::vector of final state momenta
    // relates to the parent particle t=(0,2,3), t~=(1,4,5) and is fixed at the generator level.
    // If we want the results combining each final state, we must add these together.
    // Note: Experimentally p^{x,y}_nu is equated to the MET, of course.

    // this->UpdateCutflow(c_events, true);

    m_nReco++;

    std::vector<TLorentzVector> p_nu_R, p_R(p.size());
    TLorentzVector p_l, p_nu;
    double a, b, c, k, dh, dl, mblv, mjjb, chi2, chi2min = 1.0e10;
    unsigned int imin, jmin;

    // Calculate neutrino p_z solutions
    if (Q_l == +1) {
        // this->UpdateCutflow(c_topDecays, true);
        p_l = p[2];
        p_nu = p[3];
    }
    else if (Q_l == -1) {
        // this->UpdateCutflow(c_antitopDecays, true);
        p_l = p[4];
        p_nu = p[5];
    }
    else {
        std::cout << "Error: Q_l must be ±1.\n";
    }

    double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
    double px_nu = p_nu.Px(), py_nu = p_nu.Py();
    E_l = sqrt(px_l * px_l + py_l * py_l + pz_l * pz_l);
    if (std::abs(E_l - p_l.E()) > 0.00001) std::cout << "ERROR: Lepton energy doesn't match.\n";

    k = 3218.42645 + px_l * px_nu + py_l * py_nu; // m_Wmass * m_Wmass/2 = 3218.42645
    a = px_l * px_l + py_l * py_l;
    b = -2 * k * (pz_l);
    c = (px_nu * px_nu + py_nu * py_nu) * E_l * E_l - k * k;

    double roots[2];
    int nReal = SolveP2(roots, b / a, c / a);
    p_nu_R.clear();
    if (nReal == 2) {
        // two real solutions; pick best match
        // this->UpdateCutflow(c_realSolutions, true);
        for (auto& root : roots) {
            double pz = root;
            TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu * px_nu + py_nu * py_nu + pz * pz));
            p_nu_R.push_back(p);
        }
    }
    else {
        // no real solutions; take the real part
        double pz = roots[0];
        TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu * px_nu + py_nu * py_nu + pz * pz));
        p_nu_R.push_back(p);
    }

    std::vector<TLorentzVector> p_q(4);

    p_q[0] = p[0];
    p_q[1] = p[1];
    if (Q_l == +1) {
        p_q[2] = p[4];
        p_q[3] = p[5];
    }
    else if (Q_l == -1) {
        p_q[2] = p[2];
        p_q[3] = p[3];
    }
    std::vector<std::vector<int> > q_perms;
    if (m_btags == 2) q_perms = { {0, 1, 2, 3}, {1, 0, 2, 3} };
    if (m_btags == 1) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                  {1, 0, 2, 3}, {2, 0, 1, 3}, {3, 0, 1, 2} };
    if (m_btags == 0) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                  {1, 0, 2, 3}, {2, 0, 1, 3}, {3, 0, 1, 2},
                                  {2, 0, 1, 3}, {2, 1, 0, 3}, {2, 3, 0, 1},
                                  {3, 0, 1, 2}, {3, 1, 0, 2}, {3, 2, 0, 1} };

    imin = 0;
    jmin = 0;
    for (unsigned int i = 0; i < p_nu_R.size(); i++) {
        for (unsigned int j = 0; j < q_perms.size(); j++) {
            mblv = (p_q[q_perms[j][0]] + p_l + p_nu_R[i]).M();
            mjjb = (p_q[q_perms[j][1]] + p_q[q_perms[j][2]] + p_q[q_perms[j][3]]).M();
            dh = mjjb - m_tmass;
            dl = mblv - m_tmass;
            chi2 = dh * dh + dl * dl;
            if (chi2 < chi2min) {
                chi2min = chi2;
                imin = i;
                jmin = j;
            }
        }
    }

    if (Q_l == +1) {
        p_R[0] = p_q[q_perms[jmin][0]];
        p_R[1] = p_q[q_perms[jmin][1]];
        p_R[2] = p_l;
        p_R[3] = p_nu_R[imin];
        p_R[4] = p_q[q_perms[jmin][2]];
        p_R[5] = p_q[q_perms[jmin][3]];
    }
    else if (Q_l == -1) {
        p_R[0] = p_q[q_perms[jmin][1]];
        p_R[1] = p_q[q_perms[jmin][0]];
        p_R[2] = p_q[q_perms[jmin][2]];
        p_R[3] = p_q[q_perms[jmin][3]];
        p_R[4] = p_l;
        p_R[5] = p_nu_R[imin];
    }

    // Assess b-matching performance
    unsigned int b_lep;
    if (Q_l == 1) b_lep = 0;
    if (Q_l == -1) b_lep = 1;
    if (b_lep == jmin) m_nQuarksMatched++;

    // Assess neutrino reconstruction performance.
    double pz_nu_truth = p_nu.Pz();
    double Root0MinusTruth;
    double Root1MinusTruth;
    if (nReal == 2) {
        Root0MinusTruth = std::abs(roots[0] - pz_nu_truth);
        Root1MinusTruth = std::abs(roots[1] - pz_nu_truth);
    }
    else if (nReal == 0) {
        Root0MinusTruth = std::abs(roots[0] - pz_nu_truth);
        Root1MinusTruth = std::abs(roots[0] - pz_nu_truth);
    }
    unsigned int bestRoot;
    if (Root0MinusTruth < Root1MinusTruth) bestRoot = 0;
    else if (Root1MinusTruth < Root0MinusTruth) bestRoot = 1;
    else bestRoot = 0;
    if (imin == bestRoot) m_nNeutrinoMatched++;

    if (verbose) {
        // Print reconstruction performance.
        std::cout << "True pz_nu = " << p_nu.Pz() << "\n";
        std::cout << "Possible neutrino solutions:\n";
        if (nReal == 2) {
            printf("                             %f\n", roots[0]);
            printf("                             %f\n", roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        else if (nReal == 0) {
            printf("                             %f + %fi\n", roots[0], roots[1]);
            printf("                             %f - %fi\n", roots[0], roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        if (imin == bestRoot) std::cout << "Neutrino solution: correct.\n";
        else std::cout << "Neutrino solution: incorrect.\n";
        if (b_lep == jmin) std::cout << "b-assignment: correct.\n";
        std::cout << "b-assignment: incorrect.\n";
        std::cout << "--\n";
    }
    return p_R;
}


std::vector<TLorentzVector> Analysis::ReconstructDilepton(const std::pair<TLorentzVector, TLorentzVector>& p_l,
                                                          const std::vector<TLorentzVector>& p_b, const std::vector<TLorentzVector>& p_q,
                                                          const TLorentzVector& p_miss)
{
    // Uses the Sonnenschein method algebraically solve tt dilepton equations.
    // http://arxiv.org/abs/hep-ph/0510100
    // selects solution that minimises mtt

    if (verbose) std::cout << "beginning dilepton reconstruction ...\n";
    std::vector<TLorentzVector> p_bv;

    int n_b = p_b.size();

    if (verbose) std::cout << "number of b-tagged jets = " << n_b << "\n";

    std::vector<std::vector<std::pair<TLorentzVector, TLorentzVector> > > p_vsol;

    // solve quartic for each possible b-jet combination
    for (int i = 0; i < n_b; i++)
        for (int j = 0; j < n_b; j++)
            if (i != j) p_vsol.push_back(this->KinematicReconstruction(p_l.first, p_l.second, p_b[i], p_b[j], p_miss));

    if (verbose) std::cout << "selecting best solution ...\n";
    if (verbose) std::cout << "number b-combinations: " << p_vsol.size() << "\n";

    int nsols = 0;
    for (int i = 0; i < p_vsol.size(); i++) {
        if (verbose) std::cout << "number of real solutions for b-combo " << i << ": " << p_vsol.at(i).size() << "\n";
        nsols += p_vsol.at(i).size();
    }
    if (verbose) std::cout << "total number of solutions = " << nsols << "\n";
    if (nsols == 0) return p_bv;

    // allow no events with null b-combinations
    // bool no_sols = false;
    // for (int i = 0; i < p_vsol.size(); i++) {
    //     if (p_vsol.at(i).size() == 0) no_sols = true;
    // }
    // if (no_sols) return p_bv;

    int I = 0, J = 0, K = 0, L = 0, l = 0;
    double mtt_min = DBL_MAX;
    for (int i = 0; i < n_b; i++) {
        for (int j = 0; j < n_b; j++) {
            if (i != j) {
                for (int k = 0; k < p_vsol.at(l).size(); k++) {
                    double mtt = (p_l.first + p_l.second + p_b[i] + p_b[j] + p_vsol.at(l).at(k).first + p_vsol.at(l).at(k).second).M();
                    if (mtt < mtt_min) {
                        mtt_min = mtt;
                        I = i; J = j; K = k; L = l;
                    }
                }
                l++;
            }
        }
    }

    double mW1 = (p_l.first + p_vsol.at(L).at(K).first).M();
    double mW2 = (p_l.second + p_vsol.at(L).at(K).second).M();
    double mt1 = (p_l.first + p_b[I] + p_vsol.at(L).at(K).first).M();
    double mt2 = (p_l.second + p_b[J] + p_vsol.at(L).at(K).second).M();
    double mtt = (p_l.first + p_l.second + p_b[I] + p_b[J] + p_vsol.at(L).at(K).first + p_vsol.at(L).at(K).second).M();

    if (true) {
        // std::cout << "selected a solution\n";
        if (std::abs(mW1 - m_Wmass) > 10.0) {
            std::cout << "mW+   = " << mW1 << "\n";
            p_l.first.Print();
            p_vsol.at(L).at(K).first.Print();
            std::cout << "number b-combinations: " << p_vsol.size() << "\n";
            for (int i = 0; i < p_vsol.size(); i++) {
                std::cout << "number of real solutions for b-combo " << i << ": " << p_vsol.at(i).size() << "\n";
            }
            std::cout << "total number of solutions = " << nsols << "\n";
            std::cout << "I = " << I << " J = " << J << " K = " << K << " L = "<<  L << "\n";
        }
        // if (std::abs(mW2 - m_Wmass) > 1.0) std::cout << "mW-   = " << mW2 << "\n";
        // if (std::abs(mt1 - m_tmass) > 0.1) std::cout << "mt    = " << mt1 << "\n";
        // if (std::abs(mt2 - m_tmass) > 0.1) std::cout << "mtbar = " << mt2 << "\n";
        // std::cout << "mtt   = " << mtt << "\n";
    }

    p_bv.push_back(p_b.at(I));
    p_bv.push_back(p_b.at(J));
    p_bv.push_back(p_vsol.at(L).at(K).first);
    p_bv.push_back(p_vsol.at(L).at(K).second);

    if (verbose) std::cout << "completed dilepton reconstruction.\n";

    return p_bv;
}

std::vector<std::pair<TLorentzVector, TLorentzVector> > Analysis::KinematicReconstruction(const TLorentzVector& p_l1, const TLorentzVector& p_l2,
                                                                                          const TLorentzVector& p_b1, const TLorentzVector& p_b2,
                                                                                          const TLorentzVector& p_miss)
{
    std::vector<std::pair<TLorentzVector, TLorentzVector> > p_v;

    if (verbose) std::cout << "beginning kinematic reconstruction ...\n";

    // store vector components
    double pl1x = p_l1.Px(), pl1y = p_l1.Py(), pl1z = p_l1.Pz(), El1 = p_l1.E();
    double pl2x = p_l2.Px(), pl2y = p_l2.Py(), pl2z = p_l2.Pz(), El2 = p_l2.E();
    double pb1x = p_b1.Px(), pb1y = p_b1.Py(), pb1z = p_b1.Pz(), Eb1 = p_b1.E();
    double pb2x = p_b2.Px(), pb2y = p_b2.Py(), pb2z = p_b2.Pz(), Eb2 = p_b2.E();
    double Emissx = p_miss.Px(), Emissy = p_miss.Py();

    // store on-shell pole masses
    double mt1 = m_tmass, mt2 = m_tmass;
    double mw1 = m_Wmass, mw2 = m_Wmass;
    double mb1 = m_bmass, mb2 = m_bmass;
    double ml1 = 0.0, ml2 = 0.0;

    double a1 = (Eb1 + El1) * (mw1 * mw1 - ml1 * ml1)
              - El1 * (mt1 * mt1 - mb1 * mb1 - ml1 * ml1)
              + 2 * Eb1 * El1 * El1
              - 2 * El1 * (pb1x * pl1x + pb1y * pl1y + pb1z * pl1z);

    double a2 = 2 * (Eb1 * pl1x - El1 * pb1x);
    double a3 = 2 * (Eb1 * pl1y - El1 * pb1y);
    double a4 = 2 * (Eb1 * pl1z - El1 * pb1z);

    double b1 = (Eb2 + El2) * (mw2 * mw2 - ml2 * ml2)
                - El2 * (mt2 * mt2 - mb2 * mb2 - ml2 * ml2)
                + 2 * Eb2 * El2 * El2
                - 2 * El2 * (pb2x * pl2x + pb2y * pl2y + pb2z * pl2z);

    double b2 = 2 * (Eb2 * pl2x - El2 * pb2x);
    double b3 = 2 * (Eb2 * pl2y - El2 * pb2y);
    double b4 = 2 * (Eb2 * pl2z - El2 * pb2z);

    double c22 = pow(mw1 * mw1 - ml1 * ml1, 2)
                 -4 * (El1 * El1 - pl1z * pl1z) * (a1 / a4) * (a1 / a4)
                 -4 * (mw1 * mw1 - ml1 * ml1) * pl1z * a1 / a4;

    double c21 = 4 * (mw1 * mw1 - ml1 * ml1) * (pl1x - pl1z * a2 / a4)
               -8 * (El1 * El1 - pl1z * pl1z) * a1 * a2 / (a4 * a4)
               -8 * pl1x * pl1z * a1 / a4;

    double c20 = -4 * (El1 * El1 - pl1x * pl1x)
               -4 * (El1 * El1 - pl1z * pl1z) * (a2 / a4) * (a2 / a4)
               -8 * pl1x * pl1z * a2 / a4;

    double c11 = 4 * (mw1 * mw1 - ml1 * ml1) * (pl1y - pl1z * a3 / a4)
               -8 * (El1 * El1 - pl1z * pl1z) * a1 * a3 / (a4 * a4)
               -8 * pl1y * pl1z * a1 /  a4;

    double c10 = -8 * (El1 * El1 - pl1z * pl1z) * a2 * a3 / (a4 * a4)
               +8 * pl1x * pl1y
               -8 * pl1x * pl1z * a3 / a4
               -8 * pl1y * pl1z * a2 / a4;

    double c00 = -4 * (El1 * El1 - pl1y * pl1y)
               -4 * (El1 * El1 - pl1z * pl1z) * (a3 /a4) * (a3 / a4)
               -8 * pl1y * pl1z * a3 / a4;

    c22 = c22 * a4 * a4;
    c21 = c21 * a4 * a4;
    c20 = c20 * a4 * a4;
    c11 = c11 * a4 * a4;
    c10 = c10 * a4 * a4;
    c00 = c00 * a4 * a4;

    double dd22 = pow(mw2 * mw2 - ml2 * ml2, 2)
                  -4 * (El2 * El2 - pl2z * pl2z) * (b1 / b4) * (b1 / b4)
                  -4 * (mw2 * mw2 - ml2 * ml2) * pl2z * b1 / b4;

    double dd21 = 4 * (mw2 * mw2 - ml2 * ml2) * (pl2x - pl2z * b2 / b4)
                  -8 * (El2 * El2 - pl2z * pl2z) * b1 * b2 / (b4 * b4)
                  -8 * pl2x * pl2z * b1 / b4;

    double dd20 = -4 * (El2 * El2 - pl2x * pl2x)
                  -4 * (El2 * El2 - pl2z * pl2z) * (b2 / b4) * (b2 / b4)
                  -8 * pl2x * pl2z * b2 / b4;

    double dd11 = 4 * (mw2 * mw2 - ml2 * ml2) * (pl2y - pl2z * b3 / b4)
                  -8 * (El2 * El2 -pl2z * pl2z) * b1 * b3 / (b4 * b4)
                  -8 * pl2y * pl2z * b1 / b4;

    double dd10 = -8 * (El2 * El2 - pl2z * pl2z) * b2 * b3 / (b4 * b4)
                  +8 * pl2x * pl2y
                  -8 * pl2x * pl2z * b3 / b4
                  -8 * pl2y * pl2z * b2 / b4;

    double dd00 = -4 * (El2 * El2 - pl2y * pl2y)
                  -4 * (El2 * El2 - pl2z * pl2z) * (b3 / b4) * (b3 / b4)
                  -8 * pl2y * pl2z * b3 / b4;

    dd22 = dd22 * b4 * b4;
    dd21 = dd21 * b4 * b4;
    dd20 = dd20 * b4 * b4;
    dd11 = dd11 * b4 * b4;
    dd10 = dd10 * b4 * b4;
    dd00 = dd00 * b4 * b4;

    double d22 = dd22
                 + Emissx * Emissx * dd20
                 + Emissy * Emissy * dd00
                 + Emissx * Emissy * dd10
                 + Emissx * dd21
                 + Emissy * dd11;

    double d21 = -dd21
                 -2 * Emissx * dd20
                 - Emissy * dd10;

    double d20 = dd20;

    double d11 = -dd11
                 -2 * Emissy * dd00
                 - Emissx * dd10;

    double d10 = dd10;
    double d00 = dd00;

    double h4 = c00 * c00 * d22 * d22
                + c11 * d22 * (c11 * d00 - c00 * d11)
                + c00 * c22 * (d11 * d11 - 2 * d00 * d22)
                + c22 * d00 * (c22 * d00 - c11 * d11);

    double h3 = c00 * d21 * (2 * c00 * d22 - c11 * d11)
                + c00 * d11 * (2 * c22 * d10 + c21 * d11)
                + c22 * d00 * (2 * c21 * d00 - c11 * d10)
                - c00 * d22 * (c11 * d10 + c10 * d11)
                -2 * c00 * d00 * (c22 * d21 + c21 * d22)
                - d00 * d11 * (c11 * c21 + c10 * c22)
                + c11 * d00 * (c11 * d21 + 2 * c10 * d22);

    double h2 = c00 * c00 * (2 * d22 * d20 + d21 * d21)
                - c00 * d21 * (c11 * d10 + c10 * d11)
                + c11 * d20 * (c11 * d00 - c00 * d11)
                + c00 * d10 * (c22 * d10 - c10 * d22)
                + c00 * d11 * (2 * c21 * d10 + c20 * d11)
                + (2 * c22 * c20 + c21 * c21) * d00 * d00
                - 2 * c00 * d00 * (c22 * d20 + c21 * d21 + c20 * d22)
                + c10 * d00 * (2 * c11 * d21 + c10 * d22)
                - d00 * d10 * (c11 * c21 + c10 * c22)
                - d00 * d11 * (c11 * c20 + c10 * c21);

    double h1 = c00 * d21 * (2 * c00 * d20 - c10 * d10)
                - c00 * d20 * (c11 * d10 + c10 * d11)
                + c00 * d10 * (c21 * d10 + 2 * c20 * d11)
                - 2 * c00 * d00 * (c21 * d20 + c20 * d21)
                + c10 * d00 * (2 * c11 * d20 + c10 * d21)
                + c20 * d00 * (2 * c21 * d00 - c10 * d11)
                - d00 * d10 * (c11 * c20 + c10 * c21);

    double h0 = c00 * c00 * d20 * d20
                + c10 * d20 * (c10 * d00 - c00 * d10)
                + c20 * d10 * (c00 * d10 - c10 * d00)
                + c20 * d00 * (c20 * d00 - 2 * c00 * d20);

    float a[4] = {static_cast<float>(h1 / h0), static_cast<float>(h2 / h0), static_cast<float>(h3 / h0), static_cast<float>(h4 / h0)};

    if (verbose) std::cout << "solving quartic ...\n";
    double x[4];
    int nReal = SolveP4(x, a[0], a[1], a[2], a[3]);
    if (verbose) std::cout << "found " << nReal << " real solutions.\n";

    double zero_check;
    for (int i = 0; i < nReal; i++) {
        zero_check = pow(x[i], 4) + a[0] * pow(x[i], 3) + a[1] * pow(x[i], 2) + a[2] * x[i] + a[3];
        if (zero_check > 10e-4) std::cout << "warning expression should be zero; solution " << i << " gives: " << zero_check << "\n";
    }

    if (x[0] != x[0] or x[1] != x[1] or x[2] != x[2] or x[3] != x[3]) std::cout << "warning: NaNs in quartic solutions.\n";

    bool keepComplex = false;

    std::vector<double> pv1xs;
    int nSolutions;

    if (keepComplex) {
        if (nReal == 4) {
            // they live in x[0], x[1], x[2], x[3].
            nSolutions = 4;
            for (int i = 0; i < nSolutions; i++) pv1xs.push_back(x[i]);
        }
        else if (nReal == 2) {
            // x[0], x[1] are the real roots and x[2]±i*x[3] are the complex
            nSolutions = 3;
            for (int i = 0; i < nSolutions; i++) pv1xs.push_back(x[i]);
        }
        else if (nReal == 0) {
            // the equation has two pairs of pairs of complex conjugate roots in x[0]±i*x[1] and x[2]±i*x[3].
            nSolutions = 2;
            pv1xs.push_back(x[0]);
            pv1xs.push_back(x[2]);
        }
    }
    else {
        nSolutions = nReal;
        for (int i = 0; i < nReal; i++) pv1xs.push_back(x[i]);
    }

    // Create pairs of neutrino momenta for each real root
    if (verbose) std::cout << "creating neutrino 4-momenta for each real root ...\n";
    for (int i = 0; i < nSolutions; i++) {
        double pv1x = x[i];
        double pv2x = Emissx - pv1x;

        double c2 = c22 + c21 * pv1x + c20 * pv1x * pv1x;
        double c1 = c11 + c10 * pv1x;
        double c0 = c00;
        double d2 = d22 + d21 * pv1x + d20 * pv1x * pv1x;
        double d1 = d11 + d10 * pv1x;
        double d0 = d00;

        double pv1y = (c0 * d2 - c2 * d0) / (c1 * d0 - c0 * d1);
        double pv2y = Emissy - pv1y;

        double pv1z = -(a1 + a2 * pv1x + a3 * pv1y) / a4;
        double pv2z = -(b1 + b2 * pv2x + b3 * pv2y) / b4;

        double Ev1 = sqrt(pv1x * pv1x + pv1y * pv1y + pv1z * pv1z);
        double Ev2 = sqrt(pv2x * pv2x + pv2y * pv2y + pv2z * pv2z);

        TLorentzVector pv1(pv1x, pv1y, pv1z, Ev1);
        TLorentzVector pv2(pv2x, pv2y, pv2z, Ev2);

        std::pair<TLorentzVector, TLorentzVector> pv(pv1, pv2);

        p_v.push_back(pv);
    }

    if (verbose) std::cout << "found " << p_v.size() << " pairs of neutrinos.\n";

    return p_v;
}


double Analysis::TotalAsymmetry(TH1D* h_A, TH1D* h_B) {
    double A = h_A->Integral("width");
    double B = h_B->Integral("width");
    double Atot = (A - B)/(A + B);
    return Atot;
}


void Analysis::InitialiseCutflow()
{
    m_cutflow = std::vector<int>(m_cuts, -999);
    m_cutNames = std::vector<TString>(
    m_cuts,                         "no name          ");
    m_cutNames[c_events]          = "Events           ";
    m_cutNames[c_twoElectrons]    = "Two electrons    ";
    m_cutNames[c_oppositeCharge]  = "Opposite Charge  ";
    m_cutNames[c_sufficientBtags] = "Sufficient b-tags";
    m_cutNames[c_MET]             = "MET              ";
    m_cutNames[c_realSolutions]   = "pvz real         ";
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
    std::cout << "Cutflow: " << std::endl;
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut]);

        std::cout << m_cutNames[cut].Data() << " " << m_cutflow[cut] << "\n";
    }
    h_cutflow->Write();
    m_outputFile->Close();
    delete m_outputFile;
}
