#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"

Analysis::Analysis(const TString& model, const TString& initial_state, const TString& intermediates, const TString& final_state,  const int energy, const TString& options, const int vegasIterations, const string vegasPoints, const bool add_gg, const bool add_qq, const int luminosity, const int btags, const TString& tag):
    m_model(model),
    m_initial_state(initial_state),
    m_intermediates(intermediates),
    m_channel(final_state),
    m_energy(energy),
    m_options(options),
    m_vegasIterations(vegasIterations),
    m_vegasPoints(vegasPoints),
    m_add_gg(add_gg),
    m_add_qq(add_qq),
    m_luminosity(luminosity),
    m_efficiency(1.0),
    m_btags(btags),
    m_tag(tag),
    m_reco(1), 
    m_debug(false),
    m_outputFile(nullptr),
    m_inputFiles(nullptr),
    m_weightFiles(nullptr),
    m_chainNtup(nullptr),
    m_ntup(nullptr)
{ 
    ;
}


void Analysis::Run()
{
    this->PreLoop();
    this->Loop();
    this->PostLoop();
}


void Analysis::EachEvent()
{
    // UpdateCutflow(c_entries, true);

    std::vector<TLorentzVector> p(6);
    for (int i = 0; i < 6; i++)
        p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));

    TLorentzVector P(0, 0, 0, 0);
    for (auto& l : p) P += l;

    std::vector<TLorentzVector> pcm = p;
    TVector3 v = -1 * P.BoostVector();
    for (auto& l : pcm) l.Boost(v);

    std::vector<TLorentzVector> ptop = p;
    v = -1 * (p[0] + p[2] + p[3]).BoostVector();
    for (auto& l : ptop) l.Boost(v);

    std::vector<TLorentzVector> patop = p;
    v = -1 * (p[1] + p[4] + p[5]).BoostVector();
    for (auto& l : patop) l.Boost(v);

    std::vector<TLorentzVector> p_R1(6), pcm_R1(6), ptop_R1(6), patop_R1(6);
    std::vector<TLorentzVector> p_R2(6), pcm_R2(6), ptop_R2(6), patop_R2(6);
    TLorentzVector P_R1, P_R2;
    if (m_reco == 1) {
        p_R1 = this->ReconstructDilepton(p); // both decay leptonically

        P_R1.SetPxPyPzE(0, 0, 0, 0);
        for (auto& l : p_R1) P_R1 += l;

        pcm_R1 = p_R1;
        v = -1 * P_R1.BoostVector();
        for (auto& l : pcm_R1) l.Boost(v);

        ptop_R1 = p_R1;
        v = -1 * (p_R1[0] + p_R1[2] + p_R1[3]).BoostVector();
        for (auto& l : ptop_R1) l.Boost(v);

        patop_R2 = p_R1;
        v = -1 * (p_R1[1] + p_R1[4] + p_R1[5]).BoostVector();
        for (auto& l : patop_R2) l.Boost(v);
    }
    else if (m_reco == 2) {
        p_R1 = this->ReconstructSemilepton(p, +1); // top decays leptonically
        p_R2 = this->ReconstructSemilepton(p, -1); // antitop decays leptonically

        P_R1.SetPxPyPzE(0, 0, 0, 0);
        for (auto& l : p_R1) P_R1 += l;

        P_R2.SetPxPyPzE(0, 0, 0, 0);
        for (auto& l : p_R2) P_R2 += l;

        pcm_R1 = p_R1;
        v = -1 * P_R1.BoostVector();
        for (auto& l : pcm_R1) l.Boost(v);

        pcm_R2 = p_R2;
        v = -1 * P_R2.BoostVector();
        for (auto& l : pcm_R2) l.Boost(v);

        ptop_R1 = p_R1;
        v = -1 * (p_R1[0] + p_R1[2] + p_R1[3]).BoostVector();
        for (auto& l : ptop_R1) l.Boost(v);

        patop_R2 = p_R2;
        v = -1 * (p_R2[1] + p_R2[4] + p_R2[5]).BoostVector();
        for (auto& l : patop_R2) l.Boost(v);
    }

    // top and antitop
    TLorentzVector p_t = p[0] + p[2] + p[3];
    TLorentzVector p_tb = p[1] + p[4] + p[5];
    TLorentzVector pcm_t = pcm[0] + pcm[2] + pcm[3];
    TLorentzVector pcm_tb = pcm[1] + pcm[4] + pcm[5];
    TLorentzVector p_W = p[2] + p[3];
    TLorentzVector pcm_t_R1;
    TLorentzVector pcm_tb_R1;
    TLorentzVector pcm_t_R2;
    TLorentzVector pcm_tb_R2;

    if (m_reco > 0) {
        pcm_t_R1 = pcm_R1[0] + pcm_R1[2] + pcm_R1[3];
        pcm_tb_R1 = pcm_R1[1] + pcm_R1[4] + pcm_R1[5];
    }
    if (m_reco == 2) {
        pcm_t_R2 = pcm_R2[0] + pcm_R2[2] + pcm_R2[3];
        pcm_tb_R2 = pcm_R2[1] + pcm_R2[4] + pcm_R2[5];
    }


    double mtt = P.M() / 1000;
    double mtt_R1 = P_R1.M() / 1000;
    double mtt_R2 = P_R2.M() / 1000;

    double HT = 0;
    for (auto& l : p) HT += l.Pt();
    HT = HT / 1000;
    TLorentzVector pvis = p[0] + p[1] +p[2] + p[4];
    double mvis = pvis.M();
    double pTvis = pvis.Pt();
    double KT = sqrt(mvis * mvis + pTvis * pTvis) + (p[3] + p[5]).Pt();
    KT = KT / 1000;

    double ytt = P.Rapidity();
    double cosTheta = pcm_t.CosTheta();
    double cosThetaStar = int(ytt / std::abs(ytt)) * cosTheta;

    std::vector<double> deltaRs;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            deltaRs.push_back(p[i].DeltaR(p[j]));

    double cosTheta1 = cos(ptop[2].Angle(pcm_t.Vect()));
    double cosTheta2 = cos(patop[4].Angle(pcm_tb.Vect()));
    double cos1cos2 = cosTheta1 * cosTheta2;
    double deltaPhi = p[2].DeltaPhi(p[4]) / m_pi;

    double deltaRbW = p_W.DeltaR(p[0]);
    std::vector<double> deltaRt;
    deltaRt.push_back(p_t.DeltaR(p[0]));
    deltaRt.push_back(p_t.DeltaR(p[2]));
    deltaRt.push_back(p_t.DeltaR(p[3]));
    auto deltaRmax = max_element(std::begin(deltaRt), std::end(deltaRt));

    double ytt_R1, ytt_R2;
    double cosTheta_R1, cosTheta_R2;
    double cosThetaStar_R1, cosThetaStar_R2;
    double cosTheta1_R1, cosTheta2_R2;

    if (m_reco > 0) {
        ytt_R1 = P_R1.Rapidity();
        cosTheta_R1 = pcm_t_R1.CosTheta();
        cosThetaStar_R1 = int(ytt_R1 / std::abs(ytt_R1)) * cosTheta_R1;
        cosTheta1_R1 = cos(ptop_R1[2].Angle(pcm_t_R1.Vect()));
    }
    if (m_reco == 2) {
        ytt_R2 = P_R2.Rapidity();
        cosTheta_R2 = pcm_t_R2.CosTheta();
        cosThetaStar_R2 = int(ytt_R2 / std::abs(ytt_R2)) * cosTheta_R2;
        cosTheta2_R2 = cos(patop_R2[4].Angle(pcm_tb_R2.Vect()));
    }

    double weight = 1, weight_R = 1;
    if (m_xsec) {
        double it = m_ntup->iteration();
        weight = m_ntup->weight();
        weight = 1000 * weight * m_sigma/iteration_weights[it - 1];
        if (m_reco == 2) weight_R = weight / 2;
        else weight_R = weight;
    }

    // if (this->PassCuts("truth")) {
        h_mt->Fill(p_t.M(), weight);
        h_mtbar->Fill(p_tb.M(), weight);
        h_mtt->Fill(mtt, weight);
        h_ytt->Fill(ytt, weight);
        h_cosTheta->Fill(cosTheta, weight);
        h_cosThetaStar->Fill(cosThetaStar, weight);
        h_HT->Fill(HT, weight);
        h_KT->Fill(KT, weight);
        h_pv1x->Fill(p[3].Px(), weight);
        h_pv1y->Fill(p[3].Py(), weight);
        h_pv1z->Fill(p[3].Pz(), weight);
        h_pv2x->Fill(p[5].Px(), weight);
        h_pv2y->Fill(p[5].Py(), weight);
        h_pv2z->Fill(p[5].Pz(), weight);
        h_deltaPhi->Fill(deltaPhi, weight);
        h_cosTheta1->Fill(cosTheta1, weight);
        h_cosTheta2->Fill(cosTheta2, weight);
        h_cos1cos2->Fill(cos1cos2, weight);
        h_deltaRbW->Fill(deltaRbW, weight);
        h_deltaRmax->Fill(*deltaRmax, weight);

        for (int i = 0; i < (int) deltaRs.size(); i++)
            h_deltaRs[i]->Fill(deltaRs[i], weight);

        for (int i = 0; i < (int) p.size(); i++) {
            h_eta[i]->Fill(p[i].Eta(), weight);
            h_pt[i]->Fill(p[i].Pt(), weight);
        }

        if (cosThetaStar > 0) h_mtt_F->Fill(mtt, weight);
        if (cosThetaStar < 0) h_mtt_B->Fill(mtt, weight);

        if (cosTheta1 > 0) h_mtt_Fl->Fill(mtt, weight);
        if (cosTheta1 < 0) h_mtt_Bl->Fill(mtt, weight);

        h2_mtt_cosThetaStar->Fill(mtt, cosThetaStar, weight);
        h2_mtt_deltaPhi->Fill(mtt, deltaPhi, weight);
        h2_mtt_cosTheta1->Fill(mtt, cosTheta1, weight);
        h2_mtt_cosTheta2->Fill(mtt, cosTheta2, weight);
        h2_mtt_cos1cos2->Fill(mtt, cos1cos2, weight);
        h2_HT_deltaPhi->Fill(HT, deltaPhi, weight);
        h2_KT_deltaPhi->Fill(KT, deltaPhi, weight);
    // }

    if (m_reco > 0) {
        // if (this->PassCuts("R1")) {
            h_mtt_R->Fill(mtt_R1, weight_R);
            h_ytt_R->Fill(ytt_R1, weight_R);
            h_mt_R->Fill(pcm_t_R1.M(), weight_R);
            h_mtbar_R->Fill(pcm_tb_R1.M(), weight_R);
            h_cosTheta_R->Fill(cosTheta_R1, weight_R);
            h_cosThetaStar_R->Fill(cosThetaStar_R1, weight_R);
            h_pv1x_R->Fill(p_R1[3].Px(), weight_R);
            h_pv1y_R->Fill(p_R1[3].Py(), weight_R);
            h_pv1z_R->Fill(p_R1[3].Pz(), weight_R);
            h_pv2x_R->Fill(p_R1[5].Px(), weight_R);
            h_pv2y_R->Fill(p_R1[5].Py(), weight_R);
            h_pv2z_R->Fill(p_R1[5].Pz(), weight_R);
            if (cosThetaStar_R1 > 0) h_mtt_FR->Fill(mtt_R1, weight_R);
            if (cosThetaStar_R1 < 0) h_mtt_BR->Fill(mtt_R1, weight_R);
            h2_mtt_cosThetaStar_R->Fill(mtt_R1, cosThetaStar_R1, weight_R);
            h2_mtt_cosThetal_R->Fill(mtt_R1, cosTheta1_R1, weight_R);
        // }
    }
    if (m_reco == 2) {
        // if (this->PassCuts("R2")) {
            h_mtt_R->Fill(mtt_R2, weight_R);
            h_ytt_R->Fill(ytt_R2, weight_R);
            h_mt_R->Fill(pcm_t_R2.M(), weight_R);
            h_mtbar_R->Fill(pcm_tb_R2.M(), weight_R);
            h_cosTheta_R->Fill(cosTheta_R2, weight_R);
            h_cosThetaStar_R->Fill(cosThetaStar_R2, weight_R);
            if (cosThetaStar_R2 > 0) h_mtt_FR->Fill(mtt_R2, weight_R);
            if (cosThetaStar_R2 < 0) h_mtt_BR->Fill(mtt_R2, weight_R);
            h2_mtt_cosThetaStar_R->Fill(mtt_R2, cosThetaStar_R2, weight_R);
            h2_mtt_cosThetal_R->Fill(mtt_R2, cosTheta2_R2, weight_R);
        // }
    }
}


void Analysis::SetupInputFiles()
{
    m_inputFiles = new std::vector<TString>;
    m_weightFiles = new std::vector<TString>;
    TString filename;

    string E = "";
    if (m_energy != 13) "_" + std::to_string(m_energy);

    if (m_add_gg) {
      filename = m_dataDirectory + "/SM_" + "gg-G-" + m_channel + E + "_2-4_" + std::to_string(m_vegasIterations) + "x" + m_vegasPoints;
      m_inputFiles->push_back(filename + ".root");
      m_weightFiles->push_back(filename + ".log");
    }

    if (m_add_qq) {
      filename = m_dataDirectory + "/SM_" + "qq-G-" + m_channel + E + "_2-4_" + std::to_string(m_vegasIterations) + "x" + m_vegasPoints;
      m_inputFiles->push_back(filename + ".root");
      m_weightFiles->push_back(filename + ".log");
    }

    filename = m_dataDirectory + "/" + m_model + "_" + m_initial_state + "-" + m_intermediates + m_channel + E + m_options + std::to_string(m_vegasIterations) + "x" + m_vegasPoints;
    m_inputFiles->push_back(filename + ".root");
    m_weightFiles->push_back(filename + ".log");

    // Check all input files exist
    for (auto inputFile : *m_inputFiles) {
        bool exists = false;
        struct stat buffer;
        exists = stat(inputFile.Data(), &buffer) == 0;
        if (exists == false) {
            printf("Error: %s does not exist.\n", inputFile.Data());
            exit(exists);
        }
    }
}


void Analysis::SetupOutputFiles()
{
    TString outfilename;
    TString initial_state = m_initial_state;
    TString intermediates = m_intermediates;
    string E = "";
    if (m_energy != 13) "_" + std::to_string(m_energy);

    if (m_add_gg || m_add_qq) intermediates = "G" + intermediates;
    if (m_add_gg) initial_state = "gg" + initial_state;
    outfilename = m_dataDirectory + "/" + m_model + "_" + initial_state + "-" + intermediates + m_channel + E + m_options + std::to_string(m_vegasIterations) + "x" + m_vegasPoints;
    m_outputFilename = outfilename + ".a";

    if (m_reco == 2 && m_btags != 2) m_outputFilename += ".b" + std::to_string(m_btags);
    string ytt = std::to_string(m_ytt);
    if (m_ytt > 0) m_outputFilename += ".y" + ytt.erase(ytt.find_last_not_of('0') + 1, string::npos);
    if (m_Emin >= 0 || m_Emax >= 0) m_outputFilename += ".E" + std::to_string(m_Emin) + "-" + std::to_string(m_Emax);
    string eff = std::to_string(m_efficiency);
    if (m_efficiency < 1.0) m_outputFilename += ".e" + eff.erase(eff.find_last_not_of('0') + 1, string::npos);
    if (m_luminosity > 0) m_outputFilename += ".L" + std::to_string(m_luminosity);
    if (m_fid == true) m_outputFilename += ".fid";
    m_outputFilename += m_tag;
    m_outputFilename += ".root";

    printf("-- Output --\n");
    printf("%s\n", m_outputFilename.Data());
    m_outputFile = new TFile(m_outputFilename, "RECREATE");
}

void Analysis::PostLoop()
{
    this->CheckResults();
    if (m_reco == 2) this->CheckPerformance();
    this->MakeDistributions();
    // this->PrintCutflow();
    this->WriteHistograms();
}


void Analysis::CheckResults()
{
    printf("-- Results --\n");
    double sigma = h_mtt->Integral("width");
    printf("Analysis Cross Section = %.15le [pb]\n", sigma);
    printf("\n");
}


void Analysis::CheckPerformance()
{
    double quarkRecoRatio = m_nQuarksMatched/(double)m_nReco;
    double neutrinoRecoRatio = m_nNeutrinoMatched/(double)m_nReco;
    printf("-- Performance --\n");
    printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio * 100);
    printf("pzNu assignment: %.1f%% correct\n", neutrinoRecoRatio * 100);
    printf("\n");
}


void Analysis::ApplyLuminosity(TH1D* h) 
{
    double sigma, N, dN;
    for (int i = 1; i < h->GetNbinsX() + 1; i++) {
        sigma = h->GetBinContent(i);
        N = m_luminosity * m_efficiency * sigma;
        h->SetBinContent(i, N);
        dN = sqrt(N);
        h->SetBinError(i, dN);
    }
    h->GetYaxis()->SetTitle("Expected events");
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
    if (m_luminosity > -1) this->AsymmetryUncertainty(h_numerator, h1, h2);
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
        if (N > 0) dA = sqrt((1.0 - A * A) / N);
        else dA = 0;
        hA->SetBinError(i, dA);
    }
}


void Analysis::MakeHistograms()
{
    double binWidth = 0.05;
    double Emin = 0;
    double Emax = 13;
    double nbins = (Emax - Emin) / binWidth;

    h_mtt = new TH1D("mtt", "m_{tt}", nbins, Emin, Emax);
    h_mtt->Sumw2();
    h_ytt = new TH1D("ytt", "y_{tt}", 50, -2.5, 2.5);
    h_ytt->Sumw2();
    h_mt = new TH1D("mt", "m_{t}", nbins, 0, 350);
    h_mt->Sumw2();
    h_mtbar = new TH1D("mtbar", "m_{#bar{t}}", nbins, 0, 350);
    h_mtbar->Sumw2();
    h_mtt_F = new TH1D("mtt_F", "m_{tt}^{forward}", nbins, Emin, Emax);
    h_mtt_F->Sumw2();
    h_mtt_B = new TH1D("mtt_B", "m_{tt}^{backward}", nbins, Emin, Emax);
    h_mtt_B->Sumw2();
    h_cosTheta = new TH1D("costheta", "cos#theta", nbins, -1.0, 1.0);
    h_cosTheta->Sumw2();
    h_cosThetaStar = new TH1D("costhetastar", "cos#theta^{*}", nbins, -1.0, 1.0);
    h_cosThetaStar->Sumw2();

    h_HT = new TH1D("HT", "H_{T}", nbins, 0, 4);
    h_HT->Sumw2();

    h_KT = new TH1D("KT", "K_{T}", nbins, 0, 4);
    h_KT->Sumw2();

    h_deltaPhi = new TH1D("deltaphi", "#Delta#phi", 10, 0, 1);
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
    h_cosTheta1 = new TH1D("costheta1", "cos#theta_{l+}", 10, -1.0, 1.0);
    h_cosTheta1->Sumw2();
    h_cosTheta2 = new TH1D("costheta2", "cos#theta_{l-}", 10, -1.0, 1.0);
    h_cosTheta2->Sumw2();
    h_cos1cos2 = new TH1D("cos1cos2", "cos#theta_{l+}cos#theta_{l-}", 20, -1.0, 1.0);
    h_cos1cos2->Sumw2();

    h_mtt_Fl = new TH1D("mtt_Fl", "m_{tt}^{F,l}", 19, 2.05, 3.95);
    h_mtt_Fl->Sumw2();
    h_mtt_Bl = new TH1D("mtt_Bl", "m_{tt}^{B,l}", 19, 2.05, 3.95);
    h_mtt_Bl->Sumw2();

    h2_mtt_cosThetaStar = new TH2D("mtt_costhetastar", "m_{tt} cos#theta^{*}", nbins, Emin, Emax, 2, -1.0, 1.0);
    h2_mtt_cosThetaStar->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_cosThetaStar->GetYaxis()->SetTitle("cos#theta^*");
    h2_mtt_cosThetaStar->Sumw2();

    h2_mtt_cosThetaStar_R = new TH2D("mtt_costhetastar_r", "m_{tt} cos#theta^{*} (reco)", nbins, Emin, Emax, 2, -1.0, 1.0);
    h2_mtt_cosThetaStar_R->GetXaxis()->SetTitle("m_{tt} (reco)");
    h2_mtt_cosThetaStar_R->GetYaxis()->SetTitle("cos#theta^* (reco)");
    h2_mtt_cosThetaStar_R->Sumw2();

    h2_mtt_deltaPhi = new TH2D("mtt_deltaphi", "m_{tt} #Delta#phi_{l}", nbins, Emin, Emax, 10, 0, 1);
    h2_mtt_deltaPhi->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l}");
    h2_mtt_deltaPhi->Sumw2();

    h2_mtt_cosTheta1 = new TH2D("mtt_costheta1", "m_{tt} cos#theta_{l+}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta1->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h2_mtt_cosTheta1->GetYaxis()->SetTitle("cos#theta_{l+}");
    h2_mtt_cosTheta1->Sumw2();

    h2_mtt_cosTheta2 = new TH2D("mtt_costheta2", "m_{tt} cos#theta_{l-}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta2->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h2_mtt_cosTheta2->GetYaxis()->SetTitle("cos#theta_{l-}");
    h2_mtt_cosTheta2->Sumw2();

    h2_mtt_cosThetal_R = new TH2D("mtt_costhetal_R", "m_{tt} cos#theta_{l} (reco)", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosThetal_R->GetXaxis()->SetTitle("m_{tt} (reco) [TeV]");
    h2_mtt_cosThetal_R->GetYaxis()->SetTitle("cos#theta_{l}");
    h2_mtt_cosThetal_R->Sumw2();

    h2_mtt_cos1cos2 = new TH2D("mtt_cos1cos2", "m_{tt} cos#theta_{l+}cos#theta_{l-}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cos1cos2->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_cos1cos2->GetYaxis()->SetTitle("cos#theta_{l+}cos#theta_{l-}");
    h2_mtt_cos1cos2->Sumw2();

    h2_HT_deltaPhi = new TH2D("HT_deltaphi", "H_{T} #Delta#phi_{l}", nbins, 0, 4, 10, 0, 1);
    h2_HT_deltaPhi->GetXaxis()->SetTitle("H_{T} [TeV]");
    h2_HT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_HT_deltaPhi->Sumw2();

    h2_KT_deltaPhi = new TH2D("KT_deltaphi", "K_{T} #Delta#phi_{l}", nbins, 0, 4, 10, 0, 1);
    h2_KT_deltaPhi->GetXaxis()->SetTitle("K_{T} [TeV]");
    h2_KT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
    h2_KT_deltaPhi->Sumw2();

    std::vector<string> deltaRnames, deltaRtitles;
    std::vector<string> n_eta, t_eta, n_pt, t_pt;

    std::vector<string> particles1 = {"b1", "b2", "l", "v", "q1", "q2"};
    std::vector<string> particles2 = {"b", "#bar{b}", "l+", "#nu", "q", "#bar{q}'"};

    for (int i = 0; i < 6; i++ ) {
        n_eta.push_back("eta" + particles1[i]);
        t_eta.push_back("#eta_{" + particles2[i] + "}");
    }

    for (int i = 0; i < 6; i++ ) {
        n_pt.push_back("pt" + particles1[i]);
        t_pt.push_back("#p_{T}^{" + particles2[i] + "}");
    }

    for (int i = 0; i < 6; i++ ) {
        for (int j = i + 1; j < 6; j++) {
            deltaRnames.push_back("deltaR" + particles1[i] + particles1[j]);
            deltaRtitles.push_back("#Delta R (" + particles2[i] + ", " + particles2[j] + ")");
        }
    }

    h_deltaRbW = new TH1D("deltaRbW", "#Delta#R(bW)", 100, 0, 5);
    h_deltaRbW->Sumw2();
    h_deltaRmax = new TH1D("deltaRmax", "#Delta#R(max)", 100, 0, 5);
    h_deltaRmax->Sumw2();

    for (int i = 0; i < (int) deltaRnames.size(); i++)
        h_deltaRs.push_back(new TH1D(deltaRnames[i].c_str(), deltaRtitles[i].c_str(), 100, 0, 5));

    for (int i = 0; i < (int) n_eta.size(); i++) {
        h_eta.push_back(new TH1D(n_eta[i].c_str(), t_eta[i].c_str(), 100, 0, 5));
        h_pt.push_back(new TH1D(n_pt[i].c_str(), t_pt[i].c_str(), 400, 0, 100));
    }

    if (m_reco > 0) {
        h_mtt_R = new TH1D("mtt_R", "m_{tt}^{reco}", nbins, Emin, Emax);
        h_mtt_R->Sumw2();
        h_mt_R = new TH1D("mt_R", "m_{t}^{reco}", nbins, 0, 350);
        h_mt_R->Sumw2();
        h_mtbar_R = new TH1D("mtbar_R", "m^{reco}_{#bar{t}}", nbins, 0, 350);
        h_mtbar_R->Sumw2();
        h_mtt_FR = new TH1D("mtt_FR", "m_{tt}^{forward} (reco)", nbins, Emin, Emax);
        h_mtt_FR->Sumw2();
        h_mtt_BR = new TH1D("mtt_BR", "m_{tt}^{backward} (reco)", nbins, Emin, Emax);
        h_mtt_BR->Sumw2();
        h_mtt_FD = new TH1D("mtt_FD", "m_{tt}^{forward} (reco)", nbins, Emin, Emax);
        h_mtt_FD->Sumw2();
        h_mtt_BD  = new TH1D("mtt_BD", "m_{tt}^{backward} (reco)", nbins, Emin, Emax);
        h_mtt_BD->Sumw2();
        h_ytt_R = new TH1D("ytt_R", "y_{tt}^{reco}", 50, -2.5, 2.5);
        h_ytt_R->Sumw2();
        h_cosTheta_R = new TH1D("cosTheta_R", "cos#theta_{reco}", nbins, -1.0, 1.0);
        h_cosTheta_R->Sumw2();
        h_cosThetaStar_R = new TH1D("cosThetaStar_R", "cos#theta_{reco}^{*}", nbins, -1.0, 1.0);
        h_cosThetaStar_R->Sumw2();
        h_pv1x_R = new TH1D("pv1x_r", "p_{x}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
        h_pv1x_R->Sumw2();
        h_pv1y_R = new TH1D("pv1y_r", "p_{y}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
        h_pv1y_R->Sumw2();
        h_pv1z_R = new TH1D("pv1z_r", "p_{z}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
        h_pv1z_R->Sumw2();
        h_pv2x_R = new TH1D("pv2x_r", "p_{x}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
        h_pv2x_R->Sumw2();
        h_pv2y_R = new TH1D("pv2y_r", "p_{y}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
        h_pv2y_R->Sumw2();
        h_pv2z_R = new TH1D("pv2z_r", "p_{z}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
        h_pv2z_R->Sumw2();
    }
}


void Analysis::MakeDistributions()
{
    this->MakeDistribution1D(h_mtt, "TeV");
    this->MakeDistribution1D(h_mtt_F, "TeV");
    this->MakeDistribution1D(h_mtt_B, "TeV");
    this->MakeDistribution1D(h_mtt_Fl, "TeV");
    this->MakeDistribution1D(h_mtt_Bl, "TeV");
    this->MakeDistribution1D(h_mt, "TeV");
    this->MakeDistribution1D(h_mtbar, "TeV");
    this->MakeDistribution1D(h_ytt, "");
    this->MakeDistribution1D(h_cosTheta, "");
    this->MakeDistribution1D(h_cosThetaStar, "");
    this->MakeDistribution1D(h_deltaPhi, "rad / #pi");
    this->MakeDistribution1D(h_pv1x, "GeV");
    this->MakeDistribution1D(h_pv1y, "GeV");
    this->MakeDistribution1D(h_pv1z, "GeV");
    this->MakeDistribution1D(h_pv2x, "GeV");
    this->MakeDistribution1D(h_pv2y, "GeV");
    this->MakeDistribution1D(h_pv2z, "GeV");
    this->MakeDistribution1D(h_cosTheta1, "");
    this->MakeDistribution1D(h_cosTheta2, "");
    this->MakeDistribution1D(h_cos1cos2, "");
    this->MakeDistribution1D(h_HT, "TeV");
    this->MakeDistribution1D(h_KT, "TeV");
    this->MakeDistribution1D(h_deltaRbW, "");
    this->MakeDistribution1D(h_deltaRmax, "");
    this->MakeDistribution1D(h_KT, "TeV");

    h_mtt_Fn = (TH1D*) h_mtt_F->Clone("h_mtt_Fn");
    h_mtt_Fn->Divide(h_mtt);

    h_mtt_Bn = (TH1D*) h_mtt_B->Clone("h_mtt_Bn");
    h_mtt_Bn->Divide(h_mtt);

    h_AFB = this->Asymmetry("AFB", "A_{FB}*", h_mtt_F, h_mtt_B);
    h_AFB->GetYaxis()->SetTitle(h_AFB->GetTitle());
    h_AFB->GetXaxis()->SetTitle("m_{tt} [TeV]");

    h_Ap = this->Asymmetry("Ap", "A_{P}", h_mtt_Fl, h_mtt_Bl);
    h_Ap->GetYaxis()->SetTitle(h_Ap->GetTitle());
    h_Ap->GetXaxis()->SetTitle("m_{tt} [TeV]");

    for (auto& h_deltaR : h_deltaRs) {
        h_deltaR->GetYaxis()->SetTitle("d#sigma / d #Delta R");
        h_deltaR->GetXaxis()->SetTitle("#Delta R");
    }

    for (auto& h : h_eta) {
        h->GetYaxis()->SetTitle("d#sigma / d #eta");
        h->GetXaxis()->SetTitle("#eta");
    }

    for (auto& h : h_pt) {
        h->GetYaxis()->SetTitle("d#sigma / d p_{T}");
        h->GetXaxis()->SetTitle("p_{T}");
    }

    if (m_reco > 0) {
        this->MakeDistribution1D(h_mtt_R, "TeV");
        this->MakeDistribution1D(h_mtt_FR, "TeV");
        this->MakeDistribution1D(h_mtt_BR, "TeV");
        this->MakeDistribution1D(h_mt_R, "TeV");
        this->MakeDistribution1D(h_mtbar_R, "TeV");
        this->MakeDistribution1D(h_ytt_R, "TeV");
        this->MakeDistribution1D(h_cosTheta_R, "");
        this->MakeDistribution1D(h_cosThetaStar_R, "");
        this->MakeDistribution1D(h_pv1x_R, "GeV");
        this->MakeDistribution1D(h_pv1y_R, "GeV");
        this->MakeDistribution1D(h_pv1z_R, "GeV");
        this->MakeDistribution1D(h_pv2x_R, "GeV");
        this->MakeDistribution1D(h_pv2y_R, "GeV");
        this->MakeDistribution1D(h_pv2z_R, "GeV");

        h_mtt_FRn = (TH1D*) h_mtt_FR->Clone("h_mtt_FRn");
        h_mtt_FRn->Divide(h_mtt_R);

        h_mtt_BRn = (TH1D*) h_mtt_BR->Clone("h_mtt_BRn");
        h_mtt_BRn->Divide(h_mtt_R);

        h_AFB_R = this->Asymmetry("AFB_R", "A_{FB}^{reco}", h_mtt_FR, h_mtt_BR);
        h_AFB_R->GetYaxis()->SetTitle(h_AFB_R->GetTitle());
        h_AFB_R->GetXaxis()->SetTitle("m_{tt}^{reco} [TeV]");
    }
}


void Analysis::MakeDistribution1D(TH1D* h, const TString& units)
{
    TString ytitle, yunits, xunits;
    if (m_xsec) {
        ytitle = "d#sigma / d"; //pp->t#bar{t}->b#bar{b}l^{+}l^{-}#nu#bar{#nu}
        if (units != "") {
            yunits = " [fb/" + units + "]";
            xunits = " [" + units + "]";
        }
        else{
            yunits = "";
            xunits = "";
        }
    }
    else {
        ytitle = "Generated events";
        xunits = " [" + units + "]";
        if (units != "") xunits = " [" + units + "]";
        else xunits = "";
    }
    h->GetYaxis()->SetTitle(ytitle + h->GetTitle() + yunits);
    h->GetXaxis()->SetTitle(h->GetTitle() + xunits);
    if (m_xsec && m_useLumi) this->ApplyLuminosity(h);
    else h->Scale(1,"width");
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::MakeDistribution2D(TH2D* h) {
  double sigma, N, dN;
  int k;
  if (m_xsec && m_useLumi) {
    for (int i = 1; i < h->GetNbinsX() + 1; i++) {
      for (int j = 1; j < h->GetNbinsY() + 1; j++) {
        k = h->GetBin(i, j);
        sigma = h->GetBinContent(k);
        N = sigma*m_luminosity*m_efficiency;
        h->SetBinContent(k, N);
        dN = sqrt(N);
        h->SetBinError(k, dN);
      }
    }
    h->GetZaxis()->SetTitle("Expected events");
  }
  else h->Scale(1, "width");
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
    m_outputFile->cd();
    m_outputFile->cd("/");

    h_mtt_Fn->Write();
    h_mtt_Bn->Write();
    h_AFB->Write();
    h_Ap->Write();

    for (int i = 0; i < (int) h_deltaRs.size(); i++)
        h_deltaRs[i]->Write();

    for (int i = 0; i < (int) h_eta.size(); i++)
        h_eta[i]->Write();

    for (int i = 0; i < (int) h_pt.size(); i++)
        h_pt[i]->Write();

    if (m_reco > 0) {
        h_mtt_FRn->Write();
        h_mtt_BRn->Write();
        h_AFB_R->Write();
    }

    this->MakeDistribution2D(h2_mtt_deltaPhi);
    this->MakeDistribution2D(h2_mtt_cos1cos2);
    this->MakeDistribution2D(h2_HT_deltaPhi);
    this->MakeDistribution2D(h2_KT_deltaPhi);
    this->MakeDistribution2D(h2_mtt_cosThetaStar);
    this->MakeDistribution2D(h2_mtt_cosThetaStar_R);
    this->MakeDistribution2D(h2_mtt_cosTheta1);
    this->MakeDistribution2D(h2_mtt_cosTheta2);
    this->MakeDistribution2D(h2_mtt_cosThetal_R);

    TF1 *func = new TF1("func1", "[0]*x + [1]", -1, 1);
    TObjArray slices1;
    this->NormalizeSliceY(h2_mtt_cosTheta1);
    h2_mtt_cosTheta1->FitSlicesY(func, 0, -1, 0, "QRN", &slices1);
    for (auto slice : slices1) slice->Write();
    TH1D* h_AL1 = (TH1D*) slices1[0]->Clone("AL1");
    slices1.Clear();

    TObjArray slices2;
    this->NormalizeSliceY(h2_mtt_cosTheta2);
    h2_mtt_cosTheta2->FitSlicesY(func, 0, -1, 0, "QRN", &slices2);
    for (auto slice : slices2) slice->Write();
    TH1D* h_AL2 = (TH1D*) slices2[0]->Clone("AL2");
    slices2.Clear();

    TObjArray slicesR;
    this->NormalizeSliceY(h2_mtt_cosThetal_R);
    h2_mtt_cosThetal_R->FitSlicesY(func, 0, -1, 0, "QRN", &slicesR);
    for (auto slice : slicesR) slice->Write();
    TH1D* h_AL_R = (TH1D*) slicesR[0]->Clone("AL_R");
    slicesR.Clear();

    // printf("Bin width = %f\n", 2*h2_mtt_cosTheta1->GetYaxis()->GetBinWidth(1));
    h_AL1->Scale(2/h2_mtt_cosTheta1->GetYaxis()->GetBinWidth(1));
    h_AL1->SetTitle("A_{L}");
    h_AL1->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AL1->GetYaxis()->SetTitle("A_{L}");
    h_AL1->Write();
    h_AL2->Scale(2/h2_mtt_cosTheta2->GetYaxis()->GetBinWidth(1));
    h_AL2->SetTitle("A_{L}");
    h_AL2->GetXaxis()->SetTitle("m_{tt} [TeV]");
    h_AL2->GetYaxis()->SetTitle("A_{L}");
    h_AL2->Write();
    h_AL_R->Scale(2/h2_mtt_cosThetal_R->GetYaxis()->GetBinWidth(1));
    h_AL_R->SetTitle("A_{L} (reco)");
    h_AL_R->GetXaxis()->SetTitle("m_{tt} (reco) [TeV]");
    h_AL_R->GetYaxis()->SetTitle("A_{L} (reco)");
    h_AL_R->Write();

    m_outputFile->Close();
    delete m_outputFile;
}


// bool Analysis::PassCuts(string& type) 
// {
//     if(!(type == "truth" or type == "R1" or type == "R2")) return false;
//     if(this->PassCutsET(type)) {
//         if(this->PassCutsEta(type)) {
//             if(this->PassCutsMET(type)) {
//                 if(this->PassCutsMtt(type)) {
//                     if(this->PassCutsYtt(type)) {
//                         return true;
//                     }
//                 }
//             }
//         }
//     }
//     return false;
// }

// bool Analysis::PassCutsMET(string& type)
// {
//     bool pass;
//     pass = true;

//     this->UpdateCutflow(c_MET, pass);
//     return pass;
// }

// bool Analysis::PassCutsMtt(string& type)
// {
//     double mtt;
//     if (type == "truth") mtt = P.M()/1000;
//     else if (type == "R1") mtt = P_R1.M()/1000;
//     else if (type == "R2") mtt = P_R2.M()/1000;
//     else return false;

//     if (m_Emin < 0 and m_Emax < 0) return true;

//     if (mtt > m_Emin) {
//         if (mtt < m_Emax) {
//             UpdateCutflow(c_mtt, true);
//             return true;
//         }
//     }
//     UpdateCutflow(c_mtt, false);
//     return false;
// }

// bool Analysis::PassCutsEta(string& type)
// {
//     if (m_fid == false) {
//         UpdateCutflow(c_eta, true);
//         return true;
//     }
//     if (type == "truth") {
//         for (unsigned int i = 0; i < p.size(); i++) {
//             // bool outsideCrack = abs(p[i].PseudoRapidity()) <= 1.37 || abs(p[i].PseudoRapidity()) >= 1.52;
//             bool central = std::abs(p[i].PseudoRapidity()) <= 2.5;
//             // bool passesEtaCuts = outsideCrack && central;
//             if (central == false) {
//                 UpdateCutflow(c_eta, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else if (type == "R1") {
//         for (unsigned int i = 0; i < p_R1.size(); i++) {
//             // bool outsideCrack = p_R1[i].PseudoRapidity() <= 1.37 || p_R1[i].PseudoRapidity() >= 1.52;
//             bool central = p_R1[i].PseudoRapidity() <= 2.5;
//             // bool passesEtaCuts = outsideCrack && central;
//             if (central == false) {
//                 UpdateCutflow(c_eta, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else if (type == "R2") {
//         for (unsigned int i = 0; i < p_R2.size(); i++) {
//             // bool outsideCrack = p_R2[i].PseudoRapidity() <= 1.37 || p_R2[i].PseudoRapidity() >= 1.52;
//             bool central = p_R2[i].PseudoRapidity() <= 2.5;
//             // bool passesEtaCuts = outsideCrack && central;
//             if (central == false) {
//                 UpdateCutflow(c_eta, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else return false;
//     UpdateCutflow(c_eta, true);
//     return true;
// }

// bool Analysis::PassCutsET(string& type)
// {
//     if (m_fid == false) {
//         UpdateCutflow(c_Et, true);
//         return true;
//     }
//     if (type == "truth") {
//         for (unsigned int i = 0; i < p.size(); i++) {
//             if(p[i].Et() <= 25) {
//                 UpdateCutflow(c_Et, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else if (type == "R1") {
//         for (unsigned int i = 0; i < p_R1.size(); i++) {
//             if(p_R1[i].Et() <= 25) {
//                 UpdateCutflow(c_Et, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else if (type == "R2") {
//         for (unsigned int i = 0; i < p_R2.size(); i++) {
//             if(p_R2[i].Et() <= 25) {
//                 UpdateCutflow(c_Et, false);
//                 return false;
//             }
//             else continue;
//         }
//     }
//     else return false;
//     UpdateCutflow(c_Et, true);
//     return true;
// }

// bool Analysis::PassCutsYtt(string& type)
// {
//     double ytt;
//     if (type == "truth") ytt = std::abs(P.Rapidity());
//     else if (type == "R1") ytt = std::abs(P_R1.Rapidity());
//     else if (type == "R2") ytt = std::abs(P_R2.Rapidity());
//     else return false;
//     if (ytt > m_ytt) {
//         UpdateCutflow(c_ytt, true);
//         return true;
//     }
//     UpdateCutflow(c_ytt, false);
//     return false;
// }


void Analysis::PreLoop()
{
    this->SetDataDirectory();
    this->SetupInputFiles();
    this->SetupOutputFiles();
    this->ResetCounters();
    this->InitialiseCutflow();
    this->MakeHistograms();
}


void Analysis::SetDataDirectory()
{

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    string Hostname(hostname);

    if (Hostname == "Sunder")
        m_dataDirectory = "/Users/declan/Data/zprime";
    else if ((Hostname.find("lxplus") != std::string::npos) || (Hostname.find("cern") != std::string::npos))
        m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime";
    else if ((Hostname.find("cyan") != std::string::npos) || (Hostname.find("blue") != std::string::npos) || (Hostname.find("green") != std::string::npos))
        m_dataDirectory = "/scratch/dam1g09/zprime";
    else
        printf("Hostname %s not recognised.\n", Hostname.c_str());
}


TString Analysis::GetOutputFilename()
{
    return m_outputFilename;
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


void Analysis::GetCrossSection(TString filename)
{
    filename.Replace(filename.Last('.'), 5, ".log");
    printf("%s\n", filename.Data());
    std::ifstream logstream(filename.Data());
    if (!logstream.is_open()) printf("Error: failed to open %s!\n", filename.Data());
    string line;
    string target = " Cross section";
    std::vector<string> parts;
    bool found = false;
    while(getline(logstream, line)) {
        trim(line);
        boost::split(parts, line, boost::is_any_of(":"));
        for (auto& part : parts) trim(part);
        if (parts[0] == target) {
            m_sigma = std::stod(parts[1]);
            found = true;
        }
    }
    logstream.close();
    if (!found) {
        printf("Error: Failed to read generation cross section. Check target log file: %s", filename.Data());
        exit(1);
    }
    else printf("Generation Cross section = %.15le [pb]\n", m_sigma);
}


void Analysis::GetIterationWeights(TString log)
{
    iteration_weights.clear();
    log.Replace(log.Last('.'), 5, ".log");
    std::ifstream logstream(log.Data());
    if (!logstream.is_open()) printf("Error: failed to open %s!\n", log.Data());
    string line;
    string target = " Iteration weighting";
    std::vector<string> parts;
    bool found = false;
    while(getline(logstream, line)) {
        trim(line);
        split(parts, line, boost::is_any_of(":"));
        for (auto& part: parts) trim(part);
        if (parts[0] == target) {
            iteration_weights.push_back(stod(parts[2]));
            found = true;
        }
    }
    logstream.close();
    if (!found) {
        printf("Error: Failed to read Vegas iteration weights. Check %s", log.Data());
        exit(1);
    }
}


void Analysis::Loop()
{
    for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        printf("\n-- Input %li --\n", i - m_inputFiles->begin() + 1);
        cout << (*i) << endl;
        this->SetupTreesForNewFile((*i));
        this->GetCrossSection(*i);
        this->GetIterationWeights(*i);
        this->GetChannelFactors();
        Long64_t nEntries;
        nEntries = this->TotalEvents();
        printf("Processing %lld entries...\n", nEntries);
        for (Long64_t jentry = 0; jentry < nEntries; ++jentry) {
            Long64_t ientry = this->IncrementEvent(jentry);
            if (ientry < 0) break;
            this->EachEvent();
            ProgressBar(jentry, nEntries - 1, 50);
        }
        this->CleanUp();
    }
    cout << endl;
}


Analysis::~Analysis()
{
    delete m_inputFiles;
}


Long64_t Analysis::TotalEvents()
{
    if (m_ntup != 0) {return m_ntup->totalEvents();}
    return -999;
}


Long64_t Analysis::IncrementEvent(Long64_t i)
{
    Long64_t ev(-1);
    if (m_ntup != 0) {ev = m_ntup->LoadTree(i);}
    return ev;
}


void Analysis::SetupTreesForNewFile(const TString& s)
{
    TString treeToUse = "RootTuple";

    m_chainNtup = new TChain(treeToUse,"");
    TString TStringNtuple = s + "/" + treeToUse;
    m_chainNtup->Add(TStringNtuple,0);
    m_ntup = new RootTuple(m_chainNtup);
}


void Analysis::CleanUp()
{
    delete m_chainNtup;
    delete m_ntup;
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
    if (Q_l == 1) {
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
        exit(1);
    }

    double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
    double px_nu = p_nu.Px(), py_nu = p_nu.Py();
    E_l = sqrt(px_l * px_l + py_l * py_l + pz_l * pz_l);
    if (std::abs(E_l - p_l.E()) > 0.00001) printf("ERROR: Lepton energy doesn't match.\n");

    k = 3218.42645 + px_l * px_nu + py_l * py_nu; // m_Wmass * m_Wmass/2 = 3218.42645
    a = px_l * px_l + py_l * py_l;
    b = -2 * k * (pz_l);
    c = (px_nu * px_nu + py_nu * py_nu) * E_l * E_l - k * k;

    double roots[2];
    int nRealRoots = SolveP2(roots, b/a, c/a);
    p_nu_R.clear();
    if (nRealRoots == 2) {
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
    bool oldreco = false;

    if (oldreco) {
        std::vector<TLorentzVector> p_b(2), p_q(2);

        p_b[0] = p[0];
        p_b[1] = p[1];
        if (Q_l == 1) {
            p_q[0] = p[4];
            p_q[1]= p[5];
        }
        else{
            p_q[0] = p[2];
            p_q[1]= p[3];
        }
        imin = 0;
        jmin = 0;
        for (int i = 0; i < (int) p_nu_R.size(); i++) {
            for (int j = 0; j < 2; j++) {
                mblv = (p_b[j] + p_l + p_nu_R[i]).M();
                mjjb = (p_b[1-j] + p_q[0] + p_q[1]).M();
                dh = mjjb - m_tmass;
                dl = mblv - m_tmass;
                chi2 = dh * dh + dl * dl;
                if (chi2 < chi2min) {
                    chi2min = chi2;
                    imin = i;
                    jmin = j;
                }
                // printf("i = %i, j = %i: m_bjj = %.15le, mblv = %.15le, chi2 = %.15le\n", i, j, mjjb, mblv, chi2);
            }
        }
        if (Q_l == 1) {
            p_R[0] = p_b[jmin];
            p_R[1] = p_b[1-jmin];
            p_R[2] = p[2];
            p_R[3] = p_nu_R[imin];
            p_R[4] = p[4];
            p_R[5] = p[5];
        }
        else if (Q_l == -1) {
            p_R[0] = p_b[1-jmin];
            p_R[1] = p_b[jmin];
            p_R[2] = p[2];
            p_R[3] = p[3];
            p_R[4] = p[4];
            p_R[5] = p_nu_R[imin];
        }
    }
    else {
        std::vector<TLorentzVector> p_q(4);

        p_q[0] = p[0];
        p_q[1] = p[1];
        if (Q_l == 1) {
            p_q[2] = p[4];
            p_q[3]= p[5];
        }
        else if (Q_l == -1) {
            p_q[2] = p[2];
            p_q[3]= p[3];
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
                chi2 = dh*dh + dl*dl;
                if (chi2 < chi2min) {
                    chi2min = chi2;
                    imin = i;
                    jmin = j;
                }
            }
        }

        if (Q_l == 1) {
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
    }

    // Assess b-matching performance
    // unsigned int b_lep;
    // if (Q_l == 1) b_lep = 0;
    // if (Q_l == -1) b_lep = 1;
    // if (b_lep == jmin) m_nQuarksMatched++;

    // Assess neutrino reconstruction performance.
    // double pz_nu_truth = p_nu.Pz();
    // double Root0MinusTruth = std::abs(roots[0].real() - pz_nu_truth);
    // double Root1MinusTruth = std::abs(roots[1].real() - pz_nu_truth);
    // unsigned int bestRoot;
    // if (Root0MinusTruth < Root1MinusTruth) bestRoot = 0;
    // else if (Root1MinusTruth < Root0MinusTruth) bestRoot = 1;
    // else bestRoot = 0;
    // if (imin == bestRoot) m_nNeutrinoMatched++;
    // Print reconstruction performance.
    // printf("True pz_nu = %f\n", p_nu.Pz());
    // printf("Possible neutrino solutions:\n");
    // printf("                             %f + %fi\n", roots[0].real(), roots[0].imag());
    // printf("                             %f + %fi\n", roots[1].real(), roots[1].imag());
    // printf("Chosen solution:             %f + %fi\n", roots[imin].real(), roots[imin].imag());
    // if (imin == bestRoot) printf("Neutrino solution: correct. \n");
    // else printf("Neutrino solution: incorrect. \n");
    // if (b_lep == jmin) printf("b-assignment: correct. \n");
    // else printf("b-assignment: incorrect. \n");
    // printf("--\n");
    return p_R;
}


std::vector<TLorentzVector> Analysis::ReconstructDilepton(const std::vector<TLorentzVector>& p)
{
  // Uses the Sonnenschein method algebraically solve tt dilepton equations.
  // http://arxiv.org/abs/hep-ph/0510100
  // selects solution that minimises mtt

  // this->UpdateCutflow(c_events, true);
  if (m_debug) printf("--- start dilepton reconstruction ---\n");

  m_nReco++;

  std::vector<TLorentzVector> p_R(p.size());

  TLorentzVector pb1 = p[0], pb2 = p[1], pl1 = p[2], pv1 = p[3], pl2 = p[4], pv2 = p[5];

  double pb1x = pb1.Px(), pb1y = pb1.Py(), pb1z = pb1.Pz(), Eb1 = pb1.E();
  double pb2x = pb2.Px(), pb2y = pb2.Py(), pb2z = pb2.Pz(), Eb2 = pb2.E();
  double pl1x = pl1.Px(), pl1y = pl1.Py(), pl1z = pl1.Pz(), El1 = pl1.E();
  double pl2x = pl2.Px(), pl2y = pl2.Py(), pl2z = pl2.Pz(), El2 = pl2.E();
  double pv1x = pv1.Px(), pv1y = pv1.Py(), pv1z = pv1.Pz(), Ev1 = pv1.E();
  double pv2x = pv2.Px(), pv2y = pv2.Py(), pv2z = pv2.Pz(), Ev2 = pv2.E();
  double Emissx = pv1x + pv2x;
  double Emissy = pv1y + pv2y;

  // Use on-shell pole masses
  double mt1 = m_tmass, mt2 = m_tmass;
  double mw1 = m_Wmass, mw2 = m_Wmass;
  double mb1 = m_bmass, mb2 = m_bmass;
  double ml1 = 0, ml2 = 0;

  // Use off-shell true masses
  // double mt1 = (p[0] + p[2] + p[3]).M();
  // double mt2 = (p[1] + p[4] + p[5]).M();
  // double mw1 = (p[2] + p[3]).M(), mw2 = (p[4] + p[5]).M();
  // double mb1 = p[0].M();
  // double mb2 = p[1].M();

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

  double c22 = (mw1 * mw1 - ml1 * ml1) * (mw1 * mw1 - ml1 * ml1)
               -4 * (El1 * El1 - pl1z * pl1z) * (a1 / a4) * (a1 / a4)
               -4 * (mw1 * mw1 - ml1 * ml1) * pl1z * a1 / a4;

  double c21 = 4 * (mw1*mw1 - ml1*ml1) * (pl1x - pl1z * a2 / a4)
               -8 * (El1 * El1 - pl1z*pl1z) * a1 * a2 / (a4 * a4) 
               -8 * pl1x * pl1z * a1 / a4;

  double c20 = -4 * (El1 * El1 - pl1x * pl1x)
               -4 * (El1 * El1 - pl1z * pl1z) * (a2 / a4) * (a2 / a4)
               -8 * pl1x * pl1z * a2 / a4;

  double c11 = 4 * (mw1*mw1 - ml1*ml1)*(pl1y - pl1z * a3 / a4)
               -8 * (El1 * El1 - pl1z * pl1z) * a1 * a3 / (a4 * a4) 
               -8 * pl1y * pl1z * a1 /  a4;

  double c10 = -8 * (El1*El1 - pl1z*pl1z) * a2 * a3 / (a4 * a4) 
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

  double dd22 = (mw2 * mw2 - ml2 * ml2) * (mw2 * mw2 - ml2 * ml2)
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
               - 2 * Emissx * dd20 
               - Emissy * dd10;

  double d20 = dd20;

  double d11 = -dd11 
               - 2 * Emissy * dd00 
               - Emissx * dd10;

  double d10 = dd10;
  double d00 = dd00;

  const double h4 = c00 * c00 * d22 * d22 
                  + c11 * d22 * (c11 * d00 - c00 * d11)
                  + c00 * c22 * (d11 * d11 - 2 * d00 * d22) 
                  + c22 * d00 * (c22 * d00 - c11 * d11);

  const double h3 = c00 * d21 * (2 * c00 * d22 - c11 * d11) 
                  + c00 * d11 * (2 * c22 * d10 + c21 * d11) 
                  + c22 * d00 * (2 * c21 * d00 - c11 * d10) 
                  - c00 * d22 * (c11 * d10 + c10 * d11)  
                  -2 * c00 * d00 * (c22 * d21 + c21 * d22)
                  - d00 * d11 * (c11 * c21 + c10 * c22) 
                  + c11 * d00 * (c11 * d21 + 2 * c10 * d22);

  const double h2 = c00 * c00 * (2 * d22 * d20 + d21 * d21) 
                  - c00 * d21 * (c11 * d10 + c10 * d11)  
                  + c11 * d20 * (c11 * d00 - c00 * d11) 
                  + c00 * d10 * (c22 * d10 - c10 * d22)   
                  + c00 * d11 * (2 * c21 * d10 + c20 * d11) 
                  + (2 * c22 * c20 + c21 * c21) * d00 * d00   
                  - 2 * c00 * d00 * (c22 * d20 + c21 * d21 + c20 * d22)    
                  + c10 * d00 * (2 * c11 * d21 + c10 * d22) 
                  - d00 * d10 * (c11 * c21 + c10 * c22)   
                  - d00 * d11 * (c11 * c20 + c10 * c21);

  const double h1 = c00 * d21 * (2 * c00 * d20 - c10 * d10) 
                  - c00 * d20 * (c11 * d10 + c10 * d11)  
                  + c00 * d10 * (c21 * d10 + 2 * c20 * d11) 
                  - 2 * c00 * d00 * (c21 * d20 + c20 * d21)  
                  + c10 * d00 * (2 * c11 * d20 + c10 * d21) 
                  + c20 * d00 * (2 * c21 * d00 - c10 * d11)  
                  - d00 * d10 * (c11 * c20 + c10 * c21);

  const double h0 = c00 * c00 * d20 * d20 
                  + c10 * d20 * (c10 * d00 - c00 * d10)  
                  + c20 * d10 * (c00 * d10 - c10 * d00) 
                  + c20 * d00 * (c20 * d00 - 2 * c00 * d20);

  int dig = DECIMAL_DIG;
  if (m_debug) {
    printf("h4 = %.*e\n", dig, h4);
    printf("h3 = %.*e\n", dig, h3);
    printf("h2 = %.*e\n", dig, h2);   
    printf("h1 = %.*e\n", dig, h1);
    printf("h0 = %.*e\n", dig, h0);
    printf("\n");
  }

  double a[5] = {1.0, h1/h0, h2/h0, h3/h0, h4/h0};

  if (m_debug) for (int i = 0; i < 5; i++) cout << "a(" << i << ") = "<<  a[i] << endl; 

  double x[4];
  const int nRealRoots = SolveP4(x, a[1], a[2], a[3], a[4]);
  // If nRealRoots = 4, they live in x[0], x[1], x[2], x[3].
  // If nRealRoots = 2, x[0], x[1] are the real roots and x[2]±i*x[3] are the complex.
  // If nRealRoots = 0, the equation has two pairs of pairs of complex conjugate roots in x[0]±i*x[1] and x[2]±i*x[3].

  if (m_debug) cout << "Found " << nRealRoots << " real roots" << endl;

  int nSolutions;
  std::vector<double> pv1x_Rs;
  if (nRealRoots == 4) {
    nSolutions = 4;
    for (int i = 0; i < nSolutions; i++) pv1x_Rs.push_back(x[i]);
  }
  else if (nRealRoots == 2) {
    nSolutions = 3;
    for (int i = 0; i < nSolutions; i++) pv1x_Rs.push_back(x[i]); 
  }
  else if (nRealRoots == 0) {
    nSolutions = 2;
    pv1x_Rs.push_back(x[0]);
    pv1x_Rs.push_back(x[2]); 
  }

  if (m_debug) {
    cout << "pv1x = " << pv1x << endl; 
    for (int i = 0; i < 4; i++) cout << "x(" << i << ") = " << x[i] << endl;
  }

  // find root closest to true pl1x
  // double old_diff = std::abs(pv1x - x[0]), new_diff, closest_root = x[0];
  // for (unsigned int i = 1; i < 4; i++) {
  //   new_diff = std::abs(pv1x - x[i]);
  //   if (new_diff < old_diff) {
  //     closest_root = x[i];
  //     old_diff = new_diff;
  //   }
  // }
  // std::vector<double> realRoots;
  // for (int i = 0; i < nRealRoots; i++) realRoots.push_back(x[i]);

  // Create pairs of neutrino momenta for each real root
  std::vector<TLorentzVector> pv1_Rs(nSolutions), pv2_Rs(nSolutions);
  for (int i = 0; i < nSolutions; i++) {
    double pv1x_R = x[i];
    double pv2x_R = Emissx - pv1x_R;

    double c2 = c22 + c21 * pv1x_R + c20 * pv1x_R * pv1x_R;
    double c1 = c11 + c10 * pv1x_R;
    double c0 = c00;
    double d2 = d22 + d21 * pv1x_R + d20 * pv1x_R * pv1x_R;
    double d1 = d11 + d10 * pv1x_R;
    double d0 = d00;

    double pv1y_R = (c0 * d2 - c2 * d0) / (c1 * d0 - c0 * d1);
    double pv2y_R = Emissy - pv1y_R;

    double pv1z_R = -(a1 + a2 * pv1x_R + a3 * pv1y_R) / a4;
    double pv2z_R = -(b1 + b2 * pv2x_R + b3 * pv2y_R) / b4;

    double Ev1_R = sqrt(pv1x_R * pv1x_R + pv1y_R * pv1y_R + pv1z_R * pv1z_R);
    double Ev2_R = sqrt(pv2x_R * pv2x_R + pv2y_R * pv2y_R + pv2z_R * pv2z_R);

    pv1_Rs[i].SetPxPyPzE(pv1x_R, pv1y_R, pv1z_R, Ev1_R);
    pv2_Rs[i].SetPxPyPzE(pv2x_R, pv2y_R, pv2z_R, Ev2_R);

    if (m_debug) {
      printf("pv1x   = %.*e\n", dig, pv1x);
      printf("pv1x_R = %.*e\n", dig, pv1x_R);
      printf("pv2x   = %.*e\n", dig, pv2x);
      printf("pv2x_R = %.*e\n", dig, pv2x_R);
      printf("pv1y   = %.*e\n", dig, pv1y);
      printf("pv1y_R = %.*e\n", dig, pv1y_R);
      printf("pv2y   = %.*e\n", dig, pv2y);
      printf("pv2y_R = %.*e\n", dig, pv2y_R);
      printf("pv1z   = %.*e\n", dig, pv1z);
      printf("pv1z_R = %.*e\n", dig, pv1z_R);
      printf("pv2z   = %.*e\n", dig, pv2z);
      printf("pv2z_R = %.*e\n", dig, pv2z_R);
      printf("Ev1    = %.*e\n", dig, Ev1);
      printf("Ev1_R  = %.*e\n", dig, Ev1_R);
      printf("Ev2    = %.*e\n", dig, Ev2);
      printf("Ev2_R  = %.*e\n", dig, Ev2_R);

    }
  }

  int I = 0;
  double mtt_min = DBL_MAX;
  for (int i = 0; i < nSolutions; i++) {
    double mtt = (p[0] + p[1] + p[2] + p[4] + pv1_Rs[i] + pv2_Rs[i]).M();
    if (mtt < mtt_min) {
      mtt_min = mtt;
      I = i;
    }
  }

  p_R[0] = pb1;
  p_R[1] = pb2;
  p_R[2] = pl1;
  p_R[3] = pv1_Rs[I];
  p_R[4] = pl2;
  p_R[5] = pv2_Rs[I];

  if (m_debug) printf("--- end dilepton reconstruction ---\n\n");
  return p_R;
}


void Analysis::GetChannelFactors()
{
    // scale dilepton to other classifications
    // fac_ee = 1
    // fac_emu = 2
    // if (m_reco) m_sigma = m_sigma*24; // 2 [e+ + e-] x 2 [e + mu] x 6 [3 x 2]
    // fac_qq = 36
}


void Analysis::UpdateCutflow(const int cut, const bool passed)
{
    if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
    if (passed) m_cutflow[cut] += 1;
}


void Analysis::InitialiseCutflow()
{
    m_cutflow = std::vector<int>(m_cuts, -999);
    m_cutNames = std::vector<TString>(m_cuts, "no name");
    m_cutNames[c_entries]       = "Entries         ";
    m_cutNames[c_topDecays]     = "t->be#nu        ";
    m_cutNames[c_antitopDecays] = "#bar{t}->be#nu  ";
    m_cutNames[c_events]        = "Events          ";
    m_cutNames[c_realSolutions] = "p^{#nu}_{z} real";
    m_cutNames[c_MET]           = "MET             ";
    m_cutNames[c_mtt]           = "mtt             ";
    m_cutNames[c_ytt]           = "ytt             ";
    m_cutNames[c_eta]           = "#eta            ";
    m_cutNames[c_Et]            = "E_{T}           ";

    h_cutflow = new TH1D("Cutflow", "Cutflow", m_cuts, 0, m_cuts);
}


void Analysis::PrintCutflow()
{
    printf("-- Cutflow --\n");
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut]);

        printf("%s %i pass\n", m_cutNames[cut].Data(), m_cutflow[cut]);
    }
    h_cutflow->Write();
}


void Analysis::SetYttCut(const double ytt)
{
    m_ytt = ytt;
}


void Analysis::SetXsec(const bool xsec)
{
    m_xsec = xsec;
}


void Analysis::SetFiducial(const bool fid)
{
    m_fid = fid;
}


void Analysis::SetEnergyRange(double Emin = -1, double Emax = -1)
{
    m_Emin = Emin;
    m_Emax = Emax;
}
