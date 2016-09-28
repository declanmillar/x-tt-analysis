#include "analysis.hpp"

using namespace std;
using namespace boost;

string trim(string const& str) {
    if(str.empty())
    return str;

    size_t firstScan = str.find_first_not_of(' ');
    size_t first = firstScan == string::npos ? str.length() : firstScan;
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, last-first+1);
}

AnalysisZprime::AnalysisZprime(const TString model, const TString initial_state, const TString intermediates, const TString final_state,  const int energy, const TString options, const int vegasIterations, const string vegasPoints, const bool add_ggG, const bool add_qqG, const int luminosity, const int btags, const bool discardComplex, const TString analysisLabel):
    m_model(model),
    m_initial_state(initial_state),
    m_intermediates(intermediates),
    m_channel(final_state),
    m_energy(energy),
    m_options(options),
    m_vegasIterations(vegasIterations),
    m_vegasPoints(vegasPoints),
    m_add_ggG(add_ggG),
    m_add_qqG(add_qqG),
    m_luminosity(luminosity),
    m_efficiency(1.0),
    nBtags(btags),
    m_discardComplex(discardComplex),
    m_analysisLabel(analysisLabel),
    m_reco(1),
    m_pi(3.14159265),
    m_GeV(1000.0),
    m_bmass(4.18),
    m_Wmass(80.23),
    m_tmass(173.0),
    m_ytt(0),
    m_Emin(-1),
    m_Emax(-1),
    m_discardEvent(false),
    m_inputFiles(NULL),
    m_weightFiles(NULL),
    m_ntup(NULL),
    m_chainNtup(NULL),
    m_outputFile(NULL) {
}

void AnalysisZprime::Run(){
    this->PreLoop();
    this->Loop();
    this->PostLoop();
}

inline string BoolToString(bool b) {return b ? "1" : "0";}

TString AnalysisZprime::GetOutputFilename() {
    return m_outputFilename;
}

void AnalysisZprime::EachEvent() {
    m_discardEvent = false;
    UpdateCutflow(c_entries, true);
    p = vector<TLorentzVector>(6);
    for (int i = 0; i < 6; i++)
        p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));

    P.SetPxPyPzE(0, 0, 0, 0);
    pcm = vector<TLorentzVector>(p.size());
    ptop = vector<TLorentzVector>(p.size());
    patop = vector<TLorentzVector>(p.size());
    for (int i = 0; i < 6; i++) {
        P += p[i];
        pcm[i] = p[i];
        ptop[i] = p[i];
        patop[i] = p[i];
    }

    TVector3 V = -1*P.BoostVector();
    TVector3 Vtop = -1*(p[0] + p[2] + p[3]).BoostVector();
    TVector3 Vatop = -1*(p[1] + p[4] + p[5]).BoostVector();

    Pcm.SetPxPyPzE(0, 0, 0, 0);
    Ptop.SetPxPyPzE(0, 0, 0, 0);
    Patop.SetPxPyPzE(0, 0, 0, 0);
    for (int i = 0; i < 6; i++) {
        pcm[i].Boost(V);
        ptop[i].Boost(Vtop);
        patop[i].Boost(Vatop);
        Pcm += pcm[i];
        Ptop += ptop[i];
        Patop += patop[i];
    }


    if (m_reco == 1) {
        p_R1 = this->ReconstructDilepton(p); // both decay leptonically

        P_R1.SetPxPyPzE(0, 0, 0, 0);
        for (int i = 0; i < 6; i++) P_R1 += p_R1[i];

        TVector3 V_R1 = -1*P_R1.BoostVector();
        TVector3 Vtop_R1 = -1*(p_R1[0] + p_R1[2] + p_R1[3]).BoostVector();
        TVector3 Vatop_R1 = -1*(p_R2[1] + p_R2[4] + p_R2[5]).BoostVector();
        pcm_R1 = p_R1;
        ptop_R1 = p_R1;
        patop_R2 = p_R1;
        for (unsigned int i = 0; i < p.size(); i++) {
            pcm_R1[i].Boost(V_R1);
            ptop_R1[i].Boost(Vtop_R1);
            patop_R2[i].Boost(Vatop_R1);
        }
    }
    if (m_reco == 2) {
        p_R1 = this->ReconstructSemiLeptonic(p, 1); // top decays leptonically
        p_R2 = this->ReconstructSemiLeptonic(p, -1); // top decays hadronically

        if (m_discardEvent) return;

        P_R1.SetPxPyPzE(0, 0, 0, 0);
        P_R2.SetPxPyPzE(0, 0, 0, 0);
        for (int i = 0; i < 6; i++) {
            P_R1 += p_R1[i];
            P_R2 += p_R2[i];
        }

        // reconstructed final particle parton CoM variables
        TVector3 V_R1 = -1*P_R1.BoostVector();
        TVector3 V_R2 = -1*P_R2.BoostVector();
        TVector3 Vtop_R1 = -1*(p_R1[0] + p_R1[2] + p_R1[3]).BoostVector();
        TVector3 Vatop_R2 = -1*(p_R2[1] + p_R2[4] + p_R2[5]).BoostVector();
        pcm_R1 = vector<TLorentzVector>(p.size());
        pcm_R2 = vector<TLorentzVector>(p.size());
        ptop_R1 = vector<TLorentzVector>(p.size());
        patop_R2 = vector<TLorentzVector>(p.size());
        for (unsigned int i = 0; i < p.size(); i++) {
            pcm_R1[i] = p_R1[i];
            ptop_R1[i] = p_R1[i];
            pcm_R2[i] = p_R2[i];
            patop_R2[i] = p_R2[i];
            pcm_R1[i].Boost(V_R1);
            pcm_R2[i].Boost(V_R2);
            ptop_R1[i].Boost(Vtop_R1);
            patop_R2[i].Boost(Vatop_R2);
        }
    }

    // top and antitop
    TLorentzVector p_t = p[0] + p[2] + p[3];
    TLorentzVector p_tb = p[1] + p[4] + p[5];
    TLorentzVector p_W = p[2] + p[3];
    pcm_t = pcm[0] + pcm[2] + pcm[3];
    pcm_tb = pcm[1] + pcm[4] + pcm[5];

    if (m_reco > 0) {
        p_t_R1 = pcm_R1[0] + pcm_R1[2] + pcm_R1[3];
        p_tb_R1 = pcm_R1[1] + pcm_R1[4] + pcm_R1[5];
    }
    if (m_reco == 2) {
        p_t_R2 = pcm_R2[0] + pcm_R2[2] + pcm_R2[3];
        p_tb_R2 = pcm_R2[1] + pcm_R2[4] + pcm_R2[5];
    }


    double mtt = P.M()/1000;
    double mtt_R1 = P_R1.M()/1000;
    double mtt_R2 = P_R2.M()/1000;
    double HT = 0;
    for (int i = 0; i < 5; i++) HT = HT + p[i].Pt();
    HT = HT/1000;
    double mvis = (p[0] + p[1] +p[2] + p[4]).M();
    double pTvis = (p[0] + p[1] +p[2] + p[4]).Pt();
    double mTvis = sqrt(mvis*mvis + pTvis*pTvis);
    double KT = mTvis + (p[3] + p[5]).Pt();
    KT = KT/1000;
    double mt = pcm_t.M();
    double mtb = pcm_tb.M();
    double ytt = P.Rapidity();
    double cosTheta = pcm_t.CosTheta();
    double cosThetaStar = int(ytt/abs(ytt))*cosTheta;
    double ytt_R1 = -999;
    double ytt_R2 = -999;
    double mt_R1 = -999;
    double mtb_R1 = -999;
    double mt_R2 = -999;
    double mtb_R2 = -999;
    double cosTheta_R1 = -999;
    double cosTheta_R2 = -999;
    double cosThetaStar_R1 = -999;
    double cosThetaStar_R2 = -999;
    double deltaPhi = -999;
    double cosTheta1 = -999;
    double cosTheta2 = -999;
    double cosTheta1_R1 = -999;
    double cosTheta2_R2 = -999;
    double cos1cos2 = -999;

    vector<double> deltaRs;
    for (int i = 0; i < 6; i++)
        for (int j = i + 1; j < 6; j++)
            deltaRs.push_back(p[i].DeltaR(p[j]));

    vector<double> eta, pt;
    for (int i = 0; i < 6; i++) {
        pt.push_back(p[i].Pt());
        eta.push_back(p[i].Eta());
    }

    double deltaRbW = p_W.DeltaR(p[0]);
    vector<double> deltaRt;
    deltaRt.push_back(p_t.DeltaR(p[0]));
    deltaRt.push_back(p_t.DeltaR(p[2]));
    deltaRt.push_back(p_t.DeltaR(p[3]));
    auto deltaRmax = max_element(std::begin(deltaRt), std::end(deltaRt));

    cosTheta1 = cos(ptop[2].Angle(pcm_t.Vect()));
    cosTheta2 = cos(patop[4].Angle(pcm_tb.Vect()));
    cos1cos2 = cosTheta1*cosTheta2;
    deltaPhi = p[2].DeltaPhi(p[4])/m_pi;

    if(m_reco) {
        ytt_R1 = P_R1.Rapidity();
        ytt_R2 = P_R2.Rapidity();
        mt_R1 = p_t_R1.M();
        mtb_R1 = p_tb_R1.M();
        mt_R2 = p_t_R2.M();
        mtb_R2 = p_tb_R2.M();
        cosTheta_R1 = p_t_R1.CosTheta();
        cosTheta_R2 = p_t_R2.CosTheta();
        cosThetaStar_R1 = int(ytt_R1/abs(ytt_R1))*cosTheta_R1;
        cosThetaStar_R2 = int(ytt_R2/abs(ytt_R2))*cosTheta_R2;
        cosTheta1_R1 = cos(ptop_R1[2].Angle(p_t_R1.Vect()));
        cosTheta2_R2 = cos(patop_R2[4].Angle(p_tb_R2.Vect()));
    }

    // printf("Event passed all cuts.\n");
    // re-weight for different iterations
    double weight;
    double weight_R;
    if (m_xsec) {
        double it = m_ntup->iteration();
        weight = m_ntup->weight();
        double fb = 1000;
        weight = fb*weight*m_sigma/iteration_weights[it-1];
        // printf("Sigma = %.15le\n", m_sigma);
        // printf("Iteration weight = %.15le\n", iteration_weights[it-1]);
        weight_R = weight/2;
    }
    else {
        weight = 1;
        weight_R = 1;
    }

    if (this->PassCuts("truth")) {
        h_mt->Fill(mt, weight);
        h_mtbar->Fill(mtb, weight);
        h_mtt->Fill(mtt, weight);
        h_ytt->Fill(ytt, weight);
        h_cosTheta->Fill(cosTheta, weight);
        h_cosThetaStar->Fill(cosThetaStar, weight);
        h_HT->Fill(HT, weight);
        h_KT->Fill(KT, weight);
        h_pzNu->Fill(p[3].Pz(), weight_R);
        h_deltaPhi->Fill(deltaPhi, weight);
        h_cosTheta1->Fill(cosTheta1, weight);
        h_cosTheta2->Fill(cosTheta2, weight);
        h_cos1cos2->Fill(cos1cos2, weight);
        h_deltaRbW->Fill(deltaRbW, weight);
        h_deltaRmax->Fill(*deltaRmax, weight);

        for (int i = 0; i < (int) deltaRs.size(); i++)
            h_deltaRs[i]->Fill(deltaRs[i], weight);

        for (int i = 0; i < (int) eta.size(); i++) {
            h_eta[i]->Fill(eta[i], weight);
            h_pt[i]->Fill(pt[i], weight);
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
    }

    if (m_reco > 0) {
        if (this->PassCuts("R1")) {
            h_mtt_R->Fill(mtt_R1, weight_R);
            h_ytt_R->Fill(ytt_R1, weight_R);
            h_mt_R->Fill(mt_R1, weight_R);
            h_mtbar_R->Fill(mtb_R1, weight_R);
            h_cosTheta_R->Fill(cosTheta_R1, weight_R);
            h_cosThetaStar_R->Fill(cosThetaStar_R1, weight_R);
            h_pzNu_R->Fill(p_R1[3].Pz(), weight_R);
            if (cosThetaStar_R1 > 0) h_mtt_FR->Fill(mtt_R1, weight_R);
            if (cosThetaStar_R1 < 0) h_mtt_BR->Fill(mtt_R1, weight_R);
            h2_mtt_cosThetaStar_R->Fill(mtt_R1, cosThetaStar_R1, weight_R);
            h2_mtt_cosThetal_R->Fill(mtt_R1, cosTheta1_R1, weight_R);
        }
    }
    if (m_reco == 2) {
        if (this->PassCuts("R2")) {
            h_mtt_R->Fill(mtt_R2, weight_R);
            h_ytt_R->Fill(ytt_R2, weight_R);
            h_mt_R->Fill(mt_R2, weight_R);
            h_mtbar_R->Fill(mtb_R2, weight_R);
            h_cosTheta_R->Fill(cosTheta_R2, weight_R);
            h_cosThetaStar_R->Fill(cosThetaStar_R2, weight_R);
            h_pzNu_R->Fill(p_R2[5].Pz(), weight_R);
            if (cosThetaStar_R2 > 0) h_mtt_FR->Fill(mtt_R2, weight_R);
            if (cosThetaStar_R2 < 0) h_mtt_BR->Fill(mtt_R2, weight_R);
            h2_mtt_cosThetaStar_R->Fill(mtt_R2, cosThetaStar_R2, weight_R);
            h2_mtt_cosThetal_R->Fill(mtt_R2, cosTheta2_R2, weight_R);
        }
    }
}

void AnalysisZprime::SetupInputFiles() {
    m_inputFiles = new vector<TString>;
    m_weightFiles = new vector<TString>;
    TString filename;

    string E = "";
    if (m_energy != 13) "_" + to_string(m_energy);

    if (m_add_ggG) {
      filename = m_dataDirectory + "/SM_" + "gg-G-" + m_channel + E + "_2-4_" + to_string(m_vegasIterations) + "x" + m_vegasPoints;
      m_inputFiles->push_back(filename + ".root");
      m_weightFiles->push_back(filename + ".log");
    }

    if (m_add_qqG) {
      filename = m_dataDirectory + "/SM_" + "qq-" + m_channel + E + "_2-3_" + to_string(m_vegasIterations) + "x" + m_vegasPoints;
      m_inputFiles->push_back(filename + ".root");
      m_weightFiles->push_back(filename + ".log");
      filename = m_dataDirectory + "/SM_" + "qq-" + m_channel + E + "_3-4_" + to_string(m_vegasIterations) + "x" + m_vegasPoints;
      m_inputFiles->push_back(filename + ".root");
      m_weightFiles->push_back(filename + ".log");
    }

    filename = m_dataDirectory + "/" + m_model + "_" + m_initial_state + "-" + m_intermediates + m_channel + E + m_options + to_string(m_vegasIterations) + "x" + m_vegasPoints;
    m_inputFiles->push_back(filename + ".root");
    m_weightFiles->push_back(filename + ".log");

    // filename = m_dataDirectory + "/" + "SM_qq-tt-bbllvv_2-3_5x10M";
    // m_inputFiles->push_back(filename + ".root");
    // m_weightFiles->push_back(filename + ".log");
    // filename = m_dataDirectory + "/" + "SM_qq-tt-bbllvv_3-4_5x10M";
    // m_inputFiles->push_back(filename + ".root");
    // m_weightFiles->push_back(filename + ".log");

    // Check all input files exist
    for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        bool exists;
        struct stat buffer;
        exists = stat ((*i).Data(), &buffer) == 0;
        if (exists == false) {
            printf("Error: %s does not exist.\n", (*i).Data());
            exit(exists);
        }
    }
}

void AnalysisZprime::SetupOutputFiles() {
    TString outfilename;
    TString initial_state = m_initial_state;
    TString intermediates = m_intermediates;
    string E = "";
    if (m_energy != 13) "_" + to_string(m_energy);

    if (m_add_ggG || m_add_qqG) intermediates = "G" + intermediates;
    if (m_add_ggG) initial_state = "gg" + initial_state;
    outfilename = m_dataDirectory + "/" + m_model + "_" + initial_state + "-" + intermediates + m_channel + E + m_options + to_string(m_vegasIterations) + "x" + m_vegasPoints;
    m_outputFilename = outfilename + ".a";

    if (m_reco == 2 && nBtags != 2) m_outputFilename += ".b" + to_string(nBtags) + ".c" + BoolToString(m_discardComplex);
    string ytt = to_string(m_ytt);
    if (m_ytt > 0) m_outputFilename += ".y" + ytt.erase(ytt.find_last_not_of('0') + 1, string::npos);
    if (m_Emin >= 0 || m_Emax >= 0) m_outputFilename += ".E" + to_string(m_Emin) + "-" + to_string(m_Emax);
    string eff = to_string(m_efficiency);
    if (m_efficiency < 1.0) m_outputFilename += ".e" + eff.erase(eff.find_last_not_of('0') + 1, string::npos);
    if (m_luminosity > 0) m_outputFilename += ".L" + to_string(m_luminosity);
    if (m_fid == true) m_outputFilename += ".fid";
    m_outputFilename += m_analysisLabel;
    m_outputFilename += ".root";

    printf("--- Output ---\n");
    printf("%s\n", m_outputFilename.Data());
    m_outputFile = new TFile(m_outputFilename, "RECREATE");
}

void AnalysisZprime::PostLoop () {
    this->CheckResults();
    if (m_reco == 2) this->CheckPerformance();
    this->MakeGraphs();
    this->PrintCutflow();
    this->WriteHistograms();
}

void AnalysisZprime::CheckResults() {
    printf("--- Results ---\n");
    double sigma = h_mtt->Integral("width");
    printf("Analysis Cross Section = %.15le [pb]\n", sigma);
    printf("\n");
}

void AnalysisZprime::CheckPerformance () {
    double quarkRecoRatio = m_nQuarksMatched/(double)m_nReco;
    double neutrinoRecoRatio = m_nNeutrinoMatched/(double)m_nReco;
    printf("--- Performance ---\n");
    printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio*100);
    printf("pzNu assignment: %.1f%% correct\n", neutrinoRecoRatio*100);
    printf("\n");
}

double AnalysisZprime::TotalAsymmetry(TH1D* h1, TH1D* h2) {
    double N1 = h1->Integral("width");
    double N2 = h2->Integral("width");
    return (N1 - N2)/(N1 + N2);
}

void AnalysisZprime::ApplyLuminosity(TH1D* h) {
    double sigma, N, dN;
    for (int i = 1; i < h->GetNbinsX() + 1; i++) {
        sigma = h->GetBinContent(i);
        N = m_luminosity*m_efficiency*sigma;
        h->SetBinContent(i, N);
        dN = sqrt(N);
        h->SetBinError(i, dN);
    }
    h->GetYaxis()->SetTitle("Expected events");
}

TH1D* AnalysisZprime::Asymmetry(TString name, TString title, TH1D* h1, TH1D* h2) {
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

TH1D* AnalysisZprime::Asymmetry2(TString name, TString title, TH1D* h1, TH1D* h2) {
    TH1D* h = (TH1D*) h1->Clone(name);
    h->SetTitle(title);
    h->Add(h2, -1);
    return h;
}

void AnalysisZprime::AsymmetryUncertainty(TH1D* hA, TH1D* h1, TH1D* h2) {
    double A, dA, N, N1, N2;
    for (int i = 1; i < hA->GetNbinsX() + 1; i++) {
        A = hA->GetBinContent(i);
        N1 = h1->GetBinContent(i);
        N2 = h2->GetBinContent(i);
        N = N1 + N2;
        if (N > 0) dA = sqrt((1.0 - A*A)/N);
        else dA = 0;
        hA->SetBinError(i, dA);
    }
}

void AnalysisZprime::CreateHistograms() {

    double binWidth = 0.05;
    double Emin = 2.025;
    double Emax = 3.975;
    double nbins = (Emax - Emin)/binWidth;

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
    h_pzNu = new TH1D("pzNu", "p_{z}^{#nu}", nbins, -500.0, 500.0);
    h_pzNu->Sumw2();
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

    vector<string> deltaRnames, deltaRtitles;
    vector<string> n_eta, t_eta, n_pt, t_pt;

    vector<string> particles1 = {"b1", "b2", "l", "v", "q1", "q2"};
    vector<string> particles2 = {"b", "#bar{b}", "l+", "#nu", "q", "#bar{q}'"};

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
        h_pzNu_R = new TH1D("pzNu_R", "p_{z}^{#nu} (reco)", nbins, -500.0, 500.0);
        h_pzNu_R->Sumw2();
    }
}

void AnalysisZprime::MakeGraphs() {
    this->MakeDistribution(h_mtt, "TeV");
    this->MakeDistribution(h_mtt_F, "TeV");
    this->MakeDistribution(h_mtt_B, "TeV");
    this->MakeDistribution(h_mtt_Fl, "TeV");
    this->MakeDistribution(h_mtt_Bl, "TeV");
    this->MakeDistribution(h_mt, "TeV");
    this->MakeDistribution(h_mtbar, "TeV");
    this->MakeDistribution(h_ytt, "");
    this->MakeDistribution(h_cosTheta, "");
    this->MakeDistribution(h_cosThetaStar, "");
    this->MakeDistribution(h_deltaPhi, "rad / #pi");
    this->MakeDistribution(h_pzNu, "GeV");
    this->MakeDistribution(h_cosTheta1, "");
    this->MakeDistribution(h_cosTheta2, "");
    this->MakeDistribution(h_cos1cos2, "");
    this->MakeDistribution(h_HT, "TeV");
    this->MakeDistribution(h_KT, "TeV");
    this->MakeDistribution(h_deltaRbW, "");
    this->MakeDistribution(h_deltaRmax, "");
    this->MakeDistribution(h_KT, "TeV");

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

    for (auto h_deltaR : h_deltaRs) {
        h_deltaR->GetYaxis()->SetTitle("d#sigma / d #Delta R");
        h_deltaR->GetXaxis()->SetTitle("#Delta R");
    }

    for (auto h : h_eta) {
        h->GetYaxis()->SetTitle("d#sigma / d #eta");
        h->GetXaxis()->SetTitle("#eta");
    }

    for (auto h : h_pt) {
        h->GetYaxis()->SetTitle("d#sigma / d p_{T}");
        h->GetXaxis()->SetTitle("p_{T}");
    }

    if (m_reco > 0) {
        this->MakeDistribution(h_mtt_R, "TeV");
        this->MakeDistribution(h_mtt_FR, "TeV");
        this->MakeDistribution(h_mtt_BR, "TeV");
        this->MakeDistribution(h_mt_R, "TeV");
        this->MakeDistribution(h_mtbar_R, "TeV");
        this->MakeDistribution(h_ytt_R, "TeV");
        this->MakeDistribution(h_cosTheta_R, "");
        this->MakeDistribution(h_cosThetaStar_R, "");
        this->MakeDistribution(h_pzNu_R, "GeV");

        h_mtt_FRn = (TH1D*) h_mtt_FR->Clone("h_mtt_FRn");
        h_mtt_FRn->Divide(h_mtt_R);

        h_mtt_BRn = (TH1D*) h_mtt_BR->Clone("h_mtt_BRn");
        h_mtt_BRn->Divide(h_mtt_R);

        h_AFB_R = this->Asymmetry("AFB_R", "A_{FB}^{reco}", h_mtt_FR, h_mtt_BR);
        h_AFB_R->GetYaxis()->SetTitle(h_AFB_R->GetTitle());
        h_AFB_R->GetXaxis()->SetTitle("m_{tt}^{reco} [TeV]");
    }
}

void AnalysisZprime::MakeDistribution(TH1D* h, TString units) {
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

void AnalysisZprime::Make2dDistribution(TH2D* h) {
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

void AnalysisZprime::NormalizeSliceY(TH2D* h) {
    double integral = 1;
    int k;
    for (int i = 1; i < h->GetNbinsX() + 1; i++) {
      integral = h->Integral(i, i, 1, h->GetNbinsY());
      for (int j = 1; j < h->GetNbinsY() + 1; j++) {
        k = h->GetBin(i, j);
        h->SetBinContent(k, h->GetBinContent(k)/integral);
        h->SetBinError(k, h->GetBinError(k)/integral);
      }
    }
}

void AnalysisZprime::WriteHistograms() {
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

    // this->Make2dDistribution(h2_mtt_deltaPhi);
    // this->Make2dDistribution(h2_mtt_cos1cos2);
    // this->Make2dDistribution(h2_HT_deltaPhi);
    // this->Make2dDistribution(h2_KT_deltaPhi);
    this->Make2dDistribution(h2_mtt_cosThetaStar);
    this->Make2dDistribution(h2_mtt_cosThetaStar_R);
    this->Make2dDistribution(h2_mtt_cosTheta1);
    this->Make2dDistribution(h2_mtt_cosTheta2);
    this->Make2dDistribution(h2_mtt_cosThetal_R);

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

bool AnalysisZprime::PassCuts(string type) {
    if(!(type == "truth" or type == "R1" or type == "R2")) return false;
    if(this->PassCutsET(type)) {
        if(this->PassCutsEta(type)) {
            if(this->PassCutsMET(type)) {
                if(this->PassCutsMtt(type)) {
                    if(this->PassCutsYtt(type)) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

bool AnalysisZprime::PassCutsMET(string type) {
    bool pass;
    pass = true;

    this->UpdateCutflow(c_MET, pass);
    return pass;
}

bool AnalysisZprime::PassCutsMtt(string type) {
    double mtt;
    if (type == "truth") mtt = P.M()/1000;
    else if (type == "R1") mtt = P_R1.M()/1000;
    else if (type == "R2") mtt = P_R2.M()/1000;
    else return false;

    if (m_Emin < 0 and m_Emax < 0) return true;

    if (mtt > m_Emin) {
        if (mtt < m_Emax) {
            UpdateCutflow(c_mtt, true);
            return true;
        }
    }
    UpdateCutflow(c_mtt, false);
    return false;
}

bool AnalysisZprime::PassCutsEta(string type) {
    if (m_fid == false) {
        UpdateCutflow(c_eta, true);
        return true;
    }
    if (type == "truth") {
        for (unsigned int i = 0; i < p.size(); i++) {
            // bool outsideCrack = abs(p[i].PseudoRapidity()) <= 1.37 || abs(p[i].PseudoRapidity()) >= 1.52;
            bool central = abs(p[i].PseudoRapidity()) <= 2.5;
            // bool passesEtaCuts = outsideCrack && central;
            if (central == false) {
                UpdateCutflow(c_eta, false);
                return false;
            }
            else continue;
        }
    }
    else if (type == "R1") {
        for (unsigned int i = 0; i < p_R1.size(); i++) {
            // bool outsideCrack = p_R1[i].PseudoRapidity() <= 1.37 || p_R1[i].PseudoRapidity() >= 1.52;
            bool central = p_R1[i].PseudoRapidity() <= 2.5;
            // bool passesEtaCuts = outsideCrack && central;
            if (central == false) {
                UpdateCutflow(c_eta, false);
                return false;
            }
            else continue;
        }
    }
    else if (type == "R2") {
        for (unsigned int i = 0; i < p_R2.size(); i++) {
            // bool outsideCrack = p_R2[i].PseudoRapidity() <= 1.37 || p_R2[i].PseudoRapidity() >= 1.52;
            bool central = p_R2[i].PseudoRapidity() <= 2.5;
            // bool passesEtaCuts = outsideCrack && central;
            if (central == false) {
                UpdateCutflow(c_eta, false);
                return false;
            }
            else continue;
        }
    }
    else return false;
    UpdateCutflow(c_eta, true);
    return true;
}

bool AnalysisZprime::PassCutsET(string type) {
    if (m_fid == false) {
        UpdateCutflow(c_Et, true);
        return true;
    }
    if (type == "truth") {
        for(unsigned int i = 0; i < p.size(); i++) {
            if(p[i].Et() <= 25) {
                UpdateCutflow(c_Et, false);
                return false;
            }
            else continue;
        }
    }
    else if (type == "R1") {
        for(unsigned int i = 0; i < p_R1.size(); i++) {
            if(p_R1[i].Et() <= 25) {
                UpdateCutflow(c_Et, false);
                return false;
            }
            else continue;
        }
    }
    else if (type == "R2") {
        for(unsigned int i = 0; i < p_R2.size(); i++) {
            if(p_R2[i].Et() <= 25) {
                UpdateCutflow(c_Et, false);
                return false;
            }
            else continue;
        }
    }
    else return false;
    UpdateCutflow(c_Et, true);
    return true;
}

bool AnalysisZprime::PassCutsYtt(string type) {
    double ytt;
    if (type == "truth") ytt = abs(P.Rapidity());
    else if (type == "R1") ytt = abs(P_R1.Rapidity());
    else if (type == "R2") ytt = abs(P_R2.Rapidity());
    else return false;
    if (ytt > m_ytt) {
        UpdateCutflow(c_ytt, true);
        return true;
    }
    UpdateCutflow(c_ytt, false);
    return false;
}

void AnalysisZprime::PreLoop() {
    this->GetDataDirectory();
    this->SetupInputFiles();
    this->SetupOutputFiles();
    this->ResetCounters();
    this->InitialiseCutflow();
    this->CreateHistograms();
}

void AnalysisZprime::GetDataDirectory() {
    // char hostname[HOST_NAME_MAX];
    // gethostname(hostname, HOST_NAME_MAX);
    // char *hostname = getenv("HOSTNAME");
    // printf("%s\n", hostname);
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    string Hostname(hostname);
    // printf("Hostname: %s\n", Hostname.c_str());

    if (Hostname == "Sunder")
        m_dataDirectory = "/Users/declan/Data/zprime";
    else if ((Hostname.find("lxplus") != std::string::npos) || (Hostname.find("cern") != std::string::npos))
        m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime";
    else if ((Hostname.find("cyan") != std::string::npos) || (Hostname.find("blue") != std::string::npos) || (Hostname.find("green") != std::string::npos))
        m_dataDirectory = "/scratch/dam1g09/zprime";
    else
        printf("Hostname %s not recognised.\n", Hostname.c_str());
    // #elif __APPLE__ || __MACH__
    // m_dataDirectory = "/Users/declan/Data/Zprime";
    // #endif
}

void AnalysisZprime::ResetCounters() {
    if (m_luminosity >= 0) m_useLumi = true;
    else m_useLumi = false;
    m_nQuarksMatched = 0;
    m_nNeutrinoMatched = 0;
    m_nReco = 0;
    m_nRealRoots = 0;
    m_nComplexRoots = 0;
}

void AnalysisZprime::GetCrossSection(TString log) {
    log.Replace(log.Last('.'), 5, ".log");
    printf("%s\n", log.Data());
    ifstream logstream(log.Data());
    if (!logstream.is_open()) printf("Error: failed to open %s!\n", log.Data());
    string line;
    string target = "Cross section";
    vector<string> parts;
    bool found = false;
    while(getline(logstream, line)) {
        trim(line);
        split(parts, line, is_any_of(":"));
        for(auto part: parts) {
            trim(part);
            // printf("%s\n", part.c_str());
        }
        if (parts[0] == target) {
            m_sigma = stod(parts[1]);
            found = true;
        }
    }
    logstream.close();
    if (!found) {
        printf("Error: Failed to read generation cross section. Check target log file: %s", log.Data());
        exit(1);
    }
    else printf("Generation Cross section = %.15le [pb]\n", m_sigma);
}

void AnalysisZprime::GetIterationWeights(TString log) {
    iteration_weights.clear();
    log.Replace(log.Last('.'), 5, ".log");
    ifstream logstream(log.Data());
    if (!logstream.is_open()) printf("Error: failed to open %s!\n", log.Data());
    string line;
    string target = "Iteration weighting";
    vector<string> parts;
    bool found = false;
    while(getline(logstream, line)) {
        trim(line);
        split(parts, line, is_any_of(":"));
        for(auto part: parts) {
            trim(part);
        }
        if (parts[0] == target) {
            iteration_weights.push_back(stod(parts[2]));
            found = true;
        }
    }
    logstream.close();
    if (!found) {
        printf("Error: Failed to read vegas iteration weights. Check target log file: %s", m_weightsFilename.Data());
        exit(1);
    }
}

void AnalysisZprime::Loop() {
    for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        printf("\n--- Input %li ---\n", i - m_inputFiles->begin() + 1);
        cout << (*i) << endl;
        this->SetupTreesForNewFile((*i));
        this->GetCrossSection(*i);
        this->GetIterationWeights(*i);
        this->GetChannelFactors();
        Long64_t nEntries;
        nEntries = this->TotalEvents();
        printf("Processing %lld entries...\n", nEntries);
        for (Long64_t jentry = 0; jentry < nEntries; ++jentry)
        {
            Long64_t ientry = this->IncrementEvent(jentry);
            if (ientry < 0) break;
            this->EachEvent();
            // this->ProgressBar(jentry, nEntries-1, 50);
        }
        this->CleanUp();
    }
    cout << endl;
}

AnalysisZprime::~AnalysisZprime() {delete m_inputFiles;}

Long64_t AnalysisZprime::TotalEvents() {
    if (m_ntup != 0) {return m_ntup->totalEvents();}
    return -999;
}

Long64_t AnalysisZprime::IncrementEvent(Long64_t i) {
    Long64_t ev(-1);
    if (m_ntup != 0) {ev = m_ntup->LoadTree(i);}
    return ev;
}

void AnalysisZprime::SetupTreesForNewFile(const TString& s) {
    TString treeToUse = "RootTuple";

    m_chainNtup = new TChain(treeToUse,"");
    TString TStringNtuple = s + "/" + treeToUse;
    m_chainNtup->Add(TStringNtuple,0);
    m_ntup = new RootTuple(m_chainNtup);
}

void AnalysisZprime::CleanUp() {
    delete m_chainNtup;
    delete m_ntup;
}

vector<TLorentzVector> AnalysisZprime::ReconstructSemiLeptonic(vector<TLorentzVector> p, int Q_l) {
    // Returns a vector of 4-momenta for all 6 particles in the final state with matching of b-quarks to each top
    // and matching of
    // Takes a vector of true final-state particle momenta as the argument and the charge of the final
    // state lepton: if +, t decayed leptonically; if -, t~ decayed leptonically.
    // As going from bbllnn->bblnqq/bbqqln requires only a simple reweighting for parton truth,
    // it saves on storage space and processing time to store all events as bbllnn. However,
    // when we reconstruct the neutrino, we must account for the fact either the top, or the anti-top
    // may decay hadronically. This means there are two distinguishable final states:
    // Q_l = +1 : pp -> b b~ l+ nu q q'
    // Q_l = -1 : pp -> b b~ q q' l- nu
    // Note that the order here is important, as the order of indicies in the vector of final state momenta
    // relates to the parent particle t=(0,2,3), t~=(1,4,5) and is fixed at the generator level.
    // If we want the results combining each final state, we must add these together.
    // Note: Experimentally p^{x,y}_nu is equated to the MET, of course.

    this->UpdateCutflow(c_events, true);

    m_nReco++;

    vector<TLorentzVector> p_nu_R, p_R(p.size());
    TLorentzVector p_l, p_nu;
    double a, b, c, k, dh, dl, mblv, mjjb, chi2, chi2min = 1.0e10;
    unsigned int imin, jmin;
    vector<complex<double> > roots;

    // Calculate neutrino p_z solutions
    if (Q_l == 1) {
        this->UpdateCutflow(c_topDecays, true);
        p_l = p[2];
        p_nu = p[3];
    }
    else if (Q_l == -1) {
        this->UpdateCutflow(c_antitopDecays, true);
        p_l = p[4];
        p_nu = p[5];
    }
    else{
        exit(1);
    }

    double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
    double px_nu = p_nu.Px(), py_nu = p_nu.Py();
    E_l = sqrt(px_l*px_l + py_l*py_l + pz_l*pz_l);
    if (abs(E_l - p_l.E()) > 0.00001) printf("ERROR: Lepton energy doesn't match.\n");

    k = 3218.42645 + px_l*px_nu + py_l*py_nu; // m_Wmass*m_Wmass/2 = 3218.42645
    a = px_l*px_l + py_l*py_l;
    b = -2*k*(pz_l);
    c = (px_nu*px_nu + py_nu*py_nu)*E_l*E_l - k*k;

    roots = this->SolveQuadratic(a, b, c);
    p_nu_R.clear();
    if (roots[0].imag() == 0 and roots[1].imag() == 0) {
        // two real solutions; pick best match
        this->UpdateCutflow(c_realSolutions, true);
        for (auto root: roots) {
            double pz = root.real();
            TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu*px_nu + py_nu*py_nu + pz*pz));
            p_nu_R.push_back(p);
        }
    }
    else{
        // no real solutions; take the real part of 1 (real parts are the same)
        if (m_discardComplex) m_discardEvent = true; // dump complex events
        double pz = roots[0].real();
        TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu*px_nu + py_nu*py_nu + pz*pz));
        p_nu_R.push_back(p);
    }
    bool oldreco = false;

    if (oldreco) {
        vector<TLorentzVector> p_b(2), p_q(2);

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
                chi2 = dh*dh + dl*dl;
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
    else{
        vector<TLorentzVector> p_q(4);

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
        vector<vector<int> > q_perms;
        if (nBtags == 2) q_perms = { {0, 1, 2, 3}, {1, 0, 2, 3} };
        if (nBtags == 1) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
        {1, 0, 2, 3}, {2, 0, 1, 3}, {3, 0, 1, 2} };
        if (nBtags == 0) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
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
    unsigned int b_lep;
    if (Q_l == 1) b_lep = 0;
    if (Q_l == -1) b_lep = 1;
    if (b_lep == jmin) m_nQuarksMatched++;

    // Assess neutrino reconstruction performance.
    double pz_nu_truth = p_nu.Pz();
    double Root0MinusTruth = abs(roots[0].real() - pz_nu_truth);
    double Root1MinusTruth = abs(roots[1].real() - pz_nu_truth);
    unsigned int bestRoot;
    if (Root0MinusTruth < Root1MinusTruth) bestRoot = 0;
    else if (Root1MinusTruth < Root0MinusTruth) bestRoot = 1;
    else bestRoot = 0;
    if (imin == bestRoot) m_nNeutrinoMatched++;
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
    // printf("---\n");
    return p_R;
}

vector<complex<double> > AnalysisZprime::SolveQuadratic(double a, double b, double c) {
    vector<complex<double> > roots;
    complex<double> term1;
    complex<double> term2;
    complex<double> discriminator;

    term1 = -b/(2*a);
    discriminator = b*b - 4*a*c;
    term2 = sqrt(discriminator)/(2*a);

    roots.push_back(term1 + term2);
    roots.push_back(term1 - term2);

    return roots;
}

bool solveQuadratic(double a, double b, double c, double &root) {
    if (a == 0.0 || abs(a/b) < 1.0e-4) {
        if (abs(b) < 1.0e-4) 
            return false;
        else {
            root = -c/b;
            return true;
        }
    }

    double discriminant = b*b - 4.0*a*c;
    if (discriminant >= 0.0) {
        discriminant = sqrt(discriminant);
        root = (b - discriminant)*-0.5/a;
        return true;
    }

    return false;
}

bool solveCubic(double a, double b, double c, double d, double &root) {
    if (a == 0.0 || abs(a/b) < 1.0e-6) return solveQuadratic(b, c, d, root);

    double B = b/a, C = c/a, D = d/a;
    double Q = (B*B - C*3.0)/9.0, QQQ = Q*Q*Q;
    double R = (2.0*B*B*B - 9.0*B*C + 27.0*D)/54.0, RR = R*R;
    double pi = 3.14159265;

    // 3 real roots
    if (RR < QQQ) {
        // This sqrt and division is safe, since RR >= 0, so QQQ > RR,    
        // so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and     
        // thus R/sqrt(QQQ) < 1.                                     
        double theta = acos(R/sqrt(QQQ));
        // This sqrt is safe, since QQQ >= 0, and thus Q >= 0

        double r1, r2, r3;
        r1 = r2 = r3 = -2.0*sqrt(Q);
        r1 *= cos(theta/3.0);
        r2 *= cos((theta + 2*pi)/3.0);
        r3 *= cos((theta - 2*pi)/3.0);

        r1 -= B/3.0;
        r2 -= B/3.0;
        r3 -= B/3.0; 

        root = 1000000.0;

        if (r1 >= 0.0) root = r1;
        if (r2 >= 0.0 && r2 < root) root = r2;
        if (r3 >= 0.0 && r3 < root) root = r3;

        return true;
    }
    // 1 real root
    else {
        double A2 = -pow(fabs(R) + sqrt(RR - QQQ), 1.0/3.0);
        if (A2 != 0.0) {
            if (R < 0.0) A2 = -A2; 
            root = A2 + Q/A2; 
        }
        root -= B/3.0;
        return true;
    }
}

bool solveQuartic(double a, double b, double c, double d, double e, double &root) {
    // When a or (a and b) are magnitudes of order smaller than C, D, E just ignore them entirely. 
    if (a == 0.0 || abs(a/b) < 1.0e-5 || abs(a/c) < 1.0e-5 || abs(a/d) < 1.0e-5) return solveCubic(b, c, d, e, root);

    double B = b/a, C = c/a, D = d/a, E = e/a;
    double BB = B*B;
    double I = -3.0*BB*0.125 + C;
    double J = BB*B*0.125 - B*C*0.5 + D;
    double K = -3*BB*BB/256.0 + C*BB/16.0 - B*D*0.25 + E;
    double z;
    // bool foundRoot2 = false, foundRoot3 = false, foundRoot4 = false, foundRoot5 = false;
    if (solveCubic(1.0, 2*I, I*I - 4*K, -J*J, z)) {
        // double value = z*z*z + 2*z*z*I + z*(I*I - 4*K) - J*J;

        double p = sqrt(z);
        double r = -p;
        double q = (I + z - J/p)*0.5;
        double s = (I + z + J/p)*0.5;

        bool foundRoot = false, foundARoot;
        double aRoot;
        foundRoot = solveQuadratic(1.0, p, q, root);
        root -= B/4.0;

        foundARoot = solveQuadratic(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) || root < 0.0)) || (!foundRoot && foundARoot)) {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadratic(1.0, p, q, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) || root < 0.0)) || (!foundRoot && foundARoot)) {
            root = aRoot;
            foundRoot = true;
        }

        foundARoot = solveQuadratic(1.0, r, s, aRoot);
        aRoot -= B/4.0;
        if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0) || root < 0.0)) || (!foundRoot && foundARoot)) {
            root = aRoot;
            foundRoot = true;
        }
        return foundRoot;
    }
    return false;
}

vector<TLorentzVector> AnalysisZprime::ReconstructDilepton(vector<TLorentzVector> p) {

    this->UpdateCutflow(c_events, true);

    m_nReco++;

    vector<TLorentzVector> p_R(p.size());

    vector<complex<double> > roots;

    TLorentzVector pb1 = p[0], pb2 = p[1], pl1 = p[2], pv1 = p[3], pl2 = p[4], pv2 = p[5];

    double pb1x = pb1.Px(), pb1y = pb1.Py(), pb1z = pb1.Pz(), Eb1 = pb1.E();
    double pb2x = pb2.Px(), pb2y = pb2.Py(), pb2z = pb2.Pz(), Eb2 = pb2.E();
    double pl1x = pl1.Px(), pl1y = pl1.Py(), pl1z = pl1.Pz(), El1 = pl1.E();
    double pl2x = pl2.Px(), pl2y = pl2.Py(), pl2z = pl2.Pz(), El2 = pl2.E();
    double pv1x = pv1.Px(), pv1y = pv1.Py(), pv2x = pv2.Px(), pv2y = pv2.Py();
    double Emissx = pv1x + pv2x;
    double Emissy = pv1y + pv2y;

    double ml1 = 0, ml2 = 0, mv1 = 0, mv2 = 0;

    double a1 = (Eb1 + El1)*(m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)
                - El1*(m_tmass*m_tmass - m_bmass*m_bmass - ml1*ml1 - mv1*mv1)
                + 2*Eb1*El1*El1 - 2*El1*(pb1x*pl1x + pb1y*pl1y + pb1z*pl1z);
    double a2 = 2*(Eb1*pl1x - El1*pb1x);
    double a3 = 2*(Eb1*pl1y - El1*pb1y);
    double a4 = 2*(Eb1*pl1z - El1*pb1z);
    double b1 = (Eb2 + El2)*(m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2)
                - El2*(m_tmass*m_tmass - m_bmass*m_bmass - ml2*ml2 - mv2*mv2)
                + 2*Eb2*El2*El2 - 2*El2*(pb2x*pl2x + pb2y*pl2y + pb2z*pl2z);                
    double b2 = 2*(Eb2*pl2x - El2*pb2x);
    double b3 = 2*(Eb2*pl2y - El2*pb2y);
    double b4 = 2*(Eb2*pl2z - El2*pb2z);
    double c22 = (m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)*(m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)
                 - 4*(El1*El1 - pl1z*pl1z)*a1*a1/a4/a4
                 - 4*(m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)*pl1z*a1/a4;
    double c21 = 4*(m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)*(pl1x - pl1z*a2/a4)
                 - 8*(El1*El1 - pl1z*pl1z)*a1*a2/a4/a4 
                 - 8*pl1x*pl1z*a1/a4;
    double c20 = -4*(El1*El1 - pl1x*pl1x)
                 - 4*(El1*El1 - pl1z*pl1z)*a2*a2/a4/a4
                 - 8*pl1x*pl1z*a2/a4;
    double c11 = 4*(m_Wmass*m_Wmass - ml1*ml1 - mv1*mv1)*(pl1y - pl1z*a3/a4)
                 - 8*(El1*El1 - pl1z*pl1z)*a1*a3/a4/a4
                 - 8*pl1y*pl1z*a1/a4;
    double c10 = -8*(El1*El1 - pl1z*pl1z)*a2*a3/a4/a4 + 8*pl1x*pl1y
                 - 8*pl1x*pl1z*a3/a4 - 8*pl1y*pl1z*a2/a4;
    double c00 = -4*(El1*El1 - pl1y*pl1y) - 4*(El1*El1 - pl1z*pl1z)*a3*a3/a4/a4
                 - 8*pl1y*pl1z*a3/a4;
    double dd22 = (m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2)*(m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2) - 4*(El2*El2 - pl1z*pl1z)*b1*b1/b4/b4
                 - 4*(m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2)*pl1z*b1/b4;
    double dd21 = 4*(m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2)*(pl1x - pl1z*b2/b4)
                 - 8*(El2*El2 - pl1z*pl1z)*b1*b2/b4/b4 
                 - 8*pl1x*pl1z*b1/b4;
    double d20 = -4*(El2*El2 - pl1x*pl1x)
                 - 4*(El2*El2 - pl1z*pl1z)*b2*b2/b4/b4
                 - 8*pl1x*pl1z*b2/b4; 
    double dd11 = 4*(m_Wmass*m_Wmass - ml2*ml2 - mv2*mv2)*(pl1y - pl1z*b3/b4)
                 - 8*(El2*El2 - pl1z*pl1z)*b1*b3/b4/b4
                 - 8*pl1y*pl1z*b1/b4;
    double d10 = -8*(El2*El2 - pl1z*pl1z)*b2*b3/b4/b4 + 8*pl1x*pl1y
                 - 8*pl1x*pl1z*b3/b4 - 8*pl1y*pl1z*b2/b4;
    double d00 = -4*(El2*El2 - pl1y*pl1y) - 4*(El2*El2 - pl1z*pl1z)*b3*b3/b4/b4
                 - 8*pl1y*pl1z*b3/b4;
    double d22 = dd22 + Emissx*Emissx*d20 +Emissy*Emissy*d00
                 + Emissx*Emissy*d10 + Emissx*dd21 + Emissy*dd11;
    double d21 = -dd21 -2*Emissx*d20 - Emissy*d10;
    double d11 = -dd11 -2*Emissy*d00 - Emissx*d10; 
    double h4 = c00*c00*d22*d22 + c11*d22*(c11*d00 - c00*d11)
                + c00*c22*(d11*d11 - 2*d00*d22) + c22*d00*(c22*d00 - c11*d11);
    double h3 = c00*d21*(2*c00*d22 - c11*d11) + c00*d11*(2*c22*d10 + c21*d11) 
                + c22*d00*(2*c21*d00-c11*d10)-c00*d22*(c11*d10 + c10*d11)  
                -2*c00*d00*(c22*d21 + c21*d22)-d00*d11*(c11*c21 + c10*c22) 
                + c11*d00*(c11*d21 + 2*c10*d22);
    double h2 = c00*c00*(2*d22*d20 + d21*d21) - c00*d21*(c11*d10 + c10*d11)  
                + c11*d20*(c11*d00 - c00*d11) + c00*d10*(c22*d10 - c10*d22)   
                + c00*d11*(2*c21*d10 + c20*d11) + (2*c22*c20 + c21*c21)*d00*d00   
                -2*c00*d00*(c22*d20 + c21*d21 + c20*d22)    
                + c10*d00*(2*c11*d21 + c10*d22) - d00*d10*(c11*c21 + c10*c22)   
                - d00*d11*(c11*c20 + c10*c21);
    double h1 = c00*d21*(2*c00*d20 - c10*d10) - c00*d20*(c11*d10 + c10*d11)  
                + c00*d10*(c21*d10 + 2*c20*d11) - 2*c00*d00*(c21*d20 + c20*d21)  
                + c10*d00*(2*c11*d20 + c10*d21) - c20*d00*(2*c21*d00 - c10*d11)  
                - d00*d10*(c11*c20 + c10*c21);
    double h0 = c00*c00*d20*d20 + c10*d20*(c10*d00 - c00*d10)  
                + c20*d10*(c00*d10 - c10*d00) + c20*d00*(c20*d00 - 2*c00*d20);

    double pv1x_R;
    solveQuartic(h4, h3, h2, h1, h0, pv1x_R);
    double pv2x_R = Emissx - pv1x_R;

    double c2 = c22 + c21*pv1x_R + c20*pv1x_R*pv1x_R;
    double c1 = c11 + c10*pv1x_R;
    double c0 = c00;
    double d2 = d22 + d21*pv1x_R + d20*pv1x_R*pv1x_R;
    double d1 = d11 + d10*pv1x_R;
    double d0 = c00;

    double pv1y_R = (c0*d2 - c2*d0)/(c1*d0 - c0*d1);
    double pv2y_R = Emissy - pv1y_R;

    double pv1z_R = -(a1 + a2*pv1x_R + a3*pv1y_R)/a4;
    double pv2z_R = -(b1 + b2*pv2x_R + b3*pv2y_R)/b4;

    TLorentzVector pv1_R(pv1x_R, pv1y_R, pv1z_R, sqrt(pv1x_R*pv1x_R + pv1y_R*pv1y_R + pv1z_R*pv1z_R));
    TLorentzVector pv2_R(pv2x_R, pv2y_R, pv2z_R, sqrt(pv2x_R*pv2x_R + pv2y_R*pv2y_R + pv2z_R*pv2z_R));

    p_R[0] = pb1;
    p_R[1] = pb2;
    p_R[2] = pl1;
    p_R[3] = pv1_R;
    p_R[4] = pl2;
    p_R[5] = pv2_R;

    return p_R;
}

void AnalysisZprime::GetChannelFactors() {
    // scale dilepton to other classifications
    // fac_ee = 1
    // fac_emu = 2
    // m_sigma = m_sigma*24; // 2 [e+ + e-] x 2 [e + mu] x 6 [3 x 2]
    // fac_qq = 36
}

const void AnalysisZprime::UpdateCutflow(int cut, bool passed) {
    if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
    if (passed) m_cutflow[cut] +=1;
}

void AnalysisZprime::InitialiseCutflow() {
    m_cutflow = vector<int>(m_cuts, -999);
    m_cutNames = vector<TString>(m_cuts, "no name");
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

void AnalysisZprime::PrintCutflow() {
    printf("--- Cutflow ---\n");
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut]);

        printf("%s %i pass\n", m_cutNames[cut].Data(), m_cutflow[cut]);
    }
    h_cutflow->Write();
}

inline void AnalysisZprime::ProgressBar(unsigned int x, unsigned int n, unsigned int w) {
    if ( (x != n) && (x % (n/100+1) != 0) ) return;

    float ratio = x/(float)n;
    unsigned int c = ratio * w;

    cout << setw(3) << (int)(ratio*100) << "% [";
    for (unsigned int i = 0; i < c; i++) cout << "=";
    for (unsigned int i = c; i < w; i++) cout << " ";
    if (x == n) cout << "\n" << flush;
    else cout << "]\r" << flush;
}

void AnalysisZprime::SetYttCut(const double ytt){
    m_ytt = ytt;
}

void AnalysisZprime::SetXsec(const bool xsec){
    m_xsec = xsec;
}

void AnalysisZprime::SetFiducial(const bool fid){
    m_fid = fid;
}
