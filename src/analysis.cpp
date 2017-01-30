#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"

Analysis::Analysis(const TString& model, const TString& process, const TString& options, const int energy, const int luminosity, const int reco, const TString tag):
    m_model(model),
    m_process(process),
    m_options(options),
    m_energy(energy),
    m_luminosity(luminosity),
    m_reco(reco), 
    m_tag(tag), 
    m_outputFile(nullptr),
    m_inputFiles(nullptr),
    m_weightFiles(nullptr),
    m_chainNtup(nullptr),
    m_ntup(nullptr)
{ 
    this->PreLoop();
    this->Loop();
    this->PostLoop();
}


void Analysis::EachEvent()
{
    UpdateCutflow(c_entries, true);

    // set event weight
    double weight = m_ntup->weight();
    if (weight == 0) return;
    if (m_xsec) {
        weight = weight * m_sigma / m_nevents;
        if (m_useLumi) weight = weight * m_efficiency * m_luminosity;
    }
    double weight_R = weight;
    if (m_reco == 2) weight_R = weight_R / 2;

    // get number of particles in final state
    n = m_ntup->barcode()->size();

    // no reconstruction for 2 -> 2
    if (n < 6) m_reco = 0;

    std::vector<TLorentzVector> p(n);
    for (int i = 0; i < n; i++)
        p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));

    TLorentzVector P(0, 0, 0, 0);
    for (auto& l : p) P += l;

    if (this->PassFiducialCuts(p, P)) {

        std::vector<TLorentzVector> pcm = p;
        TVector3 v = -1 * P.BoostVector();
        for (auto& l : pcm) l.Boost(v);

        double mtt = P.M() / 1000;
        double ytt = P.Rapidity();

        // top and antitop
        TLorentzVector p_t;
        TLorentzVector p_tbar;
        TLorentzVector pcm_t;
        TLorentzVector pcm_tbar;

        if (n < 6) {
            p_t = p[0];
            p_tbar = p[1];
            pcm_t = pcm[0];
            pcm_tbar = pcm[1];
        }
        else {
            p_t = p[0] + p[2] + p[3];
            p_tbar = p[1] + p[4] + p[5];
            pcm_t = pcm[0] + pcm[2] + pcm[3];
            pcm_tbar = pcm[1] + pcm[4] + pcm[5];
        }

        double costheta_tt = cos(p_t.Angle(p_tbar.Vect()));
        double delta_abs_yt = std::abs(p_t.Rapidity()) - std::abs(p_tbar.Rapidity());
        double costheta = pcm_t.CosTheta();
        double costhetastar = int(ytt / std::abs(ytt)) * costheta;
        std::vector<TLorentzVector> p_R1(n), pcm_R1(n), ptop_R1(n), patop_R1(n);
        std::vector<TLorentzVector> p_R2(n), pcm_R2(n), ptop_R2(n), patop_R2(n);


        if (this->PassCuts(p, P)) {

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

            h_mtt->Fill(mtt, weight);
            h_ytt->Fill(ytt, weight);
            h_costheta_tt->Fill(costheta_tt, weight_R);
            h_cosTheta->Fill(costheta, weight);
            h_cosThetaStar->Fill(costhetastar, weight);

            if (costhetastar > 0) h_mtt_tF->Fill(mtt, weight);
            if (costhetastar < 0) h_mtt_tB->Fill(mtt, weight);

            if (delta_abs_yt > 0) h_mtt_tCF->Fill(mtt, weight);
            if (delta_abs_yt < 0) h_mtt_tCB->Fill(mtt, weight);

            // if (n < 6) {
                // h_mtt_LL->Fill(mtt, m_ntup->weightLL());
                // h_mtt_LR->Fill(mtt, m_ntup->weightLR());
                // h_mtt_RL->Fill(mtt, m_ntup->weightRL());
                // h_mtt_RR->Fill(mtt, m_ntup->weightRR());
            // }

            if (n == 6) {

                std::vector<TLorentzVector> ptop = p;
                v = -1 * (p[0] + p[2] + p[3]).BoostVector();
                for (auto& l : ptop) l.Boost(v);

                std::vector<TLorentzVector> patop = p;
                v = -1 * (p[1] + p[4] + p[5]).BoostVector();
                for (auto& l : patop) l.Boost(v);

                double HT = 0;
                for (auto& l : p) HT += l.Pt();
                HT = HT / 1000;
                TLorentzVector pvis = p[0] + p[1] + p[2] + p[4];
                double mvis = pvis.M();
                double pTvis = pvis.Pt();
                double KT = sqrt(mvis * mvis + pTvis * pTvis) + (p[3] + p[5]).Pt();
                KT = KT / 1000;

                std::vector<double> deltaRs;
                for (int i = 0; i < 6; i++)
                    for (int j = i + 1; j < 6; j++)
                        deltaRs.push_back(p[i].DeltaR(p[j]));

                double costheta_tl1 = cos(ptop[2].Angle(pcm_t.Vect()));
                double costheta_tl2 = cos(patop[4].Angle(pcm_tbar.Vect()));
                double cos1cos2 = costheta_tl1 * costheta_tl2;
                double deltaPhi = p[2].DeltaPhi(p[4]) / m_pi;

                TLorentzVector p_W = p[2] + p[3];
                double deltaR_bW = p_W.DeltaR(p[0]);
                std::vector<double> deltaR_ts;
                deltaR_ts.push_back(p_t.DeltaR(p[0]));
                deltaR_ts.push_back(p_t.DeltaR(p[2]));
                deltaR_ts.push_back(p_t.DeltaR(p[3]));
                auto deltaR_max = max_element(std::begin(deltaR_ts), std::end(deltaR_ts));

                double phil = p[2].Phi();
                double El = p[2].E();

                double phi0 = 0.7, E0 = 80;


                h_HT->Fill(HT, weight);
                h_KT->Fill(KT, weight);
                h_pv1x->Fill(p[3].Px(), weight);
                h_pv1y->Fill(p[3].Py(), weight);
                h_pv1z->Fill(p[3].Pz(), weight);
                h_pv2x->Fill(p[5].Px(), weight);
                h_pv2y->Fill(p[5].Py(), weight);
                h_pv2z->Fill(p[5].Pz(), weight);
                h_deltaPhi->Fill(deltaPhi, weight);
                h_cosTheta1->Fill(costheta_tl1, weight);
                h_cosTheta2->Fill(costheta_tl2, weight);
                h_cos1cos2->Fill(cos1cos2, weight);
                h_deltaR_bW->Fill(deltaR_bW, weight);
                h_deltaR_max->Fill(*deltaR_max, weight);

                for (int i = 0; i < (int) deltaRs.size(); i++)
                    h_deltaRs[i]->Fill(deltaRs[i], weight);

                for (int i = 0; i < (int) p.size(); i++) {
                    h_eta[i]->Fill(p[i].Eta(), weight);
                    h_pt[i]->Fill(p[i].Pt(), weight);
                }

                if (costheta_tl1 > 0) h_mtt_tlF->Fill(mtt, weight);
                if (costheta_tl2 < 0) h_mtt_tlB->Fill(mtt, weight);

                if (costheta_tl1 > 0) h_mtt_lF->Fill(mtt, weight);
                if (costheta_tl1 < 0) h_mtt_lB->Fill(mtt, weight);

                if (phil < m_pi and phil > phi0) h_mtt_philF->Fill(mtt, weight);
                if (phil < phi0) h_mtt_philB->Fill(mtt, weight);

                if (El > E0) h_mtt_ElF->Fill(mtt, weight);
                if (El < E0) h_mtt_ElB->Fill(mtt, weight);

                h2_mtt_cosThetaStar->Fill(mtt, costhetastar, weight);
                h2_mtt_deltaPhi->Fill(mtt, deltaPhi, weight);
                h2_mtt_cosTheta1->Fill(mtt, costheta_tl1, weight);
                h2_mtt_cosTheta2->Fill(mtt, costheta_tl2, weight);
                h2_mtt_cos1cos2->Fill(mtt, cos1cos2, weight);
                h2_HT_deltaPhi->Fill(HT, deltaPhi, weight);
                h2_KT_deltaPhi->Fill(KT, deltaPhi, weight);


                TLorentzVector P_R1, P_R2;
                if (m_reco == 1) {
                    p_R1 = this->ReconstructDilepton(p); // both decay leptonically

                    int i = 0;
                    for (auto& l : p_R1) {
                        if (l != l) {
                            printf("p[%i] has a NaN!\n", i);
                            l.Print();
                        }
                        i++;
                    }

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
                
                TLorentzVector p_t_R1, p_tbar_R1, p_t_R2, p_tbar_R2;
                TLorentzVector pcm_t_R1, pcm_tbar_R1, pcm_t_R2, pcm_tbar_R2;
                if (m_reco > 0) {
                    p_t_R1 = p_R1[0] + p_R1[2] + p_R1[3];
                    p_tbar_R1 = p_R1[1] + p_R1[4] + p_R1[5];
                    pcm_t_R1 = pcm_R1[0] + pcm_R1[2] + pcm_R1[3];
                    pcm_tbar_R1 = pcm_R1[1] + pcm_R1[4] + pcm_R1[5];
                }
                if (m_reco == 2) {
                    pcm_t_R2 = pcm_R2[0] + pcm_R2[2] + pcm_R2[3];
                    pcm_tbar_R2 = pcm_R2[1] + pcm_R2[4] + pcm_R2[5];
                }
            
                double mtt_R1 = P_R1.M() / 1000;
                double mtt_R2 = P_R2.M() / 1000;

                double ytt_R1, ytt_R2;
                double costheta_R1, costheta_R2;
                double costhetastar_R1, costhetastar_R2;
                double costheta_tl_R1, costheta_tl_R2;
                double costheta_tt_R1;
                double delta_abs_yt_R1;

                double pT_t_perf, pT_tbar_perf, mttperf, costhetatt_perf, costhetastar_perf, costhetatl_perf;

                if (m_reco > 0) {
                    pT_t_perf = (p_t.Pt() - p_t_R1.Pt()) / p_t.Pt();
                    pT_tbar_perf = (p_tbar.Pt() - p_tbar_R1.Pt()) / p_tbar.Pt();
                    delta_abs_yt_R1 = std::abs(p_t_R1.Rapidity()) - std::abs(p_tbar_R1.Rapidity());
                    ytt_R1 = P_R1.Rapidity();
                    costheta_R1 = pcm_t_R1.CosTheta();
                    costhetastar_R1 = int(ytt_R1 / std::abs(ytt_R1)) * costheta_R1;
                    costheta_tl_R1 = cos(ptop_R1[2].Angle(pcm_t_R1.Vect()));
                    costheta_tt_R1 = cos(p_t_R1.Angle(p_tbar_R1.Vect()));
                    costhetatt_perf = (costheta_tt - costheta_tt_R1) / costheta_tt;
                    costhetastar_perf = (costhetastar - costhetastar_R1) / costhetastar;
                    costhetatl_perf = (costheta_tl1 - costheta_tl_R1) / costheta_tl1;
                    mttperf = (mtt - mtt_R1) / mtt;
                }
                if (m_reco == 2) {
                    ytt_R2 = P_R2.Rapidity();
                    costheta_R2 = pcm_t_R2.CosTheta();
                    costhetastar_R2 = int(ytt_R2 / std::abs(ytt_R2)) * costheta_R2;
                    costheta_tl_R2 = cos(patop_R2[4].Angle(pcm_tbar_R2.Vect()));
                }


                if (m_reco > 0) {
                    if (true) {//(this->PassCuts(p_R1, P_R1)) {

                        h_mtt_R->Fill(mtt_R1, weight_R);
                        h_ytt_R->Fill(ytt_R1, weight_R);

                        h_pxt_R->Fill(p_t_R1.Px(), weight_R);
                        h_pyt_R->Fill(p_t_R1.Py(), weight_R);
                        h_pzt_R->Fill(p_t_R1.Pz(), weight_R);
                        h_Et_R->Fill(p_t_R1.E(), weight_R);
                        h_pTt_R->Fill(p_t_R1.Pt(), weight_R);
                        h_etat_R->Fill(p_t_R1.Eta(), weight_R);
                        h_phit_R->Fill(p_t_R1.Phi() / m_pi, weight_R);
                        h_mt_R->Fill(p_t_R1.M(), weight_R);
                        h_costhetatt_perf->Fill(costhetatt_perf, weight);
                        h_costhetastar_perf->Fill(costhetastar_perf, weight);
                        h_costhetatl_perf->Fill(costhetatl_perf, weight);
                        h_m_tt_perf->Fill(mttperf, weight);
                        h_pT_t_perf->Fill(pT_t_perf, weight_R);
                        h_pT_tbar_perf->Fill(pT_tbar_perf, weight_R);
                        h2_m_tt_pT_t_perf->Fill(mttperf, pT_t_perf, weight_R);
                        h2_m_tt_costheta_tt_perf->Fill(mttperf, costhetatt_perf, weight);


                        h_pxtbar_R->Fill(p_tbar_R1.Px(), weight_R);
                        h_pytbar_R->Fill(p_tbar_R1.Py(), weight_R);
                        h_pztbar_R->Fill(p_tbar_R1.Pz(), weight_R);
                        h_Etbar_R->Fill(p_tbar_R1.E(), weight_R);
                        h_pTtbar_R->Fill(p_tbar_R1.Pt(), weight_R);
                        h_etatbar_R->Fill(p_tbar_R1.Eta(), weight_R);
                        h_phitbar_R->Fill(p_tbar_R1.Phi() / m_pi, weight_R);
                        h_mtbar_R->Fill(p_tbar_R1.M(), weight_R);

                        h_costheta_tt_R->Fill(costheta_tt_R1, weight_R);
                        h_cosTheta_R->Fill(costheta_R1, weight_R);
                        h_cosThetaStar_R->Fill(costhetastar_R1, weight_R);
                        h_pv1x_R->Fill(p_R1[3].Px(), weight_R);
                        h_pv1y_R->Fill(p_R1[3].Py(), weight_R);
                        h_pv1z_R->Fill(p_R1[3].Pz(), weight_R);
                        h_pv2x_R->Fill(p_R1[5].Px(), weight_R);
                        h_pv2y_R->Fill(p_R1[5].Py(), weight_R);
                        h_pv2z_R->Fill(p_R1[5].Pz(), weight_R);

                        if (costhetastar_R1 > 0) h_mtt_tF_R->Fill(mtt_R1, weight_R);
                        if (costhetastar_R1 < 0) h_mtt_tB_R->Fill(mtt_R1, weight_R);

                        if (delta_abs_yt > 0) h_mtt_tCF->Fill(mtt_R1, weight);
                        if (delta_abs_yt < 0) h_mtt_tCB->Fill(mtt_R1, weight);

                        if (costheta_tl1 > 0) h_mtt_tlF->Fill(mtt_R1, weight);
                        if (costheta_tl2 < 0) h_mtt_tlB->Fill(mtt_R1, weight);

                        if (costheta_tl1 > 0) h_mtt_lF->Fill(mtt_R1, weight);
                        if (costheta_tl1 < 0) h_mtt_lB->Fill(mtt_R1, weight);

                        h2_mtt_cosThetaStar_R->Fill(mtt_R1, costhetastar_R1, weight_R);
                        h2_mtt_cosThetal_R->Fill(mtt_R1, costheta_tl_R1, weight_R);
                    }
                }
                if (m_reco == 2) {
                    if (true) {//(this->PassCuts(p_R2, P_R2)) {
                        h_mtt_R->Fill(mtt_R2, weight_R);
                        h_ytt_R->Fill(ytt_R2, weight_R);
                        h_mt_R->Fill(pcm_t_R2.M(), weight_R);
                        h_mtbar_R->Fill(pcm_tbar_R2.M(), weight_R);
                        h_cosTheta_R ->Fill(costheta_R2, weight_R);
                        h_cosThetaStar_R->Fill(costhetastar_R2, weight_R);
                        if (costhetastar_R2 > 0) h_mtt_tF_R->Fill(mtt_R2, weight_R);
                        if (costhetastar_R2 < 0) h_mtt_tB_R->Fill(mtt_R2, weight_R);
                        h2_mtt_cosThetaStar_R->Fill(mtt_R2, costhetastar_R2, weight_R);
                        h2_mtt_cosThetal_R->Fill(mtt_R2, costheta_tl_R2, weight_R);
                    }
                }
            }
        }
    }
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
            // if (false);
            // if (initial == "gg") options = ".2-4";
            // else if (initial == "qq") options = ".2-4";
            options = m_options;
            std::string intermediates = "";
            if (initial == "uu" || initial == "dd") {
                intermediates = intermediates + "-AZ";
                if (m_model != "SM") intermediates = intermediates + "X";
            }
            filename = m_dataDirectory + "/" + initial + intermediates + final + "." + model + "." + E + "TeV" + "." + m_pdf + options;
            m_inputFiles->push_back(filename + ".root");
            m_weightFiles->push_back(filename + ".log");
        }
    }

    // m_inputFiles->push_back(m_dataDirectory + "/" + "dduu-AZ-tt-bbllvv.SM.13TeV.CT14LL.2-4.root");
    // m_weightFiles->push_back(m_dataDirectory + "/" + "dduu-AZ-tt-bbllvv.SM.13TeV.CT14LL.2-4.log");

    // check all input files exist
    for (auto inputFile : *m_inputFiles) {
        struct stat buffer;
        bool exists = stat(inputFile.Data(), &buffer) == 0;
        if (exists == false) {
            std::cout << "Error: no" << inputFile.Data() << std::endl;
            exit(exists);
        }
    }
}


void Analysis::SetupOutputFiles()
{
    std::string E = std::to_string(m_energy);

    m_outputFilename = m_dataDirectory + "/" + m_process + "." + m_model + "." + E + "TeV" + "." + m_pdf + m_options;
    m_outputFilename = m_outputFilename + ".r" + std::to_string(m_reco) + m_tag;

    if (m_reco == 2 && m_btags != 2) m_outputFilename += ".b" + std::to_string(m_btags);
    string ytt = std::to_string(m_ytt);
    if (m_ytt > 0) m_outputFilename += ".y" + ytt.erase(ytt.find_last_not_of('0') + 1, string::npos);
    if (m_Emin >= 0 || m_Emax >= 0) m_outputFilename += ".E" + std::to_string(m_Emin) + "-" + std::to_string(m_Emax);
    string eff = std::to_string(m_efficiency);
    if (m_efficiency < 1.0) m_outputFilename += ".e" + eff.erase(eff.find_last_not_of('0') + 1, string::npos);
    if (m_luminosity > 0) m_outputFilename += ".L" + std::to_string(m_luminosity);
    if (m_fid) m_outputFilename += ".fid";
    m_outputFilename += ".root";
    std::cout << "output: " << m_outputFilename.Data() << std::endl;
    m_outputFile = new TFile(m_outputFilename, "RECREATE");
}

void Analysis::PostLoop()
{
    this->CheckResults();
    if (n < 6) this->TotalSpinAsymmetries();
    if (m_reco == 2) this->CheckPerformance();
    this->MakeDistributions();
    this->WriteHistograms();
    this->PrintCutflow();
}


void Analysis::CheckResults()
{
    double sigma = h_mtt->Integral("width");
    std::cout << "analysis cross section: " << sigma << " [fb]" << std::endl;
}


void Analysis::CheckPerformance()
{
    double quarkRecoRatio = m_nQuarksMatched / (double) m_nReco;
    double neutrinoRecoRatio = m_nNeutrinoMatched / (double) m_nReco;
    printf("-- Performance --\n");
    printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio * 100);
    printf("pzNu assignment: %.1f%% correct\n", neutrinoRecoRatio * 100);
    printf("\n");
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
    double Emin = 2.05;
    double Emax = 3.95;
    double nbins = (Emax - Emin) / binWidth;
    std:: cout << "energy range: " << Emin << " to " << Emax << " [TeV]" << std::endl;

    h_mtt = new TH1D("m_tt", "m_{tt}", nbins, Emin, Emax);
    h_mtt->Sumw2();
    h_ytt = new TH1D("y_tt", "y_{tt}", 50, -2.5, 2.5);
    h_ytt->Sumw2();

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

    h2_mtt_cosThetaStar = new TH2D("mtt_costhetastar", "m_{tt} cos#theta^{*}", nbins, Emin, Emax, 2, -1.0, 1.0);
    h2_mtt_cosThetaStar->GetXaxis()->SetTitle("m_{tt}");
    h2_mtt_cosThetaStar->GetYaxis()->SetTitle("cos#theta^*");
    h2_mtt_cosThetaStar->Sumw2();

    h2_mtt_cosThetaStar_R = new TH2D("mtt_costhetastar_r", "m_{tt} cos#theta^{*} (reco)", nbins, Emin, Emax, 2, -1.0, 1.0);
    h2_mtt_cosThetaStar_R->GetXaxis()->SetTitle("m_{tt} (reco)");
    h2_mtt_cosThetaStar_R->GetYaxis()->SetTitle("cos#theta^* (reco)");
    h2_mtt_cosThetaStar_R->Sumw2();

    h_costheta_tt = new TH1D("costheta_tt", "cos#theta_{t,#bar{t}}", 20, -1.0, 1.0);
    h_costheta_tt->Sumw2();

    h_mtt_LL = new TH1D("mtt_LL", "m_{tt}^{LL}", nbins, Emin, Emax);
    h_mtt_LL->Sumw2();
    h_mtt_LR = new TH1D("mtt_LR", "m_{tt}^{LR}", nbins, Emin, Emax);
    h_mtt_LR->Sumw2();
    h_mtt_RL = new TH1D("mtt_RL", "m_{tt}^{RL}", nbins, Emin, Emax);
    h_mtt_RL->Sumw2();
    h_mtt_RR = new TH1D("mtt_RR", "m_{tt}^{RR}", nbins, Emin, Emax);
    h_mtt_RR->Sumw2();

    if (true) {

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

        h_HT = new TH1D("HT", "H_{T}", 40, 0, 4);
        h_HT->Sumw2();

        h_KT = new TH1D("KT", "K_{T}", 40, 0, 4);
        h_KT->Sumw2();

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

        h2_mtt_cosTheta1 = new TH2D("mtt_costheta_tl1", "m_{tt} cos#theta_{t,l+}", nbins, Emin, Emax, 20, -1.0, 1.0);
        h2_mtt_cosTheta1->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h2_mtt_cosTheta1->GetYaxis()->SetTitle("cos#theta_{l+}");
        h2_mtt_cosTheta1->Sumw2();

        h2_mtt_cosTheta2 = new TH2D("mtt_costheta_tl2", "m_{tt} cos#theta_{t,l-}", nbins, Emin, Emax, 20, -1.0, 1.0);
        h2_mtt_cosTheta2->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h2_mtt_cosTheta2->GetYaxis()->SetTitle("cos#theta_{l-}");
        h2_mtt_cosTheta2->Sumw2();

        h2_mtt_cosThetal_R = new TH2D("mtt_costhetal_R", "m_{tt} cos#theta_{t,l} (reco)", nbins, Emin, Emax, 20, -1.0, 1.0);
        h2_mtt_cosThetal_R->GetXaxis()->SetTitle("m_{tt} (reco) [TeV]");
        h2_mtt_cosThetal_R->GetYaxis()->SetTitle("cos#theta_{l}");
        h2_mtt_cosThetal_R->Sumw2();

        h2_mtt_cos1cos2 = new TH2D("mtt_cos1cos2", "m_{tt} cos#theta_{l+}cos#theta_{l-}", nbins, Emin, Emax, 20, -1.0, 1.0);
        h2_mtt_cos1cos2->GetXaxis()->SetTitle("m_{tt}");
        h2_mtt_cos1cos2->GetYaxis()->SetTitle("cos#theta_{l+}cos#theta_{l-}");
        h2_mtt_cos1cos2->Sumw2();

        h2_HT_deltaPhi = new TH2D("HT_deltaphi", "H_{T} #Delta#phi_{l}", 40, 0, 4, 10, 0, 1);
        h2_HT_deltaPhi->GetXaxis()->SetTitle("H_{T} [TeV]");
        h2_HT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
        h2_HT_deltaPhi->Sumw2();

        h2_KT_deltaPhi = new TH2D("KT_deltaphi", "K_{T} #Delta#phi_{l}", 40, 0, 4, 10, 0, 1);
        h2_KT_deltaPhi->GetXaxis()->SetTitle("K_{T} [TeV]");
        h2_KT_deltaPhi->GetYaxis()->SetTitle("#Delta#phi_{l} [rad] / #pi");
        h2_KT_deltaPhi->Sumw2();

        std::vector<string> deltaRnames, deltaRtitles;
        std::vector<string> n_eta, t_eta, n_pt, t_pt;

        std::vector<string> particles1 = {"b1", "b2", "l", "v", "q1", "q2"};
        std::vector<string> particles2 = {"b", "#bar{b}", "l+", "#nu", "q", "#bar{q}'"};

        for (int i = 0; i < 6; i++ ) {
            n_eta.push_back("eta_" + particles1[i]);
            t_eta.push_back("#eta_{" + particles2[i] + "}");
        }

        for (int i = 0; i < 6; i++ ) {
            n_pt.push_back("pT_" + particles1[i]);
            t_pt.push_back("p_{T}^{" + particles2[i] + "}");
        }

        for (int i = 0; i < 6; i++ ) {
            for (int j = i + 1; j < 6; j++) {
                deltaRnames.push_back("delta_R" + particles1[i] + particles1[j]);
                deltaRtitles.push_back("#Delta R(" + particles2[i] + "," + particles2[j] + ")");
            }
        }

        h_deltaR_tt = new TH1D("deltaR_tt", "#Delta R(t,#bar{t})", 100, 0, 5);
        h_deltaR_tt->Sumw2();

        h_deltaR_bW = new TH1D("deltaR_bW", "#Delta R(b,W)", 100, 0, 5);
        h_deltaR_bW->Sumw2();
        h_deltaR_max = new TH1D("deltaR_max", "#Delta R_{max}", 100, 0, 5);
        h_deltaR_max->Sumw2();

        for (int i = 0; i < (int) deltaRnames.size(); i++)
            h_deltaRs.push_back(new TH1D(deltaRnames[i].c_str(), deltaRtitles[i].c_str(), 100, 0, 5));

        for (int i = 0; i < (int) n_eta.size(); i++) {
            h_eta.push_back(new TH1D(n_eta[i].c_str(), t_eta[i].c_str(), 100, 0, 5));
            h_pt.push_back(new TH1D(n_pt[i].c_str(), t_pt[i].c_str(), 80, 0, 2000));
        }

        if (m_reco > 0) {
            h_mtt_R = new TH1D("m_tt_R", "m_{tt}^{reco}", nbins, Emin, Emax);
            h_mtt_R->Sumw2();

            h_pxt_R = new TH1D("px_t_R", "p_{x}_{t} (reco)", 100, 0, 5000);
            h_pxt_R->Sumw2();
            h_pyt_R = new TH1D("py_t_R", "p_{y}_{t} (reco)", 100, 0, 5000);
            h_pyt_R->Sumw2();
            h_pzt_R = new TH1D("pz_t_R", "p_{z}^{t} (reco)", 100, 0, 5000);
            h_pzt_R->Sumw2();
            h_Et_R = new TH1D("E_t_R", "E_{t} (reco)", 100, 0, 5000);
            h_Et_R->Sumw2();
            h_pTt_R = new TH1D("pT_t_R", "p_{T}^{t} (reco)", 100, 0, 5000);
            h_pTt_R->Sumw2();
            h_etat_R = new TH1D("eta_t_R", "#eta_{t} (reco)", nbins, 0, 10);
            h_etat_R->Sumw2();
            h_phit_R = new TH1D("phi_t_R", "#phi_{t} (reco)", nbins, -1, 1);
            h_phit_R->Sumw2();
            h_mt_R = new TH1D("m_t_R", "m_{t} (reco)", 40, 100, 300);
            h_mt_R->Sumw2();

            h_pxtbar_R = new TH1D("px_tbar_R", "p_{x}_{#bar{t}} (reco)", 100, 0, 5000);
            h_pxtbar_R->Sumw2();
            h_pytbar_R = new TH1D("py_tbar_R", "p_{y}_{#bar{t}} (reco)", 100, 0, 5000);
            h_pytbar_R->Sumw2();
            h_pztbar_R = new TH1D("pz_tbar_R", "p_{z}^{#bar{t}} (reco)", 100, 0, 5000);
            h_pztbar_R->Sumw2();
            h_Etbar_R = new TH1D("E_tbar_R", "E_{#bar{t}} (reco)", 100, 0, 5000);
            h_Etbar_R->Sumw2();
            h_pTtbar_R = new TH1D("pT_tbar_R", "p_{T}^{#bar{t}} (reco)", 100, 0, 5000);
            h_pTtbar_R->Sumw2();
            h_etatbar_R = new TH1D("eta_tbar_R", "#eta_{#bar{t}} (reco)", nbins, 0, 10);
            h_etatbar_R->Sumw2();
            h_phitbar_R = new TH1D("phi_tbar_R", "#phi_{#bar{t}} (reco)", nbins, -1, 1);
            h_phitbar_R->Sumw2();
            h_mtbar_R = new TH1D("m_tbar_R", "m_{#bar{t}} (reco)", 40, 100, 300);
            h_mtbar_R->Sumw2();

            h_costheta_tt_R = new TH1D("costheta_tt_R", "cos#theta_{t,#bar{t}} (reco)", 20, -1.0, 1.0);
            h_mtt_tF_R = new TH1D("mtt_FR", "m_{tt}^{forward} (reco)", nbins, Emin, Emax);
            h_mtt_tF_R->Sumw2();
            h_mtt_tB_R = new TH1D("mtt_BR", "m_{tt}^{backward} (reco)", nbins, Emin, Emax);
            h_mtt_tB_R->Sumw2();
            h_ytt_R = new TH1D("ytt_R", "y_{tt}^{reco}", 50, -2.5, 2.5);
            h_ytt_R->Sumw2();
            h_cosTheta_R = new TH1D("costheta_R", "cos#theta_{reco}", nbins, -1.0, 1.0);
            h_cosTheta_R->Sumw2();
            h_cosThetaStar_R = new TH1D("costhetastar_R", "cos#theta_{reco}^{*}", nbins, -1.0, 1.0);
            h_cosThetaStar_R->Sumw2();
            h_pv1x_R = new TH1D("pv1x_R", "p_{x}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
            h_pv1x_R->Sumw2();
            h_pv1y_R = new TH1D("pv1y_R", "p_{y}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
            h_pv1y_R->Sumw2();
            h_pv1z_R = new TH1D("pv1z_R", "p_{z}^{#nu_{1}} (reco", nbins, -500.0, 500.0);
            h_pv1z_R->Sumw2();
            h_pv2x_R = new TH1D("pv2x_R", "p_{x}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
            h_pv2x_R->Sumw2();
            h_pv2y_R = new TH1D("pv2y_R", "p_{y}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
            h_pv2y_R->Sumw2();
            h_pv2z_R = new TH1D("pv2z_R", "p_{z}^{#nu_{2}} (reco", nbins, -500.0, 500.0);
            h_pv2z_R->Sumw2();

            h_pT_t_perf = new TH1D("pT_t_perf", "p_{T}^{t} performance", 111, -10.05, 1.05);
            h_pT_t_perf->Sumw2();
            h_pT_tbar_perf = new TH1D("pT_tbar_perf", "p_{T}^{#bar{t}} performance", 111, -10.05, 1.05);
            h_pT_tbar_perf->Sumw2();
            h_costhetastar_perf = new TH1D("costhetastar_perf", "cos#theta^{*} performance", 101, -5.05, 5.05);
            h_costhetastar_perf->Sumw2();
            h_costhetatl_perf = new TH1D("costhetatl_perf", "cos#theta_{t,l} performance", 101, -5.05, 5.05);
            h_costhetatl_perf->Sumw2();
            h_costhetatt_perf = new TH1D("costheta_tt_perf", "cos#theta_{t#bar{t}} performance", 101, -5.05, 5.05);
            h_costhetatt_perf->Sumw2();
            // h_m_tt_perf = new TH1D("m_tt_perf", "m_{t#bar{t}} performance", 111, -10.05, 1.05);
            h_m_tt_perf = new TH1D("m_tt_perf", "m_{t#bar{t}} performance", 21, -1.05, 1.05);
            h_m_tt_perf->Sumw2();
            h2_m_tt_costheta_tt_perf = new TH2D("m_tt_costheta_tt_perf", "m_{t#bar{t}} cos#theta_{t#bar{t}} performance", 111, -10.05, 1.05, 101, -5.05, 5.05);
            h2_m_tt_costheta_tt_perf->Sumw2();
            h2_m_tt_pT_t_perf = new TH2D("m_tt_pT_t_perf", "m_{t#bar{t}} p_{T}^{#bar{t}} performance", 111, -10.05, 1.05, 111, -10.05, 1.05);
            h2_m_tt_pT_t_perf->Sumw2();
        }
    }
}


void Analysis::MakeDistributions()
{
    this->MakeDistribution1D(h_mtt, "TeV");
    this->MakeDistribution1D(h_ytt, "");

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

    if (n < 6) {
        this->MakeDistribution1D(h_mtt_LL, "TeV");
        this->MakeDistribution1D(h_mtt_LR, "TeV");
        this->MakeDistribution1D(h_mtt_RL, "TeV");
        this->MakeDistribution1D(h_mtt_RR, "TeV");

        h_ALL = this->MakeALL();
        h_ALL->GetYaxis()->SetTitle(h_ALL->GetTitle());
        h_ALL->GetXaxis()->SetTitle("m_{tt} [TeV]");

        h_AL = this->MakeAL();
        h_AL->GetYaxis()->SetTitle(h_AL->GetTitle());
        h_AL->GetXaxis()->SetTitle("m_{tt} [TeV]");
    }

    if (n == 6) {
 
        this->MakeDistribution1D(h_HT, "TeV");
        this->MakeDistribution1D(h_KT, "TeV");

        this->MakeDistribution1D(h_cosTheta1, "");
        this->MakeDistribution1D(h_cosTheta2, "");
        this->MakeDistribution1D(h_cos1cos2, "");
        this->MakeDistribution1D(h_deltaPhi, "rad / #pi");

        this->MakeDistribution1D(h_mtt_tlF, "TeV");
        this->MakeDistribution1D(h_mtt_tlB, "TeV");
        
        h_AtlFB = this->Asymmetry("AtlFB", "A^{tl}_{FB^*}", h_mtt_tlF, h_mtt_tlB);
        h_AtlFB->GetYaxis()->SetTitle(h_AtlFB->GetTitle());
        h_AtlFB->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h_AtlFB->Write();

        this->MakeDistribution1D(h_mtt_lF, "TeV");
        this->MakeDistribution1D(h_mtt_lB, "TeV");

        h_Ap = this->Asymmetry("Ap", "A_{P}", h_mtt_lF, h_mtt_lB);
        h_Ap->GetYaxis()->SetTitle(h_Ap->GetTitle());
        h_Ap->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h_Ap->Write();

        h_AlEl = this->Asymmetry("AlEl", "A^{l}_{E_l}", h_mtt_ElF, h_mtt_ElB);
        h_AlEl->GetYaxis()->SetTitle(h_AlEl->GetTitle());
        h_AlEl->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h_AlEl->Write();

        h_Aphil = this->Asymmetry("Aphil", "A^{l}_{#phi}", h_mtt_philF, h_mtt_philB);
        h_Aphil->GetYaxis()->SetTitle(h_Aphil->GetTitle());
        h_Aphil->GetXaxis()->SetTitle("m_{tt} [TeV]");
        h_Aphil->Write();

        this->MakeDistribution1D(h_pv1x, "GeV");
        this->MakeDistribution1D(h_pv1y, "GeV");
        this->MakeDistribution1D(h_pv1z, "GeV");
        this->MakeDistribution1D(h_pv2x, "GeV");
        this->MakeDistribution1D(h_pv2y, "GeV");
        this->MakeDistribution1D(h_pv2z, "GeV");

        for (auto& h : h_pt) {
            h->GetYaxis()->SetTitle("d#sigma / d p_{T}");
            h->GetXaxis()->SetTitle("p_{T}");
        }

        for (auto& h : h_eta) {
            h->GetYaxis()->SetTitle("d#sigma / d #eta");
            h->GetXaxis()->SetTitle("#eta");
        }

        for (auto& h_deltaR : h_deltaRs) {
            h_deltaR->GetYaxis()->SetTitle("d#sigma / d #Delta R");
            h_deltaR->GetXaxis()->SetTitle("#Delta R");
        }

        this->MakeDistribution1D(h_deltaR_bW, "");
        this->MakeDistribution1D(h_deltaR_max, "");
        this->MakeDistribution1D(h_deltaR_tt, "");

        if (m_reco > 0) {
            this->MakeDistribution1D(h_mtt_R, "TeV");
            this->MakeDistribution1D(h_mtt_tF_R, "TeV");
            this->MakeDistribution1D(h_mtt_tB_R, "TeV");
            this->MakeDistribution1D(h_ytt_R, "TeV");
            this->MakeDistribution1D(h_costheta_tt_R, "");

            this->MakeDistribution1D(h_pxt_R, "GeV");
            this->MakeDistribution1D(h_pyt_R, "GeV");
            this->MakeDistribution1D(h_pzt_R, "GeV");
            this->MakeDistribution1D(h_pTt_R, "GeV");
            this->MakeDistribution1D(h_etat_R, "");
            this->MakeDistribution1D(h_phit_R, "");
            this->MakeDistribution1D(h_Et_R, "GeV");
            this->MakeDistribution1D(h_mt_R, "GeV");

            this->MakeDistribution1D(h_pxtbar_R, "GeV");
            this->MakeDistribution1D(h_pytbar_R, "GeV");
            this->MakeDistribution1D(h_pztbar_R, "GeV");
            this->MakeDistribution1D(h_Etbar_R, "GeV");
            this->MakeDistribution1D(h_pTtbar_R, "GeV");
            this->MakeDistribution1D(h_etatbar_R, "");
            this->MakeDistribution1D(h_phitbar_R, "");
            this->MakeDistribution1D(h_mtbar_R, "GeV");

            this->MakeDistribution1D(h_cosTheta_R, "");
            this->MakeDistribution1D(h_cosThetaStar_R, "");

            this->MakeDistribution1D(h_pv1x_R, "GeV");
            this->MakeDistribution1D(h_pv1y_R, "GeV");
            this->MakeDistribution1D(h_pv1z_R, "GeV");
            this->MakeDistribution1D(h_pv2x_R, "GeV");
            this->MakeDistribution1D(h_pv2y_R, "GeV");
            this->MakeDistribution1D(h_pv2z_R, "GeV");

            this->MakeDistribution2D(h2_HT_deltaPhi, "H_{T}", "GeV", "#Delta#phi", "");
            this->MakeDistribution2D(h2_KT_deltaPhi, "K_{T}", "GeV", "#Delta#phi", "");

            this->MakeDistribution2D(h2_mtt_deltaPhi, "m_{tt}", "GeV", "#Delta#phi", "");
            this->MakeDistribution2D(h2_mtt_cos1cos2, "m_{tt}", "GeV", "cos#theta cos#theta", "");
            this->MakeDistribution2D(h2_mtt_cosThetaStar, "m_{tt}", "GeV", "cos#theta^{*}", "");
            this->MakeDistribution2D(h2_mtt_cosThetaStar_R, "m_{tt} (reco)", "GeV", "cos#theta^{*} (reco)", "");
            this->MakeDistribution2D(h2_mtt_cosTheta1, "m_{tt}", "GeV", "cos#theta_{l^{+}}", "");
            this->MakeDistribution2D(h2_mtt_cosTheta2, "m_{tt}", "GeV", "cos#theta_{l^{-}}", "");
            this->MakeDistribution2D(h2_mtt_cosThetal_R, "m_{tt}, (reco)", "GeV", "cos#theta_{l} (reco)", "");

            this->MakeDistribution1D(h_pT_t_perf, "");
            this->MakeDistribution1D(h_pT_tbar_perf, "");
            this->MakeDistribution1D(h_costhetatt_perf, "");
            this->MakeDistribution1D(h_costhetastar_perf, "");
            this->MakeDistribution1D(h_costhetatl_perf, "");
            this->MakeDistribution1D(h_m_tt_perf, "");
            this->MakeDistribution2D(h2_m_tt_pT_t_perf, "m_{tt} performance", "", "p_{T}^{t} performance", "");
            this->MakeDistribution2D(h2_m_tt_costheta_tt_perf, "m_{tt} performance", "", "cos#theta_{tt} performance", "");

            h_AtFB_R = this->Asymmetry("AtFB_R", "A^{t}_{FB^{*}} (reco)", h_mtt_tF_R, h_mtt_tB_R);
            h_AtFB_R->GetYaxis()->SetTitle(h_AtFB_R->GetTitle());
            h_AtFB_R->GetXaxis()->SetTitle("m_{tt} (reco) [TeV]");
            h_AtFB_R->Write();
        }
    }
}


void Analysis::MakeDistribution1D(TH1D* h, const TString& units)
{
    TString ytitle, yunits, xunits;
    if (m_xsec) {
        if (m_useLumi) {
            for (int i = 1; i < h->GetNbinsX() + 1; i++) {
                h->SetBinError(i, sqrt(h->GetBinContent(i)));
                std::cout << "N  = " << "" << h->GetBinContent(i) << std::endl;
                std::cout << "dN = " << "" << h->GetBinError(i) << std::endl;
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
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::MakeDistribution2D(TH2D* h, TString xtitle,  TString xunits, TString ytitle, TString yunits) {
    TString ztitle, zunits;
    if (m_xsec) {
        // h->Scale(m_sigma / m_nevents, "width");
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
    m_outputFile->cd();
    m_outputFile->cd("/");

    if (n < 6) {
        h_AL->Write();
        h_ALL->Write();
    }


    if (n == 6) {
        
        for (int i = 0; i < (int) h_deltaRs.size(); i++)
            h_deltaRs[i]->Write();

        for (int i = 0; i < (int) h_eta.size(); i++)
            h_eta[i]->Write();

        for (int i = 0; i < (int) h_pt.size(); i++)
            h_pt[i]->Write();


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

        TObjArray slicesR;
        this->NormalizeSliceY(h2_mtt_cosThetal_R);
        h2_mtt_cosThetal_R->FitSlicesY(func, 0, -1, 0, "QRN", &slicesR);
        for (auto slice : slicesR) slice->Write();
        TH1D* h_AL_R = (TH1D*) slicesR[0]->Clone("AL_R");
        slicesR.Clear();
        h_AL_R->Scale(2 / h2_mtt_cosThetal_R->GetYaxis()->GetBinWidth(1));
        h_AL_R->SetTitle("A_{L} (reco)");
        h_AL_R->GetXaxis()->SetTitle("m_{tt} (reco) [TeV]");
        h_AL_R->GetYaxis()->SetTitle("A_{L} (reco)");
        h_AL_R->Write();
    }
}


bool Analysis::PassFiducialCuts(const std::vector<TLorentzVector>& p, const TLorentzVector& P) 
{
    if (!m_fid) return true;
    if (this->PassCutsET(p, P)) {
        if (this->PassCutsEta(p, P)) {
            // if (this->PassCutsDeltaR(p, P)) {
                return true;
            // }
        }
    }
    return false;
}


bool Analysis::PassCuts(const std::vector<TLorentzVector>& p, const TLorentzVector& P) 
{
    // if (this->PassCutsET(p, P)) {
    //     if (this->PassCutsEta(p,P)) {
    //         if (this->PassCutsMET(p, P)) {
                if (this->PassCutsMtt(p, P)) {
                    return true;
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
}


void Analysis::SetDataDirectory()
{
    // sets directory based on hostname

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


void Analysis::GetGenerationCrossSection(TString filename)
{   
    std::cout << "filename = " << filename << std::endl;
    if (filename == "/Users/declan/Data/zprime/dd-AZ-tt-bbllvv.SM.13TeV.CT14LL.2-4.unweighted.root" || filename == "/Users/declan/Data/zprime/dd-AZX-tt-bbllvv.GLR-R-3.13TeV.CT14LL.2-4.unweighted.root") {
        m_sigma = 1000 * m_ptup->cross_section() / 100000; // temporary! Remove me
        std::cout << "bork" << std::endl;
    }
    else m_sigma = 1000 * m_ptup->cross_section();
    printf("Generation Cross section = %.15le [fb]\n", m_sigma);
}


void Analysis::Loop()
{
    for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
        printf("input %li: ", i - m_inputFiles->begin() + 1);
        std::cout << (*i) << std::endl;
        this->SetupTreesForNewFile((*i));
        this->GetGenerationCrossSection(*i);
        this->GetChannelFactors();
        m_nevents = this->TotalEvents();
        std::cout << "events: " << m_nevents << std::endl;
        for (Long64_t jevent = 0; jevent < m_nevents; ++jevent) {
            Long64_t ievent = this->IncrementEvent(jevent);
            if (ievent < 0) break;
            this->EachEvent();
            ProgressPercentage(jevent, m_nevents - 1, 50);
        }
        this->CleanUp();
    }
    std::cout << endl;
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
    TString treeToUse = "process";
    m_chainPtup = new TChain(treeToUse,"");
    TString TStringPtuple = s + "/" + treeToUse;
    m_chainPtup->Add(TStringPtuple,0);
    m_ptup = new process(m_chainPtup);
    m_ptup->LoadTree(0);

    treeToUse = "events";
    m_chainNtup = new TChain(treeToUse,"");
    TString TStringNtuple = s + "/" + treeToUse;
    m_chainNtup->Add(TStringNtuple,0);
    m_ntup = new events(m_chainNtup);
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
        printf("Error: Q_l must be ±1.\n");
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
    int nRealRoots = SolveP2(roots, b / a, c / a);
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
    if (nRealRoots == 2) {
        Root0MinusTruth = std::abs(roots[0] - pz_nu_truth);
        Root1MinusTruth = std::abs(roots[1] - pz_nu_truth);
    }
    else if (nRealRoots == 0) {
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
        printf("True pz_nu = %f\n", p_nu.Pz());
        printf("Possible neutrino solutions:\n");
        if (nRealRoots == 2) {
            printf("                             %f\n", roots[0]);
            printf("                             %f\n", roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        else if (nRealRoots == 0) {
            printf("                             %f + %fi\n", roots[0], roots[1]);
            printf("                             %f - %fi\n", roots[0], roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        if (imin == bestRoot) printf("Neutrino solution: correct. \n");
        else printf("Neutrino solution: incorrect. \n");
        if (b_lep == jmin) printf("b-assignment: correct. \n");
        else printf("b-assignment: incorrect. \n");
        printf("--\n");
    }
    return p_R;
}


std::vector<TLorentzVector> Analysis::ReconstructDilepton(const std::vector<TLorentzVector>& p)
{
    // Uses the Sonnenschein method algebraically solve tt dilepton equations. 
    // http://arxiv.org/abs/hep-ph/0510100
    // selects solution that minimises mtt

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
    // double mt1 = (p[0] + p[2] + p[3]).M(), mt2 = (p[1] + p[4] + p[5]).M();
    // double mw1 = (p[2] + p[3]).M(), mw2 = (p[4] + p[5]).M();
    // double mb1 = p[0].M(), mb2 = p[1].M();

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

    double c21 = 4 * (mw1 * mw1 - ml1 * ml1) * (pl1x - pl1z * a2 / a4)
               -8 * (El1 * El1 - pl1z*pl1z) * a1 * a2 / (a4 * a4) 
               -8 * pl1x * pl1z * a1 / a4;

    double c20 = -4 * (El1 * El1 - pl1x * pl1x)
               -4 * (El1 * El1 - pl1z * pl1z) * (a2 / a4) * (a2 / a4)
               -8 * pl1x * pl1z * a2 / a4;

    double c11 = 4 * (mw1 * mw1 - ml1 * ml1)*(pl1y - pl1z * a3 / a4)
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

    double x[4];
    int nRealRoots = SolveP4(x, a[0], a[1], a[2], a[3]);

    if (x[0] != x[0] && x[1] != x[1] && x[2] != x[2]) printf("error: Three NaNs in quartic solutions.\n");

    // if (x[0] != x[0] || x[1] != x[1] || x[2] != x[2] || x[3] != x[3]) printf("error: NaN in quartic solutions.\n");

    int nSolutions;
    std::vector<double> pv1x_Rs;
    if (nRealRoots == 4) {
        // they live in x[0], x[1], x[2], x[3].
        nSolutions = 4;
        for (int i = 0; i < nSolutions; i++) pv1x_Rs.push_back(x[i]);
    }
    else if (nRealRoots == 2) {
        // x[0], x[1] are the real roots and x[2]±i*x[3] are the complex
        nSolutions = 3;
        for (int i = 0; i < nSolutions; i++) pv1x_Rs.push_back(x[i]); 
    }
    else if (nRealRoots == 0) {
        // the equation has two pairs of pairs of complex conjugate roots in x[0]±i*x[1] and x[2]±i*x[3].
        nSolutions = 2;
        pv1x_Rs.push_back(x[0]);
        pv1x_Rs.push_back(x[2]); 
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

    return p_R;
}


void Analysis::GetChannelFactors()
{
    // scale dilepton to other classifications
    // fac_ee = 1
    // fac_emu = 2
    if (m_reco == 2) m_sigma = m_sigma * 24; // 2 [e+ + e-] x 2 [e + mu] x 6 [3 x 2]
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
    m_cutNames = std::vector<TString>(
    m_cuts,                       "no name         ");
    m_cutNames[c_entries]       = "entries         ";
    m_cutNames[c_topDecays]     = "t->bev          ";
    m_cutNames[c_antitopDecays] = "tbar->bev       ";
    m_cutNames[c_events]        = "events          ";
    m_cutNames[c_realSolutions] = "pvz real        ";
    m_cutNames[c_MET]           = "MET             ";
    m_cutNames[c_mtt]           = "mtt             ";
    m_cutNames[c_ytt]           = "ytt             ";
    m_cutNames[c_eta]           = "eta             ";
    m_cutNames[c_Et]            = "ET              ";
    m_cutNames[c_deltaR]        = "deltaR          ";

    h_cutflow = new TH1D("cutflow", "cutflow", m_cuts, 0, m_cuts);
}


void Analysis::PrintCutflow()
{
    std::cout << "cutflow: " << std::endl;
    for (int cut = 0; cut < m_cuts; cut++) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutNames[cut]);

        printf("%s %i pass\n", m_cutNames[cut].Data(), m_cutflow[cut]);
    }
    h_cutflow->Write();
    m_outputFile->Close();
    delete m_outputFile;
}

TH1D* Analysis::MakeALL () {
    TH1D* h_A = (TH1D*) h_mtt_LL->Clone();
    TH1D* h_B = (TH1D*) h_mtt_LR->Clone();

    h_A->Add(h_mtt_RR);
    h_B->Add(h_mtt_RL);

    double ALL = this->TotalAsymmetry(h_A, h_B);
    printf("ALL' = %f\n", ALL);

    TH1D* h_ALL = this->Asymmetry("ALL", "A_{LL}", h_A, h_B);
    delete h_A;
    delete h_B;
    return h_ALL;
}

TH1D* Analysis::MakeAL () {
    TH1D* h_A = (TH1D*) h_mtt_LL->Clone();
    TH1D* h_B = (TH1D*) h_mtt_RR->Clone();

    h_A->Add(h_mtt_LR);
    h_B->Add(h_mtt_RL);

    double AL = this->TotalAsymmetry(h_A, h_B);
    printf("AL' = %f\n", AL);

    TH1D* h_AL = this->Asymmetry("AL", "A_{L}", h_A, h_B);
    delete h_A;
    delete h_B;
    return h_AL;
}

void Analysis::TotalSpinAsymmetries () {
    double sigmaLL = h_mtt_LL->Integral();
    double sigmaLR = h_mtt_LR->Integral();
    double sigmaRL = h_mtt_RL->Integral();
    double sigmaRR = h_mtt_RR->Integral();

    double ALL = (sigmaLL + sigmaRR - sigmaRL - sigmaLR)/
                 (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

    double AL = (sigmaLL + sigmaLR - sigmaRL - sigmaRR)/
                (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

    printf("ALL = %f\n", ALL);
    printf("AL = %f\n", AL);
}

double Analysis::TotalAsymmetry(TH1D* h_A, TH1D* h_B) {
    double A = h_A->Integral("width");
    double B = h_B->Integral("width");
    double Atot = (A - B)/(A + B);
    return Atot;
}
