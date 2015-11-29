#include "analysis.hpp"

using namespace std;
using namespace boost;

string trim(string const& str)
{
    if(str.empty())
        return str;

    size_t firstScan = str.find_first_not_of(' ');
    size_t first     = firstScan == string::npos ? str.length() : firstScan;
    size_t last      = str.find_last_not_of(' ');
    return str.substr(first, last-first+1);
}


AnalysisZprime::AnalysisZprime(const TString channel, const TString model, const int energy, const TString options, const int vegasIterations, const int vegasPoints, const bool addQCD, const int luminosity, const int btags, const bool discardComplex, const TString analysisLabel):
  m_channel(channel),
  m_model(model),
  m_energy(energy),
  m_options(options),
  m_vegasIterations(vegasIterations),
  m_vegasPoints(vegasPoints),
  m_addQCD(addQCD),
  m_luminosity(luminosity),
  nBtags(btags),
  m_discardComplex(discardComplex),
  m_analysisLabel(analysisLabel),
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_Wmass(80.23),
  m_tmass(175.0),
  m_discardEvent(false),
  m_inputFiles(NULL),
  m_weightFiles(NULL),
  m_ntup(NULL),
  m_chainNtup(NULL),
  m_outputFile(NULL)
{
  this->PreLoop();
  this->Loop();
  this->PostLoop();
}

inline string BoolToString(bool b){return b ? "1" : "0";}

void AnalysisZprime::CreateFilenames(){
  TString base = m_dataDirectory + "/" + m_channel + "_" + m_model + "_" + to_string(m_energy) + m_options + to_string(m_vegasIterations) + "x" + to_string(m_vegasPoints);
  m_inputFileName = base + ".root";
  m_QCDfilename = m_dataDirectory + "/" + m_channel + "_QCD_" + to_string(m_energy) + m_options + to_string(m_vegasIterations) + "x" + to_string(m_vegasPoints);
  // m_inputFilenames.push_back(m_QCDfilename)
  // m_inputFilenames.push_back(base)
  m_QCDweightFile = m_QCDfilename + ".log";
  m_QCDfilename = m_QCDfilename + ".root";
  m_weightsFileName = base + ".log";
  TString QCDadded;
  if (m_addQCD) QCDadded = "QCD";
  else QCDadded = "";
  m_outputFileName = base + "." + QCDadded + to_string(nBtags) + BoolToString(m_discardComplex) + m_analysisLabel;
  if (m_luminosity > 0) m_outputFileName += "_" + to_string(m_luminosity);
  m_outputFileName += ".root";
  // printf("Input: '%s'.\n", m_inputFileName.Data());
  printf("Output: '%s'.\n", m_outputFileName.Data());
}

// bool exists(const TString& name) {
//   struct stat buffer;
//   bool exist(stat (name.Data(), &buffer) == 0);
//   printf("%s does not exist.", name);
//   return exist;
// }

TString AnalysisZprime::GetOutputFilename(){
  return m_outputFileName;
}

void AnalysisZprime::EachEvent () {
  m_discardEvent = false;
  UpdateCutflow(c_entries, true);
  p = vector<TLorentzVector>(6);
  for (unsigned int i = 0; i < m_ntup->E()->size(); i++) {
    p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));
  }

  P.SetPxPyPzE(0,0,0,0);
  pcm = vector<TLorentzVector>(p.size());
  for (unsigned int i = 0; i < p.size(); i++) {
    P += p[i];
    pcm[i] = p[i];
  }

  TVector3 V = -1*P.BoostVector();

  Pcm.SetPxPyPzE(0,0,0,0);
  for (unsigned int i = 0; i < p.size(); i++) {
    pcm[i].Boost(V);
    Pcm += pcm[i];
  }

  if (m_channel == "bbllnn") {
    p_R1 = this->ReconstructSemiLeptonic(p,1); // top decays leptonically
    p_R2 = this->ReconstructSemiLeptonic(p,-1); // top decays hadronically

    if (m_discardEvent) return;

    P_R1.SetPxPyPzE(0,0,0,0);
    P_R2.SetPxPyPzE(0,0,0,0);
    for (unsigned int i = 0; i < p.size(); i++) {
      P_R1 += p_R1[i];
      P_R2 += p_R2[i];
    }

    // reconstructed final particle parton CoM variables
    TVector3 V_R1 = -1*P_R1.BoostVector();
    TVector3 V_R2 = -1*P_R2.BoostVector();
    pcm_R1 = vector<TLorentzVector>(p.size());
    pcm_R2 = vector<TLorentzVector>(p.size());
    for (unsigned int i = 0; i < p.size(); i++) {
      pcm_R1[i] = p_R1[i];
      pcm_R2[i] = p_R2[i];
      pcm_R1[i].Boost(V_R1);
      pcm_R2[i].Boost(V_R2);
    }
  }

  // top and antitop
  if (m_channel == "tt" or m_channel == "ll") {
    p_t = pcm[0];
    p_tb = pcm[1];
  }
  else if (m_channel == "bbllnn") {
    p_t = pcm[0] + pcm[2] + pcm[3];
    p_tb = pcm[1] + pcm[4] + pcm[5];
    p_t_R1 = pcm_R1[0] + pcm_R1[2] + pcm_R1[3];
    p_tb_R1 = pcm_R1[1] + pcm_R1[4] + pcm_R1[5];
    p_t_R2 = pcm_R2[0] + pcm_R2[2] + pcm_R2[3];
    p_tb_R2 = pcm_R2[1] + pcm_R2[4] + pcm_R2[5];
  }

  double mtt = P.M()/1000;
  double mtt_R1 = P_R1.M()/1000;
  double mtt_R2 = P_R2.M()/1000;
  double mt = p_t.M();
  double mtb = p_tb.M();
  double y_t = p_t.Rapidity();
  double y_tb = p_tb.Rapidity();
  double dy = abs(y_t) - abs(y_tb);
  double ytt = P.Rapidity();
  double cosTheta = p_t.CosTheta();
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

  if (m_channel == "bbllnn"){
    ytt_R1 = P_R1.Rapidity();
    ytt_R2 = P_R2.Rapidity();
    mt_R1 = p_t_R1.M();
    mtb_R1 = p_tb_R1.M();
    mt_R2 = p_t_R2.M();
    mtb_R2 = p_tb_R2.M();

    // printf("Reconstructed top mass\n---\n");
    // printf("m_top = %f TeV\n", mt);
    // printf("m_antitop = %f TeV\n", mtb);
    // printf("m_top (reco)[lep] = %f TeV\n", mt_R1);
    // printf("m_antitop (reco)[had] = %f TeV\n", mtb_R1);
    // printf("m_top (reco)[had] = %f TeV\n", mt_R2);
    // printf("m_antitop (reco)[lep] = %f TeV\n", mtb_R2);
    // printf("---\n");

    cosTheta_R1 = p_t_R1.CosTheta();
    cosTheta_R2 = p_t_R2.CosTheta();
    cosThetaStar_R1 = int(ytt_R1/abs(ytt_R1))*cosTheta_R1;
    cosThetaStar_R2 = int(ytt_R2/abs(ytt_R2))*cosTheta_R2;
  }

  if (this->PassCuts())
  {

    // re-weight for different iterations
    double it = m_ntup->iteration();
    double weight = m_ntup->weight();
    weight = weight*m_sigma/iteration_weights[it-1];
    double weight_R = weight/2;

    // fill histograms (assumes fixed bin width!)
    h_mt->Fill(mt, weight/h_mt->GetXaxis()->GetBinWidth(1));
    h_mtbar->Fill(mtb, weight/h_mtbar->GetXaxis()->GetBinWidth(1));
    h_mtt->Fill(mtt, weight/h_mtt->GetXaxis()->GetBinWidth(1));
    h_ytt->Fill(ytt, weight/h_ytt->GetXaxis()->GetBinWidth(1));
    h_cosTheta->Fill(cosTheta, weight/h_cosTheta->GetXaxis()->GetBinWidth(1));
    h_cosThetaStar->Fill(cosThetaStar, weight/h_cosThetaStar->GetXaxis()->GetBinWidth(1));

    // asymmetries
    if (cosThetaStar > 0) h_mtt_F->Fill(mtt, weight/h_mtt_F->GetXaxis()->GetBinWidth(1));
    if (cosThetaStar < 0) h_mtt_B->Fill(mtt, weight/h_mtt_B->GetXaxis()->GetBinWidth(1));

    if (dy > 0) h_mtt_Fy->Fill(mtt, weight/h_mtt_Fy->GetXaxis()->GetBinWidth(1));
    if (dy < 0) h_mtt_By->Fill(mtt, weight/h_mtt_By->GetXaxis()->GetBinWidth(1));

    if (m_channel == "tt") {
      h_mtt_LL->Fill(mtt, m_ntup->weightLL()/h_mtt_LL->GetXaxis()->GetBinWidth(1));
      h_mtt_LR->Fill(mtt, m_ntup->weightLR()/h_mtt_LR->GetXaxis()->GetBinWidth(1));
      h_mtt_RL->Fill(mtt, m_ntup->weightRL()/h_mtt_RL->GetXaxis()->GetBinWidth(1));
      h_mtt_RR->Fill(mtt, m_ntup->weightRR()/h_mtt_RR->GetXaxis()->GetBinWidth(1));
    }
    else if (m_channel == "bbllnn") {
      h_mtt_R->Fill(mtt_R1, weight_R/h_mtt_R->GetXaxis()->GetBinWidth(1));
      h_mtt_R->Fill(mtt_R2, weight_R/h_mtt_R->GetXaxis()->GetBinWidth(1));

      if (cosThetaStar_R1 > 0) h_mtt_FR->Fill(mtt_R1, weight_R/h_mtt_FR->GetXaxis()->GetBinWidth(1));
      if (cosThetaStar_R1 < 0) h_mtt_BR->Fill(mtt_R1, weight_R/h_mtt_BR->GetXaxis()->GetBinWidth(1));

      if (cosThetaStar_R2 > 0) h_mtt_FR->Fill(mtt_R2, weight_R/h_mtt_FR->GetXaxis()->GetBinWidth(1));
      if (cosThetaStar_R2 < 0) h_mtt_BR->Fill(mtt_R2, weight_R/h_mtt_BR->GetXaxis()->GetBinWidth(1));

      h_ytt_R->Fill(ytt_R1, weight_R/h_ytt_R->GetXaxis()->GetBinWidth(1));
      h_ytt_R->Fill(ytt_R2, weight_R/h_ytt_R->GetXaxis()->GetBinWidth(1));

      h_mt_R->Fill(mt_R1, weight_R/h_mt_R->GetXaxis()->GetBinWidth(1));
      h_mt_R->Fill(mt_R2, weight_R/h_mt_R->GetXaxis()->GetBinWidth(1));

      h_mtbar_R->Fill(mtb_R1, weight_R/h_mtbar_R->GetXaxis()->GetBinWidth(1));
      h_mtbar_R->Fill(mtb_R2, weight_R/h_mtbar_R->GetXaxis()->GetBinWidth(1));

      h_cosTheta_R->Fill(cosTheta_R1, weight_R/h_cosTheta_R->GetXaxis()->GetBinWidth(1));
      h_cosTheta_R->Fill(cosTheta_R2, weight_R/h_cosTheta_R->GetXaxis()->GetBinWidth(1));

      h_cosThetaStar_R->Fill(cosThetaStar_R1, weight_R/h_cosThetaStar_R->GetXaxis()->GetBinWidth(1));
      h_cosThetaStar_R->Fill(cosThetaStar_R2, weight_R/h_cosThetaStar_R->GetXaxis()->GetBinWidth(1));

      h_pzNu->Fill(p[3].Pz(), weight_R/h_pzNu->GetXaxis()->GetBinWidth(1));

      h_pzNu_R->Fill(p_R1[3].Pz(), weight_R/h_pzNu_R->GetXaxis()->GetBinWidth(1));
      h_pzNu_R->Fill(p_R2[5].Pz(), weight_R/h_pzNu_R->GetXaxis()->GetBinWidth(1));
    }
  }
}

void AnalysisZprime::GetCrossSection(){
  TString log(m_weightFiles->at(m_ifile));
  ifstream logstream(log.Data());
  if (!logstream.is_open()) printf("Error: failed to open %s!\n", log.Data());
  string line;
  string target = "Cross section";
  vector<string> parts;
  bool found = false;
  while(getline(logstream, line)){
    trim(line);
    split(parts, line, is_any_of(":"));
    for(auto part: parts){
      trim(part);
    }
    // printf("Part 0: %s\n", parts[0].c_str());
    if (parts[0] == target){
      m_sigma = stod(parts[1]);
      found = true;
    }
  }
  logstream.close();
  if (!found) {
    printf("Error: Failed to read generation cross section. Check target log file: %s", m_weightsFileName.Data());
    exit(1);
  }
  else printf("Generation Cross section = %.15le [pb]\n", m_sigma);
}

void AnalysisZprime::PostLoop () {
  this->CheckResults();
  if (m_channel == "tt") this->TotalSpinAsymmetries();
  if (m_channel == "bbllnn") this->CheckPerformance();
  this->MakeGraphs();
  this->PrintCutflow();
  this->WriteHistograms();
  printf("--- Done ---\n");
}


void AnalysisZprime::CheckResults() {
  // printf("--- Results ---\n");
  double sigma = h_mtt->Integral("width");
  if ((abs(sigma - m_sigma) > 10e-11) and (m_cuts < 2)) {
    printf("Cross section from generation and analysis stage do not match!\n");
    printf("sigma_generation = %.15le\n", m_sigma);
    printf("sigma_analysis   = %.15le\n", sigma);
  }
  else printf("Analysis Cross Section = %.15le [pb]\n", sigma);
}


void AnalysisZprime::CheckPerformance () {
  // printf("m_nQuarksMatched = %i\n", m_nQuarksMatched);
  // printf("m_nNeutrinoMatched = %i\n", m_nNeutrinoMatched);
  // printf("m_nReco = %i\n", m_nReco);

  double quarkRecoRatio = m_nQuarksMatched/(double)m_nReco;
  double neutrinoRecoRatio = m_nNeutrinoMatched/(double)m_nReco;
  printf("--- Performance ---\n");
  printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio*100);
  printf("pzNu assignment: %.1f%% correct\n", neutrinoRecoRatio*100);
}


TH1D* AnalysisZprime::MakeALL () {
  TH1D* h_A = (TH1D*) h_mtt_LL->Clone();
  TH1D* h_B = (TH1D*) h_mtt_LR->Clone();

  h_A->Add(h_mtt_RR);
  h_B->Add(h_mtt_RL);

  double ALL = this->TotalAsymmetry(h_A,h_B);
  printf("ALL' = %f\n", ALL);

  TH1D* h_ALL = this->Asymmetry("ALL", "A_{LL}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_ALL;
}


TH1D* AnalysisZprime::MakeAL () {
  TH1D* h_A = (TH1D*) h_mtt_LL->Clone();
  TH1D* h_B = (TH1D*) h_mtt_RR->Clone();

  h_A->Add(h_mtt_LR);
  h_B->Add(h_mtt_RL);

  TH1D* h_AL = this->Asymmetry("AL", "A_{L}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_AL;
}


void AnalysisZprime::TotalSpinAsymmetries () {
  double sigmaLL = h_mtt_LL->Integral("width");
  double sigmaLR = h_mtt_LR->Integral("width");
  double sigmaRL = h_mtt_RL->Integral("width");
  double sigmaRR = h_mtt_RR->Integral("width");

  double ALL = (sigmaLL + sigmaRR - sigmaRL - sigmaLR)/
               (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

  double AL = (sigmaRR + sigmaRL - sigmaLR - sigmaLL)/
              (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

  printf("ALL = %f\n", ALL);
  printf("AL = %f\n", AL);
}


double AnalysisZprime::TotalAsymmetry(TH1D* h_A, TH1D* h_B) {
  double A = h_A->Integral("width");
  double B = h_B->Integral("width");
  double Atot = (A - B)/(A + B);
  return Atot;
}


void AnalysisZprime::ApplyLuminosity(TH1D* h) {
  if (!m_useLumi) return;
  // printf("Name: %s\n", h->GetTitle());
  // printf("Luminosity: %f\n", m_luminosity);
  double sigma = -999, N = -999, dN = -999;
  double efficiency = 1.0, pb = 1000;
  for (int i = 1; i < h->GetNbinsX(); i++) {
    sigma = h->GetBinContent(i)*h->GetBinWidth(i);
    N = m_luminosity*pb*efficiency*sigma;
    h->SetBinContent(i, N);
    dN = sqrt(N);
    h->SetBinError(i, dN);
    // printf("sigma = %f, N = %f, dN = %f\n", sigma, N, dN);
  }
  TString events = "Events";
  h->GetYaxis()->SetTitle(events);
}


TH1D* AnalysisZprime::Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B) {
  TH1D* h_numerator = (TH1D*) h_A->Clone(name);
  TH1D* h_denominator = (TH1D*) h_A->Clone();
  h_numerator->SetTitle(title);
  h_numerator->Add(h_B, -1);
  h_denominator->Add(h_B, 1);
  h_numerator->Divide(h_denominator);
  delete h_denominator;
  if (m_luminosity > -1) this->AsymmetryUncertainty(h_numerator, h_A, h_B);
  return h_numerator;
}



void AnalysisZprime::AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B) {
  double A, deltaA, N, N_A, N_B;
  for (int i = 1; i < h_Asymmetry->GetNbinsX(); i++) {
    A = h_Asymmetry->GetBinContent(i);
    N_A = h_A->GetBinContent(i);
    N_B = h_B->GetBinContent(i);
    N = N_A + N_B;
    if (N > 0) deltaA = sqrt((1.0 - A*A)/N);
    else deltaA = 0;
    // printf("A = %f, dA= %f, N= %f\n", A, deltaA, N);
    h_Asymmetry->SetBinError(i, deltaA);
  }
}


void AnalysisZprime::CreateHistograms() {

  double binWidth = 0.05;
  double Emin = 0.025;
  double Emax = 13.025;
  double nbins = (Emax-Emin)/binWidth;

  if (m_channel == "ll") {
    h_mtt = new TH1D("mll", "m_{ll}", nbins, Emin, Emax);
    h_mtt->Sumw2();
  }

  if (m_channel == "tt" or m_channel == "bbllnn") {
    h_mtt = new TH1D("mtt", "m_{tt}", nbins, Emin, Emax);
    h_mtt->Sumw2();
    h_ytt = new TH1D("ytt", "y_{tt}", nbins, -2.5, 2.5);
    h_ytt->Sumw2();
    h_mt = new TH1D("mt", "m_{t}", nbins, 0, 350);
    h_mt->Sumw2();
    h_mtbar = new TH1D("mtbar", "m_{#bar{t}}", nbins, 0, 350);
    h_mtbar->Sumw2();
    h_mtt_F = new TH1D("mtt_F", "m_{tt}^{forward}", nbins, Emin, Emax);
    h_mtt_F->Sumw2();
    h_mtt_B = new TH1D("mtt_B", "m_{tt}^{backward}", nbins, Emin, Emax);
    h_mtt_B->Sumw2();
    h_mtt_Fy = new TH1D("mtt_Fy", "m_{tt}^{F(y)}", nbins, Emin, Emax);
    h_mtt_Fy->Sumw2();
    h_mtt_By = new TH1D("mtt_By", "m_{tt}^{B(y)}", nbins, Emin, Emax);
    h_mtt_By->Sumw2();
    h_cosTheta = new TH1D("cosTheta", "cos#theta", nbins, -1.0, 1.0);
    h_cosTheta->Sumw2();
    h_cosThetaStar = new TH1D("cosThetaStar", "cos#theta^{*}", nbins, -1.0, 1.0);
    h_cosThetaStar->Sumw2();
  }

  if (m_channel == "tt") {
    h_mtt_LL = new TH1D("mtt_LL", "m_{tt}^{LL}", nbins, Emin, Emax);
    h_mtt_LL->Sumw2();
    h_mtt_LR = new TH1D("mtt_LR", "m_{tt}^{LR}", nbins, Emin, Emax);
    h_mtt_LR->Sumw2();
    h_mtt_RL = new TH1D("mtt_RL", "m_{tt}^{RL}", nbins, Emin, Emax);
    h_mtt_RL->Sumw2();
    h_mtt_RR = new TH1D("mtt_RR", "m_{tt}^{RR}", nbins, Emin, Emax);
    h_mtt_RR->Sumw2();
  }

  if (m_channel == "bbllnn") {
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
    h_mtt_BD = new TH1D("mtt_BD", "m_{tt}^{backward} (reco)", nbins, Emin, Emax);
    h_mtt_BD->Sumw2();

    h_mtt_Fl = new TH1D("mtt_Fl", "m_{tt}^{F,l}", nbins, Emin, Emax);
    h_mtt_Fl->Sumw2();
    h_mtt_Bl = new TH1D("mtt_Bl", "m_{tt}^{B,l}", nbins, Emin, Emax);
    h_mtt_Bl->Sumw2();

    h_ytt_R = new TH1D("ytt_R", "y_{tt}^{reco}", nbins, -2.5, 2.5);
    h_ytt_R->Sumw2();
    h_cosTheta_R = new TH1D("cosTheta_R", "cos#theta_{reco}", nbins, -1.0, 1.0);
    h_cosTheta_R->Sumw2();
    h_cosThetaStar_R = new TH1D("cosThetaStar_R", "cos#theta_{reco}^{*}", nbins, -1.0, 1.0);
    h_cosThetaStar_R->Sumw2();
    h_pzNu = new TH1D("pzNu", "p_{z}^{#nu}", nbins, -500.0, 500.0);
    h_pzNu->Sumw2();
    h_pzNu_R = new TH1D("pzNu_R", "p_{z}^{#nu} (reco)", nbins, -500.0, 500.0);
    h_pzNu_R->Sumw2();
  }
}


void AnalysisZprime::MakeGraphs() {
  printf("Making Graphs...\n");

  this->MakeDistribution(h_mtt, "TeV");
  this->MakeDistribution(h_mtt_F, "TeV");
  this->MakeDistribution(h_mtt_B, "TeV");
  this->MakeDistribution(h_mt, "TeV");
  this->MakeDistribution(h_mtbar, "TeV");
  this->MakeDistribution(h_ytt, "TeV");
  this->MakeDistribution(h_cosTheta, "");
  this->MakeDistribution(h_cosThetaStar, "");

  h_mtt_Fn = (TH1D*) h_mtt_F->Clone("h_mtt_Fn");
  h_mtt_Fn->Divide(h_mtt);

  h_mtt_Bn = (TH1D*) h_mtt_B->Clone("h_mtt_Bn");
  h_mtt_Bn->Divide(h_mtt);

  h_AFB = this->Asymmetry("AFB", "A^{*}_{FB}", h_mtt_F, h_mtt_B);
  h_AFB->GetYaxis()->SetTitle(h_AFB->GetTitle());
  h_AFB->GetXaxis()->SetTitle("m_{tt} [TeV]");

  h_AC = this->Asymmetry("AC", "A_{C}", h_mtt_Fy, h_mtt_By);
  h_AC->GetYaxis()->SetTitle(h_AC->GetTitle());
  h_AC->GetXaxis()->SetTitle("m_{tt} [TeV]");

  if (m_channel == "tt") {
    this->MakeDistribution(h_mtt_LL, "TeV");
    this->MakeDistribution(h_mtt_LR, "TeV");
    this->MakeDistribution(h_mtt_RL, "TeV");
    this->MakeDistribution(h_mtt_RR, "TeV");

    h_ALL = this->MakeALL();
    h_ALL->GetYaxis()->SetTitle(h_ALL->GetTitle());
    h_ALL->GetXaxis()->SetTitle("m_{tt} [TeV]");

    h_AL = this->MakeAL();
    h_AL->GetYaxis()->SetTitle(h_AL->GetTitle());
    h_AL->GetXaxis()->SetTitle("m_{tt} [TeV]");
  }

  if (m_channel == "bbllnn") {
    this->MakeDistribution(h_mtt_R, "TeV");
    this->MakeDistribution(h_mtt_FR, "TeV");
    this->MakeDistribution(h_mtt_BR, "TeV");
    this->MakeDistribution(h_mt_R, "TeV");
    this->MakeDistribution(h_mtbar_R, "TeV");
    this->MakeDistribution(h_ytt_R, "TeV");
    this->MakeDistribution(h_cosTheta_R, "");
    this->MakeDistribution(h_cosThetaStar_R, "");
    this->MakeDistribution(h_pzNu, "GeV");
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
  if (m_channel == "tt") ytitle = "d#sigma(pp->t#bar{t}) / d";
  else ytitle = "d#sigma / d"; //pp->t#bar{t}->b#bar{b}l^{+}l^{-}#nu#bar{#nu}
  if (units != ""){
    yunits = " [pb/" + units + "]";
    xunits = " [" + units + "]";
  }
  else {
    yunits = "";
    xunits = "";
  }
  h->GetYaxis()->SetTitle(ytitle + h->GetTitle() + yunits);
  h->GetXaxis()->SetTitle(h->GetTitle() + xunits);
  this->ApplyLuminosity(h);
  m_outputFile->cd();
  m_outputFile->cd("/");
  h->Write();
}


void AnalysisZprime::WriteHistograms() {
  m_outputFile->cd();
  m_outputFile->cd("/");

  if (m_channel == "ll" or m_channel == "tt" or m_channel == "bbllnn") {
    h_mtt_Fn->Write();
    h_mtt_Bn->Write();
    h_AFB->Write();
    h_AC->Write();
  }

  if (m_channel == "ll" or m_channel == "tt") {
    h_mtt_LL->Write();
    h_mtt_LR->Write();
    h_mtt_RL->Write();
    h_mtt_RR->Write();
    h_ALL->Write();
    h_AL->Write();
  }

  if (m_channel == "bbllnn") {
    h_mtt_FRn->Write();
    h_mtt_BRn->Write();
    h_AFB_R->Write();
  }
  m_outputFile->Close();
  // delete h_cutflow;
  delete m_outputFile;

}


bool AnalysisZprime::PassCuts() {
  // if (this->PassCuts_mtt())
  // {
  //   if (this->PassCuts_MET())
  //     {
  //       if (this->PassCutsFiducial())
  //       {
        return true;
  //       }
  //     }
  // }
  // return false;
}


bool AnalysisZprime::PassCutsMET () {
  bool pass;
  if (m_channel == "ll") pass = true;
  else if (m_channel == "tt") pass = true;
  else if (m_channel == "bbllnn") pass = true;
  else pass = true;

  this->UpdateCutflow(c_MET, pass);
  return pass;
}


bool AnalysisZprime::PassCutsMtt () {
  // if (mtt > 1200.0)
  //   {
  //     if (mtt < 2200.0)
  //       {
          UpdateCutflow(c_mtt, true);
          return true;
  //       }
  //   }
  // UpdateCutflow(c_mtt, false);
  // return false;
}


bool AnalysisZprime::PassCutsFiducial () {
  for (unsigned int i = 0; i < p.size(); i++) {
    bool outsideCrack = p[i].PseudoRapidity() <= 1.37 || p[i].PseudoRapidity() >= 1.52;
    bool central      = p[i].PseudoRapidity() <= 2.47;
    bool passesFiducialCuts = outsideCrack && central;
    if (passesFiducialCuts == false){
      UpdateCutflow(c_fiducial, false);
      return false;
    }
    else continue;
  }
  UpdateCutflow(c_fiducial, true);
  return true;
}


bool AnalysisZprime::PassCutsYtt () {
  if (abs(P.Rapidity()) > 0)
  {
    UpdateCutflow(c_ytt, true);
    return true;
  }
  UpdateCutflow(c_mtt, false);
  return false;
}


void AnalysisZprime::PreLoop () {
  printf("--- Setup ---\n");
  this->GetDataDirectory();
  this->CreateFilenames();
  this->SetupInputFiles();
  this->SetupOutputFiles();
  this->ResetCounters();
  this->GetCrossSection();
  this->GetIterationWeights();
  this->GetChannelFactors();
  // this->CheckFiles();
  this->InitialiseCutflow();
  this->CreateHistograms();
}

void AnalysisZprime::GetDataDirectory(){
  #if __linux || __linux__
    m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime";
  #elif __APPLE__ || __MACH__
    m_dataDirectory = "/Users/declan/Data/Zprime";
  #endif
}


void AnalysisZprime::ResetCounters () {
  if (m_luminosity >= 0) m_useLumi = true;
  else m_useLumi = false;
  m_nQuarksMatched = 0;
  m_nNeutrinoMatched = 0;
  m_nReco = 0;
  m_nRealRoots = 0;
  m_nComplexRoots = 0;
}


void AnalysisZprime::GetIterationWeights() {
  iteration_weights.clear();
  TString log(m_weightFiles->at(m_ifile));
  ifstream logstream(log.Data());
  if (!logstream.is_open()) printf("Error: failed to open %s!\n", log.Data());

  string line;
  string target = "Iteration weighting";
  vector<string> parts;
  bool found = false;
  while(getline(logstream, line)){
    trim(line);
    split(parts, line, is_any_of(":"));
    for(auto part: parts){
      trim(part);
    }
    // printf("Part 0: %s\n", parts[0].c_str());
    if (parts[0] == target){
      // printf("parts[2] = %s\n", parts[2].c_str());
      iteration_weights.push_back(stod(parts[2]));
      found = true;
    }
  }
  logstream.close();
  if (!found) {
    printf("Error: Failed to read vegas iteration weights. Check target log file: %s", m_weightsFileName.Data());
    exit(1);
  }
  // for (auto iteration_weight: iteration_weights) {
  //   printf("VEGAS iteration weight = %f\n", iteration_weight);
  // }
}


void AnalysisZprime::Loop () {
  // Loop over all files
  for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i) {
    cout << "Processing:  '" << (*i) << "'." << endl;
    this->SetupTreesForNewFile((*i));
    Long64_t nEntries;
    nEntries = this->TotalEvents();
    printf("--- Event Loop ---\n");
    printf("File contains %lld entries.\n", nEntries);
    for (Long64_t jentry = 0; jentry < nEntries; ++jentry)
    {
      Long64_t ientry = this->IncrementEvent(jentry);
      if (ientry < 0) break;
      this->EachEvent();
      this->ProgressBar(jentry, nEntries-1, 50);
    }
    this->CleanUp();
    m_ifile++;
  }
}


AnalysisZprime::~AnalysisZprime() { delete m_inputFiles; }


void AnalysisZprime::SetupOutputFiles() {
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
}


void AnalysisZprime::SetupInputFiles () {
  m_inputFiles = new vector<TString>;
  m_weightFiles = new vector<TString>;
  if (m_addQCD) {
    m_inputFiles->push_back(m_QCDfilename);
    m_weightFiles->push_back(m_QCDweightFile);
  }
  m_inputFiles->push_back(m_inputFileName);
  // printf("m_weightsFileName = %s\n",m_weightsFileName.Data());
  // printf("size = %i\n", m_weightFiles->size());
  m_weightFiles->push_back(m_weightsFileName);

}


Long64_t AnalysisZprime::TotalEvents () {
  if (m_ntup != 0){return m_ntup->totalEvents();}
  return -999;
}


Long64_t AnalysisZprime::IncrementEvent(Long64_t i) {
  Long64_t ev(-1);
  if (m_ntup != 0){ev = m_ntup->LoadTree(i);}
  return ev;
}


void AnalysisZprime::SetupTreesForNewFile(const TString& s) {
  TString treeToUse = "RootTuple";

  m_chainNtup = new TChain(treeToUse,"");
  TString TStringNtuple = s + "/" + treeToUse;
  m_chainNtup->Add(TStringNtuple,0);
  m_ntup = new RootTuple(m_chainNtup);
}


void AnalysisZprime::CleanUp () {
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

  // printf("---\n");
  this->UpdateCutflow(c_events, true);
  m_nReco++;
  int nNoTag = 4 - nBtags;

  vector<TLorentzVector> p_R(p.size()); // I am returned!
  TLorentzVector p_l, p_nu;
  vector<TLorentzVector> p_b(nBtags), p_q(nNoTag);

  p_b[0] = p[0];
  p_b[1] = p[1];
  if (Q_l == 1) {
    this->UpdateCutflow(c_topDecays, true);
    p_l = p[2];
    p_nu = p[3];
    p_q[0] = p[4];
    p_q[1]= p[5];
  }
  else if (Q_l == -1) {
    this->UpdateCutflow(c_antitopDecays, true);
    p_l = p[4];
    p_nu = p[5];
    p_q[0] = p[2];
    p_q[1]= p[3];
  }
  else {
    printf("ERROR: Invalid lepton charge.\n");
  }

  // Calculate neutrino pz solutions
  double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
  double px_nu = p_nu.Px(), py_nu = p_nu.Py();
  vector<complex<double> > root;
  double a = -999, b = -999, c = -999, k = -999;

  E_l = sqrt(px_l*px_l + py_l*py_l + pz_l*pz_l);
  if (abs(E_l - p_l.E()) > 0.00001) printf("ERROR: Lepton energy doesn't match.\n");

  k = m_Wmass*m_Wmass/2 + px_l*px_nu + py_l*py_nu;
  a = px_l*px_l + py_l*py_l;
  b = -2*k*(pz_l);
  c = (px_nu*px_nu + py_nu*py_nu)*E_l*E_l - k*k;

  root = this->SolveQuadratic(a, b, c);

  // select single solution and match 'jets'
  double X2 = -999, X2min = -999;
  TLorentzVector p_nu_R;
  double dh = -999, dl = -999, E_nu_R = -999, mblv = -999, mjjb = -999;
  int imin = -999, jmin = -999, it = 0;
  vector<double> rootR(root.size());
  unsigned int nReal;

  // re-weight for different iterations NOT SURE WHY IM DOING THIS HERE?
  // double iteration = m_ntup->iteration();
  // double weight = m_ntup->weight();
  // weight = weight*m_sigma/iteration_weights[iteration-1];

  if (root[0].imag() == 0 and root[1].imag() == 0) {
    // two real solutions; pick best match
    nReal = 2;
    m_nRealRoots++;
    this->UpdateCutflow(c_realSolutions, true);
    if (Q_l == +1) m_R1solutionIsReal = true;
    if (Q_l == -1) m_R2solutionIsReal = true;
  }
  else {
    if (m_discardComplex) m_discardEvent = true;
    if (Q_l == +1) m_R1solutionIsReal = false;
    if (Q_l == -1) m_R2solutionIsReal = false;
    // no real solutions; take the real part of 1 (real parts are the same)
    nReal = 1;
    m_nComplexRoots++;
  }

  for (unsigned int i = 0; i < nReal; i++) {
    rootR[i] = root[i].real();
    E_nu_R = sqrt(px_nu*px_nu + py_nu*py_nu + rootR[i]*rootR[i]);
    p_nu_R.SetPxPyPzE(px_nu, py_nu, rootR[i], E_nu_R);
    for (int j = 0; j < 2; j++) {
      mblv = (p_b[abs(j)] + p_l + p_nu_R).M();
      // printf("Lepton-reconstructed top mass = %f\n", mblv);
      mjjb = (p_b[abs(j-1)] + p_q[0] + p_q[1]).M();
      // printf("Hadron-reconstructed top mass = %f\n", mjjb);
      // printf("For i = %i, j = %i: m_bjj = %.15le, mblv = %.15le\n", i, j, mjjb, mblv);
      dh = mjjb - m_tmass;
      dl = mblv - m_tmass;
      X2 = dh*dh + dl*dl;
      // printf("For i = %i, j = %i: X2 = %.15le\n", i, j, X2);
      if (it == 0) {
        X2min = X2;
        imin = i;
        jmin = j;
      }
      if (X2 < X2min) {
        X2min = X2;
        imin = i;
        jmin = j;
      }
      it++;
    }
  }

  // printf("Chosen solution: imin = %i, jmin = %i\n", imin, jmin);

  // Assess neutrino reconstruction performance.

  double pz_nu_truth = p_nu.Pz();
  double Root0MinusTruth = abs(root[0].real() - pz_nu_truth);
  double Root1MinusTruth = abs(root[1].real() - pz_nu_truth);
  int bestRoot;
  if (Root0MinusTruth < Root1MinusTruth) bestRoot = 0;
  else if (Root1MinusTruth < Root0MinusTruth) bestRoot = 1;
  else bestRoot = 0;
  if (imin == bestRoot) m_nNeutrinoMatched++;

  // Access q-matching performance
  int b_lep = -999;
  if (Q_l == 1) b_lep = 0;
  if (Q_l == -1) b_lep = 1;

  bool b_match;
  if (b_lep == jmin) b_match = true;
  else b_match = false;
  if (b_match) m_nQuarksMatched++;

  // Print reconstruction performance.
  // printf("True pz_nu = %f\n", p_nu.Pz());
  // printf("Possible neutrino solutions:\n");
  // printf("                             %f + %fi\n", root[0].real(), root[0].imag());
  // printf("                             %f + %fi\n", root[1].real(), root[1].imag());
  // printf("Chosen solution:             %f + %fi\n", root[imin].real(), root[imin].imag());
  // if (imin == bestRoot) printf("Neutrino solution: correct. \n");
  // else printf("Neutrino solution: incorrect. \n");
  // if (b_match) printf("b-assignment: correct. \n");
  // else printf("b-assignment: incorrect. \n");
  // printf("---\n");

  // Fill output vector
  double pz_nu_R = rootR[imin];
  E_nu_R = sqrt(px_nu*px_nu + py_nu*py_nu + pz_nu_R*pz_nu_R);
  p_nu_R.SetPxPyPzE(px_nu, py_nu, pz_nu_R, E_nu_R);
  if (Q_l == 1){
    p_R[0] = p_b[jmin];
    p_R[1] = p_b[abs(jmin-1)];
    p_R[2] = p[2];
    p_R[3] = p_nu_R;
    p_R[4] = p[4];
    p_R[5] = p[5];
  }
  else if (Q_l == -1) {
    p_R[0] = p_b[abs(jmin-1)];
    p_R[1] = p_b[jmin];
    p_R[2] = p[2];
    p_R[3] = p[3];
    p_R[4] = p[4];
    p_R[5] = p_nu_R;
  }
  return p_R;
  printf("finished reconstruction\n");
}


vector<complex<double> > AnalysisZprime::SolveQuadratic(double a, double b, double c) {
    // solves quadratic for both roots
    // returns both as complex values in a complex vector x(2)

    vector<complex<double> > roots;
    complex<double> term1;
    complex<double> term2;
    complex<double> discriminator;

    term1 = -b/(2*a);
    discriminator = b*b - 4*a*c;
    term2 = sqrt(discriminator)/(2*a);

    // print terms
    // printf("term1 = %f + %fi\n", term1.real(), term1.imag());
    // printf("term2 = %f + %fi\n", term2.real(), term2.imag());

    roots.push_back(term1 + term2);
    roots.push_back(term1 - term2);

    return roots;
}

void AnalysisZprime::GetChannelFactors() {

  // branching ratio for t->benu=bmuv=btaumu (1/9 with QCD corrections)
  const double brtbln = 0.10779733;
  // branching ratio for t->bqq is leftover after 3 l generations
  const double brtbqq = 1 - brtbln*3;

  if (m_channel == "tt") {
    // multiply by branching ratios
    // fac_ee = brtbln*brtbln;
    // fac_emu = 2*brtbln*brtbln
    m_sigma = m_sigma*2*brtbln*brtbqq;
    // fac_qq = brtbqq*brtbqq
  }
  else if (m_channel == "bbllnn") {
    // scale dilepton to other classifications
    // fac_ee = 1
    // fac_emu = 2
    m_sigma = m_sigma*12;
    // fac_qq = 36
  }
  // sigma_ee = sigma*fac_ee
  // error_sigma_ee = error_sigma*fac_ee
  // sigma_emu = sigma*fac_emu
  // error_sigma_emu = error_sigma*fac_emu
  // sigma_eq = sigma*fac_eq
  // error_sigma_eq = error_sigma*fac_eq
  // sigma_qq = sigma*fac_qq
  // error_sigma_qq = error_sigma*fac_qq
}


const void AnalysisZprime::UpdateCutflow(int cut, bool passed) {
  if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
  if (passed) m_cutflow[cut] +=1;
}


void AnalysisZprime::InitialiseCutflow() {
  m_cutflow = vector<int>(m_cuts, -999);
  m_cutNames = vector<TString>(m_cuts, "no name");
  m_cutNames[c_entries] = "Entries";
  m_cutNames[c_topDecays] = "t->be#nu";
  m_cutNames[c_antitopDecays] = "#bar{t}->be#nu";
  m_cutNames[c_events] = "Events";
  m_cutNames[c_realSolutions] = "p^{#nu}_{z} real";
  m_cutNames[c_MET] = "MET";
  m_cutNames[c_mtt] = "mtt";
  m_cutNames[c_fiducial] = "Fiducial";

  h_cutflow = new TH1D("Cutflow", "Cutflow", m_cuts, 0, m_cuts);
}


void AnalysisZprime::PrintCutflow() {
  h_cutflow->Write();
  printf("--- Cutflow ---\n");
  for (int cut = 0; cut < m_cuts; cut++) {
    if (m_cutflow[cut] == -999) continue;

    h_cutflow->SetBinContent(cut+1, m_cutflow[cut]);
    h_cutflow->GetXaxis()->SetBinLabel(cut+1, m_cutNames[cut]);

    printf("%s cut: %i pass\n", m_cutNames[cut].Data(), m_cutflow[cut]);
  }
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
