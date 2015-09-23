#include "analysis.h"


AnalysisZprime::AnalysisZprime(const TString channel, const TString model, const double luminosity, const TString& inputFileName, const TString& weightsFileName, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_luminosity(luminosity),
  m_Wmass(80.23),
  m_tmass(175.0),
  m_channel(channel),
  m_model(model),
  m_inputFileName(inputFileName),
  m_weightsFileName(weightsFileName),
  m_outputFileName(outputFileName),
  m_inputFiles(NULL),
  m_ntup(NULL),
  m_chainNtup(NULL),
  m_outputFile(NULL)
{  
  this->PreLoop();
  this->Loop();
  this->PostLoop();
}


void AnalysisZprime::EachEvent () {
  UpdateCutflow(c_Event, true);

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
    p_r1 = this->ReconstructSemiLeptonic(p,1);
    p_r2 = this->ReconstructSemiLeptonic(p,-1);

    P_r1.SetPxPyPzE(0,0,0,0);
    P_r2.SetPxPyPzE(0,0,0,0);
    for (unsigned int i = 0; i < p.size(); i++) {
      P_r1 += p_r1[i];
      P_r2 += p_r2[i];
    }  
    
    // reconstructed final particle parton CoM variables
    TVector3 V_r1 = -1*P_r1.BoostVector();
    TVector3 V_r2 = -1*P_r2.BoostVector();
    pcm_r1 = vector<TLorentzVector>(p.size());
    pcm_r2 = vector<TLorentzVector>(p.size());
    for (unsigned int i = 0; i < p.size(); i++) {
      pcm_r1[i] = p_r1[i];
      pcm_r2[i] = p_r2[i];
      pcm_r1[i].Boost(V_r1);
      pcm_r2[i].Boost(V_r2);
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
    p_t_r1 = pcm_r1[0] + pcm_r1[2] + pcm_r1[3];
    p_tb_r1 = pcm_r1[1] + pcm_r1[4] + pcm_r1[5];
    p_t_r2 = pcm_r2[0] + pcm_r2[2] + pcm_r2[3];
    p_tb_r2 = pcm_r2[1] + pcm_r2[4] + pcm_r2[5];
  }

  double y_t = p_t.Rapidity();
  double y_tb = p_tb.Rapidity();
  double dy = std::abs(y_t) - std::abs(y_tb);
  double ytt = P.Rapidity();
  double CosTheta = p_t.CosTheta();
  double CosThetaStar = int(ytt/std::abs(ytt))*CosTheta;
  double ytt_r1 = -999;
  double ytt_r2 = -999;
  double CosTheta_r1 = -999;
  double CosTheta_r2 = -999;
  double CosThetaStar_r1 = -999;
  double CosThetaStar_r2 = -999;
  if (m_channel == "bbllnn"){
    ytt_r1 = P_r1.Rapidity();
    ytt_r2 = P_r2.Rapidity();
    CosTheta_r1 = p_t_r1.CosTheta();
    CosTheta_r2 = p_t_r2.CosTheta();
    CosThetaStar_r1 = int(ytt_r1/std::abs(ytt_r1))*CosTheta_r1;
    CosThetaStar_r2 = int(ytt_r2/std::abs(ytt_r2))*CosTheta_r2;
  }
    
  if (this->PassCuts())
  {    
    // re-weight for different iterations
    double it = m_ntup->iteration();
    double weight = m_ntup->weight_eq();
    weight = weight*m_sigma/m_weights[it-1];
    // printf("weight = %.15le\n", weight);

    // convert to TeV
    double Mff = P.M()/1000;
    double Mtt_r1 = P_r1.M()/1000;
    double Mtt_r2 = P_r2.M()/1000;

    // fill histograms (assumes fixed bin width!)
    h_Mff->Fill(Mff, weight*2/h_Mff->GetXaxis()->GetBinWidth(1));
    h_ytt->Fill(ytt, weight*2/h_ytt->GetXaxis()->GetBinWidth(1));
    h_CosTheta->Fill(CosTheta, weight*2/h_CosTheta->GetXaxis()->GetBinWidth(1));
    h_CosThetaStar->Fill(CosThetaStar, weight*2/h_CosThetaStar->GetXaxis()->GetBinWidth(1));

    // asymmetries
    if (CosThetaStar > 0) {
      h_AFBstarF->Fill(Mff, 2*weight/h_AFBstarF->GetXaxis()->GetBinWidth(1));
    }

    if (CosThetaStar < 0) {
      h_AFBstarB->Fill(Mff, 2*weight/h_AFBstarB->GetXaxis()->GetBinWidth(1));
    }

    if (dy > 0) {
      h_AttCF->Fill(Mff, weight/h_AttCF->GetXaxis()->GetBinWidth(1));
    }

    if (dy < 0) {
      h_AttCB->Fill(Mff, weight/h_AttCB->GetXaxis()->GetBinWidth(1));
    }

    if (m_channel == "tt") {
      h_MttLL->Fill(Mff, m_ntup->weightLL()/h_MttLL->GetXaxis()->GetBinWidth(1));
      h_MttLR->Fill(Mff, m_ntup->weightLR()/h_MttLR->GetXaxis()->GetBinWidth(1));
      h_MttRL->Fill(Mff, m_ntup->weightRL()/h_MttRL->GetXaxis()->GetBinWidth(1));
      h_MttRR->Fill(Mff, m_ntup->weightRR()/h_MttRR->GetXaxis()->GetBinWidth(1));
    }    
    else if (m_channel == "bbllnn") {
      h_Pz_nu->Fill(p[3].Pz(), weight*2/h_Pz_nu->GetXaxis()->GetBinWidth(1));
      
      h_ytt_r->Fill(ytt_r1, weight/h_ytt_r->GetXaxis()->GetBinWidth(1));
      h_ytt_r->Fill(ytt_r2, weight/h_ytt_r->GetXaxis()->GetBinWidth(1));

      h_CosTheta_r->Fill(CosTheta_r1, weight/h_CosTheta_r->GetXaxis()->GetBinWidth(1));
      h_CosTheta_r->Fill(CosTheta_r2, weight/h_CosTheta_r->GetXaxis()->GetBinWidth(1));

      h_CosThetaStar_r->Fill(CosThetaStar_r1, weight/h_CosThetaStar_r->GetXaxis()->GetBinWidth(1));
      h_CosThetaStar_r->Fill(CosThetaStar_r2, weight/h_CosThetaStar_r->GetXaxis()->GetBinWidth(1));

      h_Pz_nu_r->Fill(p_r1[3].Pz(), weight/h_Pz_nu_r->GetXaxis()->GetBinWidth(1));
      h_Pz_nu_r->Fill(p_r2[5].Pz(), weight/h_Pz_nu_r->GetXaxis()->GetBinWidth(1));

      h_Mtt_r->Fill(Mtt_r1, weight/h_Mtt_r->GetXaxis()->GetBinWidth(1));
      h_Mtt_r->Fill(Mtt_r2, weight/h_Mtt_r->GetXaxis()->GetBinWidth(1));

      if (CosThetaStar_r1 > 0) {
        h_AFBstar_rF->Fill(Mtt_r1, weight/h_AFBstar_rF->GetXaxis()->GetBinWidth(1));
      }
      if (CosThetaStar_r2 > 0) {
        h_AFBstar_rF->Fill(Mtt_r2, weight/h_AFBstar_rF->GetXaxis()->GetBinWidth(1));
      }

      if (CosThetaStar_r1 < 0) {
        h_AFBstar_rB->Fill(Mtt_r1, weight/h_AFBstar_rB->GetXaxis()->GetBinWidth(1));
      }
      if (CosThetaStar_r2 < 0) {
        h_AFBstar_rB->Fill(Mtt_r2, weight/h_AFBstar_rB->GetXaxis()->GetBinWidth(1));
      }
    }
  }
}


void AnalysisZprime::PostLoop () {
  this->GetResults();
  this->CheckPerformance();
  this->MakeGraphs();
  this->PrintCutflow();
  this->WriteHistograms();
  printf("--- Done ---\n");
}


void AnalysisZprime::GetResults () {
  printf("--- Results ---\n");
  double sigma = h_Mff->Integral("width");
  if (std::abs(sigma - m_sigma) > 10e-11) {
    printf("Cross section from generation and analysis stage do not match!\n");
    printf("sigma_generation = %.15le\n", m_sigma);
    printf("sigma_analysis   = %.15le\n", sigma);
  }
  else printf("sigma = %f [pb]\n", sigma);
  if (m_channel == "tt") this->TotalSpinAsymmetries();
}


void AnalysisZprime::CheckPerformance () {
  // printf("m_nQuarksMatched = %i\n", m_nQuarksMatched);
  // printf("m_nNeutrinoMatched = %i\n", m_nNeutrinoMatched);
  // printf("m_nReco = %i\n", m_nReco);  

  double quarkRecoRatio = m_nQuarksMatched/(double)m_nReco;
  double neutrinoRecoRatio = m_nNeutrinoMatched/(double)m_nReco;
  printf("--- Performance ---\n");
  printf("Quark assignment: %.1f%% correct\n", quarkRecoRatio*100);
  printf("Pz_nu assignment: %.1f%% correct\n", neutrinoRecoRatio*100);

  if (m_channel == "bbllnn") {
    h_AlL = this->Asymmetry("AL", "A_{L}", h_AlLF, h_AlLB);
    h_AllC = this->Asymmetry("AllC", "A^{ll}_{C}", h_AllCF, h_AllCB);
    h_AFBstar_r = this->Asymmetry("AFBstar_r", "A_{FB}^{*} (reco)", h_AFBstar_rF, h_AFBstar_rB);
  }
}


TH1D* AnalysisZprime::MakeALL () {
  TH1D* h_A = (TH1D*) h_MttLL->Clone();
  TH1D* h_B = (TH1D*) h_MttLR->Clone();

  h_A->Add(h_MttRR);
  h_B->Add(h_MttRL);

  double ALL = this->TotalAsymmetry(h_A,h_B);
  printf("ALL' = %f\n", ALL);

  TH1D* h_ALL = this->Asymmetry("ALL", "A_{LL}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_ALL;
}


TH1D* AnalysisZprime::MakeAL () {
  TH1D* h_A = (TH1D*) h_MttLL->Clone();
  TH1D* h_B = (TH1D*) h_MttRR->Clone();

  h_A->Add(h_MttLR);
  h_B->Add(h_MttRL);

  TH1D* h_AL = this->Asymmetry("AL", "A_{L}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_AL;
}


void AnalysisZprime::TotalSpinAsymmetries () {
  double sigmaLL = h_MttLL->Integral("width");
  double sigmaLR = h_MttLR->Integral("width");
  double sigmaRL = h_MttRL->Integral("width");
  double sigmaRR = h_MttRR->Integral("width");

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


TH1D* AnalysisZprime::Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B) {
  TH1D* h_numerator = (TH1D*) h_A->Clone(name);
  TH1D* h_denominator = (TH1D*) h_A->Clone();
  h_numerator->SetTitle(title);
  h_numerator->Add(h_B, -1);
  h_denominator->Add(h_B, 1);
  h_numerator->Divide(h_denominator);
  delete h_denominator;
  if (m_luminosity > 0) this->AsymmetryUncertainty(h_numerator, h_A, h_B);
  return h_numerator;
}


void AnalysisZprime::ApplyLuminosity(TH1D* h) {
  if (!m_useLumi) return;
  // printf("Name: %s\n", h->GetTitle());
  // printf("Luminosity: %f\n", m_luminosity);
  double sigma = -999, N = -999, dN = -999;
  double efficiency = 1.0;
  for (int i = 1; i < h->GetNbinsX(); i++) {
    sigma = h->GetBinContent(i)*h->GetBinWidth(i);
    N = m_luminosity*efficiency*sigma;
    h->SetBinContent(i, N);
    dN = std::sqrt(N);
    h->SetBinError(i, dN);
    // printf("sigma = %f, N = %f, dN = %f\n", sigma, N, dN);
  }
  TString events = "Events";
  h->GetYaxis()->SetTitle(events);
}


void AnalysisZprime::AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B) {
  double efficiency = 1.0;
  double A;
  double sigmaA, sigmaB, sigma;
  double deltaA;
  double N;
  for (int i = 1; i < h_Asymmetry->GetNbinsX(); i++) {
    A = h_Asymmetry->GetBinContent(i);
    sigmaA = h_A->GetBinContent(i);
    sigmaB = h_B->GetBinContent(i);
    sigma = (sigmaA + sigmaB)*(h_Asymmetry->GetBinWidth(i));
    N = m_luminosity*efficiency*sigma;
    if (N > 0) deltaA = std::sqrt((1.0 - A*A)/N);
    else deltaA = 0;
    // printf("A = %f, dA= %f, N= %f\n", A, deltaA, N);
    h_Asymmetry->SetBinError(i, deltaA);
  }
}


void AnalysisZprime::CreateHistograms() {
  h_AFBstarF = new TH1D("AFstar", "m_{tt}^{F*}", 50, 2.0, 4.0);
  h_AFBstarB = new TH1D("ABstar", "m_{tt}^{B*}", 50, 2.0, 4.0);
  h_AttCF = new TH1D("AttCF", "m_{tt}^{CF}", 50, 2.0, 4.0);
  h_AttCB = new TH1D("AttCB", "m_{tt}^{CB}", 50, 2.0, 4.0);
  h_CosTheta = new TH1D("CosTheta", "cos#theta", 50, -1.0, 1.0);
  h_CosThetaStar = new TH1D("CosThetaStar", "cos#theta^{*}", 50, -1.0, 1.0);

  if (m_channel == "ll") {
    h_Mff = new TH1D("Mff", "m_{ll}", 50, 2.0, 4.0);
  }

  if (m_channel == "tt") {
    h_Mff = new TH1D("Mff", "m_{tt}", 25, 2.0, 4.0);
    h_ytt = new TH1D("ytt", "y_{tt}", 25, -2.5, 2.5);
    h_MttLL = new TH1D("MttLL", "m_{tt}^{LL}", 25, 2.0, 4.0);
    h_MttLR = new TH1D("MttLR", "m_{tt}^{LR}", 25, 2.0, 4.0);
    h_MttRL = new TH1D("MttRL", "m_{tt}^{RL}", 25, 2.0, 4.0);
    h_MttRR = new TH1D("MttRR", "m_{tt}^{RR}", 25, 2.0, 4.0);
  }

  if (m_channel == "bbllnn") {
    h_Mff = new TH1D("Mff", "m_{tt}", 25, 2.0, 4.0);
    h_ytt = new TH1D("ytt", "y_{tt}", 25, -2.5, 2.5);
    h_Pz_nu = new TH1D("Pz_nu", "p_{z}^{#nu}", 25,-250.0, 250.0);
    h_CosTheta_r = new TH1D("CosTheta_r", "cos#theta_{reco}", 25, -1.0, 1.0);
    h_CosThetaStar_r = new TH1D("CosThetaStar_r", "cos#theta_{reco}^{*}", 25, -1.0, 1.0);
    h_ytt_r = new TH1D("ytt_r", "y_{tt}^{_r}", 25, -2.5, 2.5);
    h_Pz_nu_r = new TH1D("Pz_nu_r", "p_{z}^{#nu} (reco)", 25, -250.0, 250.0);
    h_Mtt_r = new TH1D("Mtt_r", "M^{reco}_{tt}", 25, 2.0, 4.0);
    h_AlLF = new TH1D("AlLF", "AlLF", 25, 2.0, 4.0);
    h_AlLB = new TH1D("AlLB", "AlLB", 25, 2.0, 4.0);
    h_AllCF = new TH1D("AllCF", "AllCF", 25, 2.0, 4.0);
    h_AllCB = new TH1D("AllCB", "AllCB", 25, 2.0, 4.0);
    h_AFBstar_rF = new TH1D("AFBstarNu_r1", "AFBstarNu_r1", 25, 2.0, 4.0);
    h_AFBstar_rB = new TH1D("AFBstarNu_r2", "AFBstarNu_r2", 25, 2.0, 4.0);
    h_real = new TH1D("real", "m_{tt}", 25, 2.0, 4.0);
    h_imaginary = new TH1D("imaginary", "m_{tt}", 25, 2.0, 4.0);
  }
}



void AnalysisZprime::MakeGraphs() {
  // printf("Making Graphs...\n");
  TString numBase;
  if (m_channel == "tt") numBase = "d#sigma(pp->t#bar{t}) / d"; 
  if (m_channel == "bbllnn") numBase = "d#sigma / d"; //pp->t#bar{t}->b#bar{b}l^{+}l^{-}#nu#bar{#nu}
  TString units = "pb";
  TString TeV = "[TeV]";

  h_Mff->GetXaxis()->SetTitle(h_Mff->GetTitle() + TeV);
  h_Mff->GetYaxis()->SetTitle(numBase + h_Mff->GetTitle() + " [" + units +"/TeV]");
  this->ApplyLuminosity(h_Mff);



  h_AFBstar = this->Asymmetry("AFBstar", "A^{*}_{FB}", h_AFBstarF, h_AFBstarB);
  h_AttC = this->Asymmetry("AttC", "A_{C}", h_AttCF, h_AttCB);

  h_CosTheta->GetYaxis()->SetTitle(numBase + h_CosTheta->GetTitle() + " [" + units +"]");
  h_CosTheta->GetXaxis()->SetTitle(h_CosTheta->GetTitle());
  this->ApplyLuminosity(h_CosTheta);

  h_CosThetaStar->GetYaxis()->SetTitle(numBase + h_CosThetaStar->GetTitle() + " [" + units +"]");
  h_CosThetaStar->GetXaxis()->SetTitle(h_CosThetaStar->GetTitle());
  this->ApplyLuminosity(h_CosThetaStar);

  if (m_channel == "tt") {
    h_MttLL->GetYaxis()->SetTitle(numBase + h_MttLL->GetTitle() + " [" + units +"/GeV]");
    h_MttLR->GetYaxis()->SetTitle(numBase + h_MttLR->GetTitle() + " [" + units +"/GeV]");
    h_MttRL->GetYaxis()->SetTitle(numBase + h_MttRL->GetTitle() + " [" + units +"/GeV]");
    h_MttRR->GetYaxis()->SetTitle(numBase + h_MttRR->GetTitle() + " [" + units +"/GeV]");
    h_ALL = this->MakeALL();
    h_AL = this->MakeAL();
    h_ALL->GetYaxis()->SetTitle(h_ALL->GetTitle());
    h_AL->GetYaxis()->SetTitle(h_AL->GetTitle());
    h_AFBstar->GetYaxis()->SetTitle(h_AFBstar->GetTitle());
    h_AttC->GetYaxis()->SetTitle(h_AttC->GetTitle());
    h_MttLL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttLR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_ALL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AFBstar->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AttC->GetXaxis()->SetTitle("M_{tt} [GeV]");
  }

  if (m_channel == "bbllnn") {

    h_ytt_r->GetYaxis()->SetTitle(numBase + h_ytt_r->GetTitle() + " [" + units +"]");
    h_ytt_r->GetXaxis()->SetTitle(h_ytt_r->GetTitle());
    this->ApplyLuminosity(h_ytt_r);

    h_ytt->GetYaxis()->SetTitle(numBase + h_ytt->GetTitle() + " [" + units +"]");
    h_ytt->GetXaxis()->SetTitle(h_ytt->GetTitle());
    this->ApplyLuminosity(h_ytt);

    h_AFBstar->GetYaxis()->SetTitle(h_AFBstar->GetTitle());
    h_AFBstar->GetXaxis()->SetTitle("M_{tt} [GeV]");

    h_CosTheta_r->GetYaxis()->SetTitle(numBase + h_CosTheta_r->GetTitle() + " [" + units +"]");
    h_CosTheta_r->GetXaxis()->SetTitle(h_CosTheta_r->GetTitle());
    this->ApplyLuminosity(h_CosTheta_r);

    h_CosThetaStar_r->GetYaxis()->SetTitle(numBase + h_CosThetaStar_r->GetTitle() + " [" + units +"]");
    h_CosThetaStar_r->GetXaxis()->SetTitle(h_CosThetaStar_r->GetTitle());
    this->ApplyLuminosity(h_CosThetaStar_r);

    h_AFBstar_r->GetYaxis()->SetTitle(h_AFBstar_r->GetTitle());
    h_AFBstar_r->GetXaxis()->SetTitle("M_{tt}^{reco} [GeV]");

    h_Pz_nu->GetYaxis()->SetTitle(numBase + h_Pz_nu->GetTitle() + " [" + units +"GeV]");
    h_Pz_nu->GetXaxis()->SetTitle(h_Pz_nu->GetTitle());
    this->ApplyLuminosity(h_Pz_nu);

    h_Pz_nu_r->GetYaxis()->SetTitle(numBase + h_Pz_nu_r->GetTitle() + " [" + units +"Gev]");
    h_Pz_nu_r->GetXaxis()->SetTitle(h_Pz_nu_r->GetTitle());
    this->ApplyLuminosity(h_Pz_nu_r);

    h_Mtt_r->GetYaxis()->SetTitle(numBase + h_Mtt_r->GetTitle() + " [" + units +"/TeV]");
    h_Mtt_r->GetXaxis()->SetTitle(h_Mtt_r->GetTitle() + TeV);
    this->ApplyLuminosity(h_Mtt_r);

    h_AllC->GetYaxis()->SetTitle(h_AllC->GetTitle());
    h_AllC->GetXaxis()->SetTitle("M_{tt} [GeV]");
  }
}


void AnalysisZprime::WriteHistograms() {
  m_outputFile->cd();
  m_outputFile->cd("/");

  if (m_channel == "tt" or "bbllnn" or "ll" ) {
    h_Mff->Write();
  }

  if (m_channel == "tt" || m_channel == "bbllnn") {
    h_AFBstar->Write();
    h_AttC->Write();
    h_ytt->Write();
    h_CosTheta->Write();
    h_CosThetaStar->Write();
  }
  
  if (m_channel == "tt") {
    h_MttLL->Write();
    h_MttLR->Write();
    h_MttRL->Write();
    h_MttRR->Write();
    h_ALL->Write();
    h_AL->Write();
  }

  if (m_channel == "bbllnn") {
    // h_AllC->Write();
    h_AFBstar_r->Write();
    h_Pz_nu->Write();
    h_Pz_nu_r->Write();
    h_Mtt_r->Write();
    h_ytt_r->Write();
    h_CosTheta_r->Write();
    h_CosThetaStar_r->Write();
  }
  h_cutflow->Write();
  m_outputFile->Close();
  delete m_outputFile;
}


bool AnalysisZprime::PassCuts() {
  // if (this->PassCuts_Mtt())
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
  // if (Mff > 1200.0)
  //   {
  //     if (Mff < 2200.0)
  //       {
          UpdateCutflow(c_Mtt, true);
          return true;
  //       }
  //   }
  // UpdateCutflow(c_Mtt, false);
  // return false;
}


bool AnalysisZprime::PassCutsFiducial () { 
  for (unsigned int i = 0; i < p.size(); i++) {
    bool outsideCrack = p[i].PseudoRapidity() <= 1.37 || p[i].PseudoRapidity() >= 1.52;
    bool central      = p[i].PseudoRapidity() <= 2.47;
    bool passesFiducialCuts = outsideCrack && central;
    if (passesFiducialCuts == false){
      UpdateCutflow(c_Fiducial, false);
      return false;
    }
    else continue;
  }
  UpdateCutflow(c_Fiducial, true);
  return true;
}


bool AnalysisZprime::PassCutsYtt () { 
  if (std::abs(P.Rapidity()) > 0)
  {
    UpdateCutflow(c_Ytt, true);
    return true;
  }
  UpdateCutflow(c_Mtt, false);
  return false;
}


void AnalysisZprime::PreLoop () {
  printf("--- Setup ---\n");
  this->ResetCounters();
  this->SetupInputFiles();
  this->SetupWeightsFiles();
  this->InitialiseCutflow();
  this->SetupOutputFiles();
  this->CreateHistograms();
}


void AnalysisZprime::ResetCounters () {
  if (m_luminosity >= 0) m_useLumi = true;
  else m_luminosity = false;
  m_nQuarksMatched = 0;
  m_nNeutrinoMatched = 0;
  m_nReco = 0;
  m_nRealRoots = 0;
  m_nComplexRoots = 0;
}


void AnalysisZprime::SetupWeightsFiles () {
  TString weightsName(m_weightsFileName);
  ifstream weights(weightsName.Data());
  if (!weights.is_open()) printf("Error: failed to open %s\n", weightsName.Data());
  weights >> m_sigma;
  string line;
  double weight;
  while ( weights >> weight ) m_weights.push_back(weight);
  weights.close();
}


void AnalysisZprime::Loop () {
  // Loop over all files
  for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i)
  {
    cout << "Input:  '" << (*i) << "'." << endl;
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
  }
}


AnalysisZprime::~AnalysisZprime() { delete m_inputFiles; }


void AnalysisZprime::SetupOutputFiles() {
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  printf("Output: '%s'.\n", m_outputFileName.Data());
}


void AnalysisZprime::SetupInputFiles () {
  m_inputFiles = new vector<TString>;  
  m_inputFiles->push_back(m_inputFileName);
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


std::vector<TLorentzVector> AnalysisZprime::ReconstructSemiLeptonic(std::vector<TLorentzVector> p, int Q_l) {
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
  m_nReco++;

  std::vector<TLorentzVector> p_r(p.size()); // I am returned!
  TLorentzVector p_l, p_nu;
  std::vector<TLorentzVector> p_b(2), p_q(2);

  p_b[0] = p[0];
  p_b[1] = p[1];
  if (Q_l == 1) {
    p_l = p[2];
    p_nu = p[3];
    p_q[0] = p[4];
    p_q[1]= p[5];
  }
  else if (Q_l == -1) {
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
  std::vector<std::complex<double> > root;
  double a = -999, b = -999, c = -999, k = -999;

  E_l = std::sqrt(px_l*px_l + py_l*py_l + pz_l*pz_l);
  if (std::abs(E_l - p_l.E()) > 0.00001) printf("ERROR: Lepton energy doesn't match.\n");

  k = m_Wmass*m_Wmass/2 + px_l*px_nu + py_l*py_nu;
  a = px_l*px_l + py_l*py_l;
  b = -2*k*(pz_l);
  c = (px_nu*px_nu + py_nu*py_nu)*E_l*E_l - k*k;

  root = this->SolveQuadratic(a, b, c);

  // select single solution and match 'jets'
  double X2 = -999, X2min = -999;
  TLorentzVector p_nu_r;
  double dh = -999, dl = -999, E_nu_r = -999, mblv = -999, mjjb = -999;
  int imin = -999, jmin = -999, it = 0;
  std::vector<double> rootR(root.size());
  unsigned int nReal;
  double Mtt_r1 = P_r1.M()/1000;
  double Mtt_r2 = P_r2.M()/1000;
  if (root[0].imag() == 0 and root[1].imag() == 0) {
    nReal = 2; // Two real solutions: pick best match.
    m_nRealRoots++;
    h_real_r1->Fill(Mtt_r1, m_ntuple->weight()/h_Mtt_r1->GetXaxis()->GetBinWidth(1));
    h_real_r2->Fill(Mtt_r2, m_ntuple->weight()/h_Mtt_r2->GetXaxis()->GetBinWidth(1));
  }
  else {
    nReal = 1; // No real solutions: take the real part of 1 (real parts are the same)
    m_nComplexRoots++;
    h_imaginary_r1->Fill(Mtt_r1, m_ntuple->weight()/h_Mtt_r1->GetXaxis()->GetBinWidth(1));
    h_imaginary_r2->Fill(Mtt_r2, m_ntuple->weight()/h_Mtt_r2->GetXaxis()->GetBinWidth(1));
  }

  for (unsigned int i = 0; i < nReal; i++) {
    rootR[i] = root[i].real();
    E_nu_r = sqrt(px_nu*px_nu + py_nu*py_nu + rootR[i]*rootR[i]);
    p_nu_r.SetPxPyPzE(px_nu, py_nu, rootR[i], E_nu_r);
    for (int j = 0; j < 2; j++) {
      mblv = (p_b[std::abs(j)] + p_l + p_nu_r).M();
      mjjb = (p_b[std::abs(j-1)] + p_q[0] + p_q[1]).M();
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
  double Root0MinusTruth = std::abs(root[0].real() - pz_nu_truth);
  double Root1MinusTruth = std::abs(root[1].real() - pz_nu_truth);
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
  double pz_nu_r = rootR[imin];
  E_nu_r = std::sqrt(px_nu*px_nu + py_nu*py_nu + pz_nu_r*pz_nu_r);
  p_nu_r.SetPxPyPzE(px_nu, py_nu, pz_nu_r, E_nu_r);
  if (Q_l == 1){
    p_r[0] = p_b[jmin];
    p_r[1] = p_b[std::abs(jmin-1)];
    p_r[2] = p[2];
    p_r[3] = p_nu_r;
    p_r[4] = p[4];
    p_r[5] = p[5];
  }
  else if (Q_l == -1) {
    p_r[0] = p_b[std::abs(jmin-1)];
    p_r[1] = p_b[jmin];
    p_r[2] = p[2];
    p_r[3] = p[3];
    p_r[4] = p[4];
    p_r[5] = p_nu_r;
  }
  return p_r;
}


std::vector<std::complex<double> > AnalysisZprime::SolveQuadratic(double a, double b, double c) {
    // solves quadratic for both roots 
    // returns both as complex values in a complex vector x(2)

    std::vector<std::complex<double> > roots;
    std::complex<double> term1;
    std::complex<double> term2;
    std::complex<double> discriminator;

    term1 = -b/(2*a);
    discriminator = b*b - 4*a*c;
    term2 = std::sqrt(discriminator)/(2*a);

    // print terms
    // printf("term1 = %f + %fi\n", term1.real(), term1.imag());
    // printf("term2 = %f + %fi\n", term2.real(), term2.imag());

    roots.push_back(term1 + term2);
    roots.push_back(term1 - term2);

    return roots;
}


const void AnalysisZprime::UpdateCutflow(int cut, bool passed) {
  if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
  if (passed) m_cutflow[cut] +=1;
}


void AnalysisZprime::InitialiseCutflow() {
  m_cutflow = std::vector<int>(m_cuts, -999);
  m_cutNames = std::vector<TString>(m_cuts, "no name");
  m_cutNames[c_Event] = "Event";
  m_cutNames[c_MET] = "MET";
  m_cutNames[c_Mtt] = "Mff";
  m_cutNames[c_Fiducial] = "Fiducial";

  h_cutflow = new TH1D("Cutflow", "Cutflow", m_cuts, 0, m_cuts);
}


void AnalysisZprime::PrintCutflow() {
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

  cout << std::setw(3) << (int)(ratio*100) << "% [";
  for (unsigned int i = 0; i < c; i++) cout << "=";
  for (unsigned int i = c; i < w; i++) cout << " ";
  if (x == n) cout << "\n" << std::flush;
  else cout << "]\r" << std::flush;
}
