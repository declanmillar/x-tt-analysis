#include "superpose2d.h"

TCanvas* superpose2d(//int ndimensions,
                   TString histName1, TString histTitle1, TString fileName1, 
                   TString histName2, TString histTitle2, TString fileName2, 
                   TString histName3, TString histTitle3, TString fileName3) {

  TFile *f1, *f2, *f3;
  TH2D *h1, *h2, *h3;
  TString name1, name2, name3;
  TString epsFileName;



  // gStyle->SetOptStat(0);
  // gStyle->SetLegendFont(132);
  // gStyle->SetLegendBorderSize(0);
  // gStyle->SetOptTitle(0);
  
  TCanvas *canvas = new TCanvas("test", "test", 1280, 751);
  canvas->cd();

  f1 = new TFile(fileName1, "READ");
  if (!f1->IsOpen()) printf("Failed to open %s\n", fileName1.Data());
  h1 = (TH2D*) f1->Get(histName1);
  name1 = histName1 + '@' + fileName1;
  if (histTitle1 == "NULL") h1->SetTitle(name1);
  else h1->SetTitle(histTitle1);
  h1->SetMarkerStyle(20);
  h1->Draw("LEGO");//("E1x0P");
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->SetTitleOffset(1.2);

  if (histName2 != "NULL") {
    f2 = new TFile(fileName2, "READ");
    if (!f2->IsOpen()) printf("Failed to open %s\n", fileName2.Data());
    name2 = histName2 + '@' + fileName2;
    h2 = (TH2D*) f2->Get(histName2);
    if (histTitle2 == "NULL") h2->SetTitle(name2);
    else h2->SetTitle(histTitle2);
    h2->SetMarkerStyle(20);
    h2->Draw("SAMELEGO");//("E1x0PSAME");
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
  }

  if (histName3 != "NULL") {
    f3 = new TFile(fileName3, "READ");
    if (!f3->IsOpen()) printf("Failed to open %s\n", fileName3.Data());
    h3 = (TH2D*) f3->Get(histName3);
    name3 = histName3 + '@' + fileName3;
    if (histTitle3 == "NULL") h3->SetTitle(name3);
    else h3->SetTitle(histTitle3);
    h3->SetMarkerStyle(20);
    h3->Draw("SAMELEGO");//("E1x0PSAME");
    h3->SetLineColor(kBlue);
    h3->SetMarkerColor(kBlue);
  }

  // normalize histograms
  bool normalize = true;
  if (normalize == true) { 
    h1->Scale(1.0/h1->Integral("width"));
    if (histName2 != "NULL") h2->Scale(1.0/h2->Integral("width"));
    if (histName3 != "NULL") h3->Scale(1.0/h3->Integral("width"));
  }

  // find overlapping volume (histograms must be have the same user ranges and same number of bins)
  double sigPerOverlap = 0;
  if (histName2 != "NULL") {
    double overlap = 0;
    for (int i = 0; i < h1->GetNbinsX(); i++) {
      for (int j = 0; j < h1->GetNbinsY(); j++) {
        int bin1 = h1->GetBin(i, j);
        int bin2 = h2->GetBin(i, j);
        double one = h1->GetBinContent(bin1);
        double two = h2->GetBinContent(bin2);
        if (one <= two) overlap += one;
        else if (one > two) overlap += two;
      }
    }
    printf("Overlapping volume = %f\n", overlap);
    sigPerOverlap = overlap/h2->Integral()*100;
    printf("Signal in overlapping volume: %f\%\n", sigPerOverlap);
  }

  // set range user
  // if(rangeMax == rangeMax && rangeMax == rangeMax) {
  //   h1->GetXaxis()->SetRangeUser(rangeMin, rangeMax); 
  //   if (nFiles > 1) h2->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  //   if (nFiles > 2) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  // }

  // if (logY == true) canvas->SetLogy();

  TLegend* legend = new TLegend(0.50, 0.60, 0.88, 0.88, "");
  legend->SetBorderSize(0);
  legend->AddEntry(h1, h1->GetTitle());
  if (histName2 != "NULL") legend->AddEntry(h2, h2->GetTitle());
  if (histName3 != "NULL") legend->AddEntry(h3, h3->GetTitle());
  legend->Draw();
  std::ostringstream strs;
  strs << sigPerOverlap;
  std::string str = strs.str();
  TString tstr(str);
  TString SigOverlap = "Signal in overlapping volume = " + tstr + "\%";
  TLatex* texBox = new TLatex(0.5,0.5, SigOverlap);
  texBox->Draw();

  return canvas;
}