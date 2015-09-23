#include "overlay.h"

bool BinsMatch(TH1D* h1, TH1D* h2){
  if (h1->GetNbinsX() != h2->GetNbinsX()) return false;
  if (h1->GetMinimum() != h2->GetMinimum()) return false;
  if (h1->GetMaximum() != h2->GetMaximum()) return false;
  return true;
}

TH1D *Significance(TH1D *h1, TH1D *h2){
  if (!BinsMatch(h1, h2)){
    printf("Bins do not match\n");
  }
  TString name = (TString) h1->GetName() + "_sig";
  TH1D* h = (TH1D*) h1->Clone(name);
  h->Add(h2, -1);
  double error1, error2, error;
  for (int i = 0; i < h->GetNbinsX(); i++){
    error1 = h1->GetBinError(i);
    error2 = h2->GetBinError(i);
    error = std::sqrt(error1*error1 + error2*error2);
    h->SetBinContent(i, h->GetBinContent(i)/error);
  }
  h->GetYaxis()->SetTitle("Significance");
  double labelSize = h1->GetXaxis()->GetLabelSize();
  double titleSize = h1->GetXaxis()->GetTitleSize();
  double titleOffset = h1->GetXaxis()->GetTitleOffset();
  // printf("x label size: %f \n", xLabelSize);
  h->GetXaxis()->SetLabelSize(labelSize*3.2);
  h->GetXaxis()->SetTitleSize(titleSize*3.2);
  h->GetXaxis()->SetTitleOffset(titleOffset/1.5);

  h->GetYaxis()->SetLabelSize(labelSize*2);
  h->GetYaxis()->SetTitleSize(titleSize*3.2);
  h->GetYaxis()->SetTitleOffset(titleOffset/3.2);
  // h->GetYaxis()->SetNdivisions(3);


  return h;
}

TCanvas* Overlay(const bool normalise, const bool findSignificance,
                 TString histName1, TString histTitle1, TString fileName1, 
                 TString histName2, TString histTitle2, TString fileName2, 
                 TString histName3, TString histTitle3, TString fileName3,
                 TString histName4, TString histTitle4, TString fileName4) {

  bool plot2(histName2 != "NULL");
  bool plot3(histName3 != "NULL");
  bool plot4(histName4 != "NULL");

  TFile *f1, *f2, *f3, *f4;
  TH1D *h1, *h2, *h3, *h4;
  TH1D *h2sig, *h3sig, *h4sig;
  TString name1, name2, name3, name4;
  TString epsFileName;

  const bool showErrors(true);
  TString drawOption;
  if (showErrors) drawOption = "E1x0P";
  else drawOption = "HIST";

  TCanvas *canvas = new TCanvas(histName1, histTitle1, 1280, 751);
  canvas->cd();

  double leftMargin = 0.2, rightMargin = 0.1; 

  double minValue, newMinValue, maxValue, newMaxValue;
  TPad *lowerPad, *upperPad;

  if (findSignificance) {
    // double 
    upperPad = new TPad("upperPad", "upperPad", 0, 0.2, 1, 1);
    upperPad->SetFillColor(-1);
    upperPad->Draw();
    upperPad->SetRightMargin(rightMargin);
    upperPad->SetLeftMargin(leftMargin);

    lowerPad = new TPad("lowerPad", "lowerPad", 0, 0.05, 1, 0.3);
    lowerPad->SetFillColor(-1);
    lowerPad->SetTopMargin(0);
    lowerPad->Draw();
    lowerPad->SetBottomMargin(0.3);
    lowerPad->SetTopMargin(0);
    lowerPad->SetRightMargin(rightMargin);
    lowerPad->SetLeftMargin(leftMargin);

    upperPad->cd();
  }

  f1 = new TFile(fileName1, "READ");
  if (!f1->IsOpen()) printf("Failed to open %s\n", fileName1.Data());
  h1 = (TH1D*) f1->Get(histName1);
  h1->SetTitle(histTitle1);
  h1->SetMarkerStyle(20);
  h1->Draw(drawOption);
  h1->SetMarkerColor(kCyan+4);
  h1->SetLineColor(kCyan+4);
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->SetTitleOffset(1.2);
  maxValue = h1->GetBinContent(h1->GetMaximumBin());
  minValue = h1->GetBinContent(h1->GetMinimumBin());

  double xLabelSize = h1->GetXaxis()->GetLabelSize();
  printf("x label size: %f \n", xLabelSize);


  if (plot2) {
    f2 = new TFile(fileName2, "READ");
    if (!f2->IsOpen()) printf("Error: failed to open %s.\n", fileName2.Data());
    h2 = (TH1D*) f2->Get(histName2);
    h2->SetTitle(histTitle2);
    h2->SetMarkerStyle(20);
    h2->Draw(drawOption + "SAME");
    h2->SetLineColor(kPink-8);
    h2->SetMarkerColor(kPink-8);
    newMinValue = h2->GetBinContent(h2->GetMinimumBin());
    if (newMinValue < minValue) minValue = newMinValue;
    newMaxValue = h2->GetBinContent(h2->GetMaximumBin());
    if (newMaxValue > maxValue) maxValue = newMaxValue;
  }

  if (plot3) {
    f3 = new TFile(fileName3, "READ");
    if (!f3->IsOpen()) printf("Error: failed to open %s.\n", fileName3.Data());
    h3 = (TH1D*) f3->Get(histName3);
    name3 = histName3 + '@' + fileName3;
    if (plot3) h3->SetTitle(name3);
    else h3->SetTitle(histTitle3);
    h3->SetMarkerStyle(20);
    h3->Draw(drawOption + "SAME");
    h3->SetLineColor(kAzure+5);
    h3->SetMarkerColor(kAzure+5);
    newMinValue = h3->GetBinContent(h3->GetMinimumBin());
    if (newMinValue < minValue) minValue = newMinValue;
    newMaxValue = h3->GetBinContent(h3->GetMaximumBin());
    if (newMaxValue > maxValue) maxValue = newMaxValue;
  }

  if (plot4) {
    f4 = new TFile(fileName4, "READ");
    if (!f4->IsOpen()) printf("Error: failed to open %s\n", fileName4.Data());
    h4 = (TH1D*) f4->Get(histName4);
    h4->SetTitle(histTitle4);
    h4->SetMarkerStyle(20);
    h4->Draw(drawOption + "SAME");
    h4->SetLineColor(kGreen+2);
    h4->SetMarkerColor(kGreen+2);
    newMinValue = h4->GetBinContent(h4->GetMinimumBin());
    if (newMinValue < minValue) minValue = newMinValue;
    newMaxValue = h4->GetBinContent(h4->GetMaximumBin());
    if (newMaxValue > maxValue) maxValue = newMaxValue;
  }

  // adjust y-axis range
  maxValue *= 1.2;
  if (minValue < 0) minValue *= 1.2;
  else if (minValue > 0) minValue = minValue - minValue*0.2;
  else minValue = 0;

  h1->GetYaxis()->SetRangeUser(minValue, maxValue);

  if (findSignificance) {
    if (plot2) h2sig = Significance(h2, h1);
    if (plot3) h3sig = Significance(h3, h1);
    if (plot4) h4sig = Significance(h4, h1);
  }

  // normalize histograms
  if (normalise) { 
    TString yTitle;
    yTitle = h1->GetYaxis()->GetTitle();
    yTitle = "1/#sigma #times " + yTitle;
    h1->GetYaxis()->SetTitle(yTitle);
    h1->Scale(1.0/std::abs(h1->Integral()));
    if (plot2) h2->Scale(std::abs(1.0/h2->Integral()));
    if (plot3) h3->Scale(std::abs(1.0/h3->Integral()));
    if (plot4) h4->Scale(std::abs(1.0/h4->Integral()));
  }

  bool overlap = false;
  double sigPerOverlap = 0;
  // find overlapping area (histograms must be have the same user ranges and same number of bins)
  if (overlap == true) {
    if (plot3) {
      double overlapPerBin = 0;
      double overlap = 0;
      for (int i = 0; i < h1->GetNbinsX(); i++) {
        int bin1 = h1->GetBin(i);
        int bin2 = h2->GetBin(i);
        double one = h1->GetBinContent(bin1);
        double two = h2->GetBinContent(bin2);
        if (one <= two) overlapPerBin = one;
        else if (one > two) overlapPerBin = two;
        if (one <= two) overlap += one;
        else if (one > two) overlap += two;
        printf("Overlap per bin = %f\n", overlapPerBin);
      }
      printf("Overlapping area = %f\n", overlap);
      sigPerOverlap = overlap/h2->Integral()*100;
      printf("Signal in overlapping area: %f%%\n", sigPerOverlap);
    }
  }

  // double rangeMin = -999;
  // double rangeMax = -999;
  // rangeMin = 2000;
  // rangeMax = 4000;
  // if(rangeMin != -999 && rangeMax != -999) {
  //   h1->GetXaxis()->SetRangeUser(rangeMin, rangeMax); 
  //   if (plot3) h2->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  //   if (plot3) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  //   if (plot4) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  // }

  // double yRangeMin = -999;
  // double yRangeMax = -999;
  // yRangeMin = -1;
  // yRangeMax = 1;
  // if(yRangeMin != -999 && yRangeMax != -999) {
  //   h1->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax); 
  //   if (plot3) h2->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
  //   if (plot3) h3->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
  //   if (plot4) h3->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
  // }


  TLegend* legend = new TLegend(0.70, 0.70, 0.88, 0.88, "");
  legend->SetBorderSize(0);
  legend->AddEntry(h1, h1->GetTitle());
  if (plot2) legend->AddEntry(h2, h2->GetTitle());
  if (plot3) legend->AddEntry(h3, h3->GetTitle());
  if (plot4) legend->AddEntry(h4, h4->GetTitle());
  legend->Draw();
  if (overlap) {
    std::ostringstream strs;
    strs << sigPerOverlap;
    std::string str = strs.str();
    TString tstr(str);
    TString SigOverlap = "Signal in overlapping area = " + tstr + "%%";
    TLatex* texBox = new TLatex(0.5,0.5, SigOverlap);
    texBox->Draw();
  }

  if (findSignificance) {
    lowerPad->cd();
    if (plot2) h2sig->Draw("HIST");
    if (plot3) h3sig->Draw("HIST SAME");
    if (plot4) h4sig->Draw("HIST SAME");
   h1->GetXaxis()->SetLabelSize(0);
  }

  return canvas;
}