def bins_match(TH1D* h1, TH1D* h2):
  return False if h1.GetNbinsX() != h2.GetNbinsX()
  return False if h1.GetMinimum() != h2.GetMinimum()
  return False if h1.GetMaximum() != h2.GetMaximum()
  return true

def significance(TH1D *h1, TH1D *h2):
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
