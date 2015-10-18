#include "analysis.h"
#include "overlay.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString channel("bbllnn"), model("GSM-SM"), options("_xc_");
  int energy = 13; //TeV
  int it = 5, points = 2000000;
  const double luminosity = 300000;
  const int btags = 2;
  const bool discardComplex(false);

  AnalysisZprime analysis(channel, model, energy, options, it, points, luminosity, btags, discardComplex);
  AnalysisZprime analysis2(channel, model, energy, options, it, points, luminosity, btags, true);

  TString output1 = analysis.GetOutputFilename();
  TString output2 = analysis2.GetOutputFilename();
  TString label1("Real(complex)");
  TString label2("Discard complex");

  TApplication* RootApp = new TApplication("RootApp", &argc, argv);

  if (channel == "bbllnn") {

    TCanvas* c_cutflow = Overlay(false, false, "Cutflow", "Cutflow", output1, "Cutflow", "Cutflow", output2);
    c_cutflow->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_Mtt = Overlay(false, false, "Mff", "Real(complex)", output1, "Mff", "Discard complex", output2);
    c_Mtt->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_Mtt_r = Overlay(false, false, "Mtt_r", label1, output1, "Mtt_r", label2, output2);
    c_Mtt_r->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AFBstar = Overlay(false, false, "AFBstar", label1, output1, "AFBstar", label2, output2);
    c_AFBstar->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AFBstar_r = Overlay(false, false, "AFBstar_r", label1, output1, "AFBstar_r", label2, output2);
    c_AFBstar_r->Draw();
    RootApp->Run(kTRUE);

    // TCanvas* c_PzNuTruthVsReco = Overlay(false, false, "Pz_nu", "Truth", output,
    //                                                    "Pz_nu_r", "Reconstructed", output);
    // c_PzNuTruthVsReco->Draw();
    // RootApp->Run(kTRUE);
    //
    // TCanvas* c_MttReco = Overlay(false, false, "Mff", "Truth", output,
    //                           "Mtt_r", "Reconstructed", output);
    // c_MttReco->Draw();
    // RootApp->Run(kTRUE);
    //
    // TCanvas* c_yttReco = Overlay(false, false, "ytt", "Truth", output,
    //                           "ytt_r", "Reconstructed", output);
    // c_yttReco->Draw();
    // RootApp->Run(kTRUE);
    //
    // TCanvas* c_CosThetaReco = Overlay(false, false, "CosTheta", "Truth", output,
    //                           "CosTheta_r", "Reconstructed", output);
    // c_CosThetaReco->Draw();
    // RootApp->Run(kTRUE);
    //
    // TCanvas* c_CosThetaStarReco = Overlay(false, false, "CosThetaStar", "Truth", output,
    //                           "CosThetaStar_r", "Reconstructed", output);
    // c_CosThetaStarReco->Draw();
    // RootApp->Run(kTRUE);
  }

  return 0;
}
