#include "analysis.h"
#include "overlay.h"
#include "superpose2d.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString dataDirectory;
  #if __linux || __linux__
    dataDirectory = "/afs/cern.ch/work/d/demillar/Zp-tt_pheno";
  #elif __APPLE__ || __MACH__
    dataDirectory = "/Users/declan/Data/Zp-tt_pheno";
  #endif

  TString ntupleDirectory = dataDirectory + "/NTuples";
  TString weightsDirectory = dataDirectory + "/Weights";
  TString histogramDirectory = dataDirectory + "/Histograms";

  const double luminosity(300000);
  TString channel("bbllnn"); 
  TString model("SM");
  TString options("_xc_");
  TString energy("13");
  TString points("5x5000");
  TString filename(channel + "_" + model + "_" + energy + options + points);
  TString inputFilename(ntupleDirectory + "/" + channel + "/" + filename + ".root");
  TString weightFilename(weightsDirectory + "/" + channel + "/" + filename + ".txt");
  TString outputFilename(histogramDirectory + "/" + channel + "/" + filename + "_hist.root");

  // TString channel2("bbllnn"); 
  // TString model2("GLR-R");
  // TString options2("_xc_");
  // TString points2("5x5000");
  // TString filename2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  // TString inputFilename2(ntupleDirectory + "/" + channel2 + "/" + filename2 + ".root");
  // TString weightFilename2(weightsDirectory + "/" + channel + "/" + filename2 + ".txt");
  // TString outputFilename2(histogramDirectory + "/" + channel + "/" + filename2 + "_hist.root");

  // TString channel3("2to6"); 
  // TString model3("GLR-R");
  // TString options3("_xc_");
  // TString points3("5x1000000");
  // TString filename3(channel3 + "_" + model3 + "_" + energy + options3 + points3);
  // TString inputFilename3(ntupleDirectory + "/" + channel3 + "/" + filename3);
  // TString outputFilename3(histogramDirectory + filename3 + "_hist.root");

  // TString channel4("2to6"); 
  // TString model4("E6-eta");
  // TString options4("_xc_");
  // TString points4("5x5000000");
  // TString filename4(channel4 + "_" + model4 + "_" + energy + options4 + points4);
  // TString inputFilename4(ntupleDirectory + "/" + channel4 + "/" + filename4);
  // TString outputFilename4(histogramDirectory + "/"  filename4 + "_hist.root");

  AnalysisZprime analysis(channel, model, luminosity, inputFilename, weightFilename, outputFilename);
  // AnalysisZprime analysis2(channel2, model2, luminosity, inputFilename2, weightFilename2, outputFilename2);
  // AnalysisZprime analysis3(channel3, model3, inputFilename3, outputFilename3);
  // AnalysisZprime analysis4(channel4, model4, inputFilename4, outputFilename4);

  // bool truthVsReco(true);

  // TApplication* RootApp = new TApplication("RootApp", &argc, argv);

  // if (channel == "tt") {
  //   TCanvas* c_Mtt = Overlay(false, true, "Mff", model, outputFilename, 
  //                                               "Mff", model2, outputFilename2);
  //                                               // "MttNuReco", model3, outputFilename3);
  //                                               // "MttNuReco", model4, outputFilename4);
  //   c_Mtt->Draw();
  //   RootApp->Run(kTRUE);
  // }

  // if (channel == "bbllnn") {

  //   TCanvas* c_Mtt = Overlay(false, true, "Mff", model, outputFilename, 
  //                                               "Mff", model2, outputFilename2);
  //                                               // "MttNuReco", model3, outputFilename3);
  //                                               // "MttNuReco", model4, outputFilename4);
  //   c_Mtt->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_Mtt_r = Overlay(false, true, "Mtt_r", model, outputFilename, 
  //                                           "Mtt_r", model2, outputFilename2);
  //                                           // "Mtt_r", model3, outputFilename3);
  //                                           // "Mtt_r", model4, outputFilename4);
  //   c_Mtt_r->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_AFBstar = Overlay(false, true, "AFBstar", model, outputFilename, 
  //                                             "AFBstar", model2, outputFilename2);
  //                                             // "AFBstar", model3, outputFilename3);
  //                                             // "AFBstar", model4, outputFilename4);
  //   c_AFBstar->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_AFBstar_r = Overlay(false, true, "AFBstar_r", model, outputFilename, 
  //                                               "AFBstar_r", model2, outputFilename2);
  //                                               // "AFBstar_r", model3, outputFilename3);
  //                                               // "AFBstar_r", model4, outputFilename4);
  //   c_AFBstar_r->Draw();
  //   RootApp->Run(kTRUE);

  //   if (truthVsReco) {
  //     TCanvas* c_PzNuTruthVsReco = Overlay(false, false, "Pz_nu", "Truth", outputFilename, 
  //                                                        "Pz_nu_r", "Reconstructed", outputFilename);
  //     c_PzNuTruthVsReco->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_MttReco = Overlay(false, false, "Mff", "Truth", outputFilename, 
  //                               "Mtt_r", "Reconstructed", outputFilename);
  //     c_MttReco->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_yttReco = Overlay(false, false, "ytt", "Truth", outputFilename, 
  //                               "ytt_r", "Reconstructed", outputFilename);
  //     c_yttReco->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_CosThetaReco = Overlay(false, false, "CosTheta", "Truth", outputFilename, 
  //                               "CosTheta_r", "Reconstructed", outputFilename);
  //     c_CosThetaReco->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_CosThetaStarReco = Overlay(false, false, "CosThetaStar", "Truth", outputFilename, 
  //                               "CosThetaStar_r", "Reconstructed", outputFilename);
  //     c_CosThetaStarReco->Draw();
  //     RootApp->Run(kTRUE);
  //   }    
  // }

  return 0;  
}
