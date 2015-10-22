#include "analysis.h"

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
}
