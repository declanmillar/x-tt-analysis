#include "analysis.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  const TString channel("bbllnn"), model("GSM-SM"), options("_xc_");
  const int energy = 13, it = 5, points = 2000000, btags = 2;
  const double luminosity = 300000;
  const bool discardComplex(false);

  AnalysisZprime analysis(channel, model, energy, options, it, points, luminosity, btags, discardComplex);
  AnalysisZprime analysis2(channel, model, energy, options, it, points, luminosity, btags, true);
}
