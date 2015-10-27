#include "analysis.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  const TString channel("bbllnn"), model("SM"), options("_xc_"), analysisLabel("");
  const int energy = 13, it = 1, points = 5000, btags = 2;
  const double luminosity = 10; // fb-1 (set negative to ignore luminosity)
  const bool discardComplex(false);

  AnalysisZprime analysis(channel, model, energy, options, it, points, luminosity, btags, discardComplex, analysisLabel);
}
