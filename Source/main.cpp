#include "analysis.hpp"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  const TString channel("tt"), model("SM"), options("_"), analysisLabel("");
  const int energy = 13, it = 2, points = 5000, btags = 2;
  const double luminosity = -1; // fb-1 (set negative to ignore luminosity)
  const bool discardComplex(false), addQCD(false);

  AnalysisZprime analysis(channel, model, energy, options, it, points, addQCD, luminosity, btags, discardComplex, analysisLabel);
}
