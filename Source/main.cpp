#include "analysis.hpp"

// namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  const TString channel("bbllnn"), model("SM"), options("_CFxc_"), analysisLabel("");
  const int energy = 13, it = 5, points = 2000000, btags = 2;
  const double luminosity = -1; // fb-1 (set negative to ignore luminosity)
  const bool discardComplex(false), addQCD(false);

  AnalysisZprime analysis(channel, model, energy, options, it, points, addQCD, luminosity, btags, discardComplex, analysisLabel);
}
