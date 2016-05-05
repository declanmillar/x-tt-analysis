#include "analysis.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("channel,f", po::value<string>()->default_value("tt-bbllvv"), "final state")
      ("model,m", po::value<string>()->default_value("SM"), "model")
      ("options,o", po::value<string>()->default_value("_"), "options")
      ("label,l", po::value<string>()->default_value(""), "labeTS")
      ("energy,E", po::value<int>()->default_value(13), "energy")
      ("itmx,N", po::value<int>()->default_value(5), "vegas iterations")
      ("points,n", po::value<string>()->default_value("5M"), "vegas points")
      ("btags,b", po::value<int>()->default_value(2), "btags")
      ("luminosity,L", po::value<double>()->default_value(-1), "luminosity")
      ("discard,d", po::value<bool>()->default_value(false), "discard complex")
      ("qcd,q", po::value<bool>()->default_value(true), "add QCD diagrams")
      ("verbose,v", po::value<bool>()->default_value(false), "run in verbose mode")
  ;
  po::variables_map opts;
  po::store(po::parse_command_line(argc, argv, desc), opts);

  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  const TString channel = opts["channel"].as<string>();
  const TString model = opts["model"].as<string>();
  const TString options = opts["options"].as<string>();
  const TString analysisLabel = opts["label"].as<string>();
  const int energy = opts["energy"].as<int>();
  const int it = opts["itmx"].as<int>();
  const string points = opts["points"].as<string>();
  const int btags = opts["btags"].as<int>();
  const int luminosity = opts["luminosity"].as<double>();
  const bool discardComplex = opts["discard"].as<bool>();
  bool addQCD = opts["qcd"].as<bool>();

  AnalysisZprime analysis(channel, model, energy, options, it, points, addQCD, luminosity, btags, discardComplex, analysisLabel);
}
