#include "analysis.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("channel,c", po::value<string>()->default_value("bbllnn"), "final state")
      ("model,m", po::value<string>()->default_value("SM"), "model")
      ("options,o", po::value<string>()->default_value("_xc_"), "options")
      ("label,l", po::value<string>()->default_value(""), "labeTS")
      ("energy,e", po::value<int>()->default_value(13), "energy")
      ("it,i", po::value<int>()->default_value(5), "vegas iterations")
      ("points,p", po::value<int>()->default_value(2000000), "vegas points")
      ("btags,b", po::value<int>()->default_value(2), "btags")
      ("luminosity,L", po::value<double>()->default_value(-1), "luminosity")
      ("discard,d", po::value<bool>()->default_value(false), "discard complex")
      ("qcd,q", po::value<bool>()->default_value(true), "add QCD")
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
  const int it = opts["it"].as<int>();
  const int points = opts["points"].as<int>();
  const int btags = opts["btags"].as<int>();
  const int luminosity = opts["luminosity"].as<double>();
  const bool discardComplex = opts["discard"].as<bool>();
  const bool addQCD = opts["qcd"].as<bool>();

  printf("channel = %s\n", channel.Data());

  AnalysisZprime analysis(channel, model, energy, options, it, points, addQCD, luminosity, btags, discardComplex, analysisLabel);
}
