#include "analysis.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("final_state,f", po::value<string>()->default_value("tt-bbllvv"), "final state")
      ("initial_state,i", po::value<string>()->default_value("qq"), "initial state")
      ("intermediates,I", po::value<string>()->default_value("AZX"), "intermediate particles")
      ("model,m", po::value<string>()->default_value("SM"), "model")
      ("options,o", po::value<string>()->default_value("_"), "options")
      ("label,l", po::value<string>()->default_value(""), "labeTS")
      ("energy,E", po::value<int>()->default_value(13), "energy")
      ("itmx,N", po::value<int>()->default_value(5), "vegas iterations")
      ("points,n", po::value<string>()->default_value("5M"), "vegas points")
      ("btags,b", po::value<int>()->default_value(2), "btags")
      ("luminosity,L", po::value<double>()->default_value(-1), "luminosity")
      ("discard,d", po::value<bool>()->default_value(false), "discard complex")
      ("ggG,g", po::value<bool>()->default_value(true), "add qqG")
      ("qqG,q", po::value<bool>()->default_value(true), "add qqG")
      ("verbose,v", po::value<bool>()->default_value(false), "run in verbose mode")
      ("ytt,y", po::value<double>()->default_value(0), "set ytt cut")
      ("xsec,x", po::value<bool>()->default_value(true), "calculate differential cross section")
  ;
  po::variables_map opts;
  po::store(po::parse_command_line(argc, argv, desc), opts);

  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString final_state = opts["final_state"].as<string>();
  TString initial_state = opts["initial_state"].as<string>();
  TString intermediates = opts["intermediates"].as<string>();
  TString model = opts["model"].as<string>();
  TString options = opts["options"].as<string>();
  TString analysisLabel = opts["label"].as<string>();
  int energy = opts["energy"].as<int>();
  int it = opts["itmx"].as<int>();
  string points = opts["points"].as<string>();
  int btags = opts["btags"].as<int>();
  int luminosity = opts["luminosity"].as<double>();
  bool discardComplex = opts["discard"].as<bool>();
  bool add_ggG = opts["ggG"].as<bool>();
  bool add_qqG = opts["qqG"].as<bool>();
  double ytt = opts["ytt"].as<double>();
  bool xsec = opts["xsec"].as<bool>();

  if (model == "SM" && intermediates == "AZX") intermediates = "AZ";

  AnalysisZprime* analysis = new AnalysisZprime(model, initial_state, intermediates, final_state, energy, options, it, points, add_ggG, add_qqG, luminosity, btags, discardComplex, analysisLabel);
  analysis->SetYttCut(ytt);
  analysis->SetXsec(xsec);
  analysis->Run();
}
