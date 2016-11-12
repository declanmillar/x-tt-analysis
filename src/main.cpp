#include "analysis.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  po::options_description desc("Allowed options");
  desc.add_options()
    ("btags,b", po::value<int>()->default_value(2))
    ("energy,E", po::value<int>()->default_value(13))
    ("fid,F", po::value<bool>()->default_value(false)->implicit_value(true))
    ("final_state,f", po::value<string>()->default_value("tt-bbllvv"))
    ("ggG,g", po::value<bool>()->default_value(false)->implicit_value(true))
    ("qqG,q", po::value<bool>()->default_value(false)->implicit_value(true))
    ("intermediates,I", po::value<string>()->default_value("AZX-"))
    ("initial_state,i", po::value<string>()->default_value("qq"))
    ("luminosity,L", po::value<double>()->default_value(-1))
    ("model,m", po::value<string>()->default_value("SM"))
    ("itmx,N", po::value<int>()->default_value(5))
    ("points,n", po::value<string>()->default_value("10M"))
    ("options,o", po::value<string>()->default_value("_"))
    ("tag,t", po::value<string>()->default_value(""))
    ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true))
    ("xsec,x", po::value<bool>()->default_value(true)->implicit_value(false))
    ("ytt,y", po::value<double>()->default_value(0))
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
  TString analysisLabel = opts["tag"].as<string>();
  int energy = opts["energy"].as<int>();
  int it = opts["itmx"].as<int>();
  string points = opts["points"].as<string>();
  int btags = opts["btags"].as<int>();
  int luminosity = opts["luminosity"].as<double>();
  bool add_ggG = opts["ggG"].as<bool>();
  bool add_qqG = opts["qqG"].as<bool>();
  double ytt = opts["ytt"].as<double>();
  bool xsec = opts["xsec"].as<bool>();
  bool fid = opts["fid"].as<bool>();

  if (model == "SM" && intermediates == "AZX-") intermediates = "AZ-";

  Analysis* analysis = new Analysis(model, initial_state, intermediates, final_state, energy, options, it, points, add_ggG, add_qqG, luminosity, btags, analysisLabel);
  analysis->SetYttCut(ytt);
  analysis->SetXsec(xsec);
  analysis->SetFiducial(fid);
  analysis->Run();
}
