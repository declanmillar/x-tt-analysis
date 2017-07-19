#include "analysis.hpp"
#include "boost/program_options.hpp"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    po::options_description desc("options");
    desc.add_options()
    ("model,m",          po::value<string>()->default_value("SM"))
    ("gg,g",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("qq,q",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("uu,u",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("dd,d",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("final_state,f",    po::value<string>()->default_value("tt-bbeevv"))
    ("add_gg,G",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("add_qq,Q",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("energy,E",         po::value<int>()->default_value(13))
    ("luminosity,L",     po::value<double>()->default_value(-1))
    ("minimumBtags,b",   po::value<int>()->default_value(2))
    ("options,o",        po::value<string>()->default_value(""))
    ("reconstruction,r", po::value<string>()->default_value("NuW"))
    ("tag,t",            po::value<string>()->default_value(""))
    ;
    po::variables_map opt;
    po::store(po::parse_command_line(argc, argv, desc), opt);

    string model = opt["model"].as<string>();
    string initial_state;
    auto gg = opt["gg"].as<bool>();
    auto qq = opt["qq"].as<bool>();
    auto dd = opt["dd"].as<bool>();
    auto uu = opt["uu"].as<bool>();
    if (gg) initial_state += "gg";
    if (qq) initial_state += "qq";
    if (dd) initial_state += "dd";
    if (uu) initial_state += "uu";
    auto final_state = opt["final_state"].as<string>();
    auto options = opt["options"].as<string>();
    auto energy = opt["energy"].as<int>();
    auto luminosity = opt["luminosity"].as<double>();
    auto add_gg = opt["add_gg"].as<bool>();
    auto add_qq = opt["add_qq"].as<bool>();
    auto minimumBtags = opt["minimumBtags"].as<int>();
    auto reconstruction = opt["reconstruction"].as<const string>();
    auto tag = opt["tag"].as<string>();

    AtlasROOTStyle atlasStyle;
    atlasStyle.SetStyle();

    string intermediates;
    if (boost::contains(initial_state, "uu") or boost::contains(initial_state, "dd")) {
        if (model == "SM") intermediates = "-AZ-";
        else intermediates = "-AZX-";
    }
    else {
        intermediates = "-";
    }

    auto process = initial_state + intermediates + final_state;

    cout << "Settings\n";
    cout << "Energy:         " << energy << " [TeV]\n";
    cout << "Process:        " << process << "\n";
    if (options != "") cout << "Options:        " << options << "\n";
    if (options != -1) cout << "Luminosity:     " << luminosity << " [fb-1]\n";
    cout << "Reconstruction: " << reconstruction << "\n";

    auto analysis = new Analysis(model, process, options, energy, luminosity, minimumBtags, reconstruction, tag);
    analysis->Run();
}
