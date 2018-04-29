#include "analysis.hpp"
#include "boost/program_options.hpp"
#include "atlas-style.hpp"
// #include <boost/timer/timer.hpp>

using namespace std;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    
    // boost::timer::auto_cpu_timer timer;
    
    po::options_description desc("options");
    desc.add_options()
    ("model,m",          po::value<string>()->default_value("SM"))
    ("gg,g",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("qq,q",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("uu,u",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("dd,d",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("final_state,F",    po::value<string>()->default_value("tt-bbeevv"))
    ("add_gg,G",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("add_qq,Q",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("energy,E",         po::value<int>()->default_value(13))
    ("luminosity,L",     po::value<double>()->default_value(-1))
    ("minimumBtags,b",   po::value<int>()->default_value(1))
    ("options,o",        po::value<string>()->default_value(""))
    ("reconstruction,r", po::value<string>()->default_value("NuW"))
    ("nevents_max,n",  po::value<Long64_t>()->default_value(0))
    ("tag,t",            po::value<string>()->default_value(""))
    ("slice,s",          po::value<bool>()->default_value(false)->implicit_value(true))
    ("inputFileName,f",  po::value<string>()->default_value(""))
    ("processFileName,p",po::value<string>()->default_value(""))
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
    auto inputFileName = opt["inputFileName"].as<string>();
    auto processFileName = opt["processFileName"].as<string>();
    auto options = opt["options"].as<string>();
    auto energy = opt["energy"].as<int>();
    auto luminosity = opt["luminosity"].as<double>();
    auto add_gg = opt["add_gg"].as<bool>();
    auto add_qq = opt["add_qq"].as<bool>();
    auto minimumBtags = opt["minimumBtags"].as<int>();
    auto reconstruction = opt["reconstruction"].as<const string>();
    auto nevents_max = opt["nevents_max"].as<const Long64_t>();
    auto tag = opt["tag"].as<string>();
    auto slice = opt["slice"].as<bool>();

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

    if (reconstruction != "TRN" and reconstruction != "KIN" and reconstruction != "NuW") {
        cout << "ERROR: invalid reconstruction \"" << reconstruction << "\"\n";
        return 0;
    }
    
    string reconstruction_long = "";
    if (reconstruction == "TRN") reconstruction_long = "Transverse variables only; no reconstruction";
    else if (reconstruction == "KIN") reconstruction_long = "Minimise invariant top pair mass";
    else if (reconstruction == "NuW") reconstruction_long = "Neutrino Weighting";

    cout << "SETTINGS\n";
    if (inputFileName == "") cout << "Process:          " << process << "\n";
    else cout                     << "Input file:       " << inputFileName << "\n";
    cout << "Minimum b-tags:   " << minimumBtags << "\n";
    cout << "Reconstruction:   " << reconstruction_long << "\n";
    if (options != "") cout << "Options:          " << options << "\n";
    if (luminosity != -1) cout << "Luminosity:       " << luminosity << " [fb-1]\n";
    cout << "Energy:           " << energy << " [TeV]\n";

    if (inputFileName == "") {
        auto analysis = new Analysis(model, process, options, energy, luminosity, minimumBtags, reconstruction, tag, slice);
        analysis->Run();
    }
    else {
        auto analysis = new Analysis(inputFileName, processFileName, luminosity, minimumBtags, reconstruction, nevents_max, tag, slice);
        analysis->Run();
    }

}
