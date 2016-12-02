#include "analysis.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    po::options_description desc("options");
    desc.add_options()
    ("model,m",          po::value<std::string>()->default_value("SM"))
    ("gg,g",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("qq,q",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("uu,u",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("dd,d",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("final_state,f",    po::value<std::string>()->default_value("tt-bbllvv"))
    ("add_gg,G",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("add_qq,Q",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("energy,E",         po::value<int>()->default_value(13))
    ("luminosity,L",     po::value<double>()->default_value(-1))
    ("options,o",        po::value<std::string>()->default_value(""))
    ("reco,t",           po::value<int>()->default_value(1))
    ;
    po::variables_map opt;
    po::store(po::parse_command_line(argc, argv, desc), opt);

    std::string model = opt["model"].as<std::string>();
    std::string initial_state;
    bool gg = opt["gg"].as<bool>();
    bool qq = opt["qq"].as<bool>();
    bool dd = opt["dd"].as<bool>();
    bool uu = opt["uu"].as<bool>();
    if (gg) initial_state += "gg";
    if (qq) initial_state += "qq";
    if (dd) initial_state += "dd";
    if (uu) initial_state += "uu";
    std::string final_state = opt["final_state"].as<std::string>();
    std::string options = opt["options"].as<std::string>();
    int energy = opt["energy"].as<int>();
    int luminosity = opt["luminosity"].as<double>();
    bool add_gg = opt["add_gg"].as<bool>();
    bool add_qq = opt["add_qq"].as<bool>();
    int reco = opt["reco"].as<int>();

    AtlasROOTStyle atlasStyle;
    atlasStyle.SetStyle();

    std::string intermediates;
    if (boost::contains(initial_state,"uu") or boost::contains(initial_state,"dd")) {
        if (model == "SM") intermediates = "-AZ-";
        else intermediates = "-AZX-";
    }
    else {
        intermediates = "-";
    }

    std::string process = initial_state + intermediates + final_state;

    std::cout << "energy: " << energy << " [TeV]" << std::endl;
    std::cout << "process: " << process << std::endl;
    std::cout << "options: " << options << std::endl;
    std::cout << "luminosity: " << luminosity << std::endl;
    std::cout << "reco: " << reco << std::endl;

    Analysis* analysis = new Analysis(model, process, options, energy, luminosity, reco);
}
