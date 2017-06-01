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
    ("final_state,f",    po::value<std::string>()->default_value("tt-bbmumuvv"))
    ("add_gg,G",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("add_qq,Q",         po::value<bool>()->default_value(false)->implicit_value(true))
    ("energy,E",         po::value<int>()->default_value(13))
    ("luminosity,L",     po::value<double>()->default_value(-1))
    ("options,o",        po::value<std::string>()->default_value(""))
    ("reco,r",           po::value<std::string>()->default_value("NuW"))
    ("tag,t",            po::value<std::string>()->default_value(""))
    ;
    po::variables_map opt;
    po::store(po::parse_command_line(argc, argv, desc), opt);

    std::string model = opt["model"].as<std::string>();
    std::string initial_state;
    auto gg = opt["gg"].as<bool>();
    auto qq = opt["qq"].as<bool>();
    auto dd = opt["dd"].as<bool>();
    auto uu = opt["uu"].as<bool>();
    if (gg) initial_state += "gg";
    if (qq) initial_state += "qq";
    if (dd) initial_state += "dd";
    if (uu) initial_state += "uu";
    auto final_state = opt["final_state"].as<std::string>();
    auto options = opt["options"].as<std::string>();
    auto energy = opt["energy"].as<int>();
    auto luminosity = opt["luminosity"].as<double>();
    auto add_gg = opt["add_gg"].as<bool>();
    auto add_qq = opt["add_qq"].as<bool>();
    auto reco = opt["reco"].as<std::string>();
    auto tag = opt["tag"].as<std::string>();

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

    intermediates = "-X-";

    auto process = initial_state + intermediates + final_state;
    std::cout << "Energy:         " << energy << " [TeV]\n";
    std::cout << "Process:        " << process << "\n";
    std::cout << "Options:        " << options << "\n";
    std::cout << "Luminosity:     " << luminosity << " [fb-1]\n";
    std::cout << "Reconstruction: " << reco << "\n";

    auto analysis = new Analysis(model, process, options, energy, luminosity, reco, tag);
    analysis->Run();
}
