#include "analysis.hpp"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    po::options_description desc("options");
    desc.add_options()
    ("model,m",          po::value<string>()->default_value("SM"))
    ("initial_state,i",  po::value<string>()->default_value("qq"))
    ("final_state,f",    po::value<string>()->default_value("tt-bbllvv"))
    ("gg,g",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("qq,q",             po::value<bool>()->default_value(false)->implicit_value(true))
    ("energy,E",         po::value<int>()->default_value(13))
    ("luminosity,L",     po::value<double>()->default_value(-1))
    ("options,o",        po::value<string>()->default_value(""))
    ("tag,t",            po::value<string>()->default_value(""))
    ;
    po::variables_map opt;
    po::store(po::parse_command_line(argc, argv, desc), opt);

    std::string model = opt["model"].as<string>();
    std::string initial_state = opt["initial_state"].as<string>();
    std::string final_state = opt["final_state"].as<string>();
    std::string options = opt["options"].as<string>();
    int energy = opt["energy"].as<int>();
    int luminosity = opt["luminosity"].as<double>();
    bool gg = opt["gg"].as<bool>();
    bool qq = opt["qq"].as<bool>();
    std::string tag = opt["tag"].as<string>();

    AtlasROOTStyle atlasStyle;
    atlasStyle.SetStyle();

    std::string intermediates;
    if (boost::contains(initial_state,"uu") or boost::contains(initial_state,"dd")) {
        if (model == "SM") intermediates = "AZ-";
        else intermediates = "AZX-";
    }
    else {
        intermediates = "-";
    }

    std::string process = initial_state + intermediates + final_state;

    std::cout << "energy: " << energy << " [TeV]" << std::endl;
    std::cout << "process: " << process << std::endl;
    std::cout << "options: " << options << std::endl;
    std::cout << "luminosity: " << luminosity << std::endl;
    std::cout << "tag: " << tag << std::endl;

    Analysis* analysis = new Analysis(model, process, options, gg, qq, energy, luminosity, tag);
}
