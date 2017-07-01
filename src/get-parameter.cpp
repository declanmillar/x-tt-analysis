#include "get-parameter.hpp"

double get_parameter( std::ifstream* file )
{
    if ( file->is_open() )
    {
        std::string temp;
        getline( *file, temp );
        size_t pos = temp.find('!');
        std::istringstream read( pos == std::string::npos ? temp : temp.substr( 0, pos ) );
        std::string line = read.str();
        double parameter = std::stod( line );
        return parameter;
    }
    throw std::invalid_argument( "file" );
}
