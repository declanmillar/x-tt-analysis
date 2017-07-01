#include "trim.hpp"

std::string trim( std::string const& str )
{
    if ( str.empty() ) return str;
    size_t firstScan = str.find_first_not_of(' ');
    size_t first = firstScan == std::string::npos ? str.length() : firstScan;
    size_t last = str.find_last_not_of(' ');
    return str.substr( first, last - first + 1 );
}
