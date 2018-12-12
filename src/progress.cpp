#include "progress.hpp"

void ProgressBar( unsigned int x, unsigned int n, unsigned int w)
{
    if ( ( x != n ) && ( x % ( n / 100 + 1 ) != 0 ) ) return;

    double ratio = x / ( double ) n;
    unsigned int c = ratio * w;

    std::cout << "[";
    for ( unsigned int i = 0; i < c; i++ ) std::cout << "#";
    for ( unsigned int i = c; i < w; i++ ) std::cout << "-";
    std::cout << "] (" << std::setw(3) << ( int )( ratio * 100 ) << "%)";
    if ( x == n or x == 0 ) std::cout << "\n" << std::flush;
    else std::cout << "\r" << std::flush;
}

void ProgressBar2( unsigned int x, unsigned int n, unsigned int w)
{
    if ( ( x != n ) && ( x % ( n / 100 + 1 ) != 0 ) ) return;

    double ratio = x / ( double ) n;
    unsigned int c = ratio * w;

    std::cout << "[";
    for ( unsigned int i = 0; i < c; i++ ) std::cout << "#";
    for ( unsigned int i = c; i < w; i++ ) std::cout << "-";
    std::cout << "] (" << std::setw(3) << ( int )( ratio * 100 ) << "%)";
    if ( x == n or x == 0 ) std::cout << "\n" << std::flush;
    else std::cout << "\n" << std::flush;
}

void ProgressPercentage( unsigned int x, unsigned int n, unsigned int w )
{
    if ( ( x != n ) && ( x % ( n / 100 + 1 ) != 0 ) ) return;

    double ratio = x / ( double ) n;
    unsigned int c = ratio * w;

    std::cout << "Progress: " << std::setw(3) << ( int )( ratio * 100 ) << "%";
    std::cout << "\r" << std::flush;
}
