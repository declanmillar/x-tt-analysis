#include "progress-bar.hpp"
#include <iostream>
#include <iomanip>

void ProgressBar(unsigned int x, unsigned int n, unsigned int w)
{
    if ((x != n) && (x % (n / 100 + 1) != 0)) return;

    float ratio = x / (float) n;
    unsigned int c = ratio * w;

    std::cout << "progress: " << std::setw(3) << (int)(ratio * 100) << '%' << '[';
    for (unsigned int i = 0; i < c; i++) std::cout << '=';
    for (unsigned int i = c; i < w; i++) std::cout << ' ';
    if (x == n) std::cout << '\n' << std::flush;
    else std::cout << ']' << '\r' << std::flush;
}

void ProgressPercentage(unsigned int x, unsigned int n, unsigned int w)
{
    if ((x != n) && (x % (n / 100 + 1) != 0)) return;

    float ratio = x / (float) n;
    unsigned int c = ratio * w;

    std::cout << "progress: " << std::setw(3) << (int)(ratio * 100) << '%';
    std::cout << '\r' << std::flush;
}