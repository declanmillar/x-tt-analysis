#include "match-bjets-to-leps.hpp"
#include "iostream"

using namespace std;

pair<TLorentzVector, TLorentzVector> MatchBjetsToLeps(const pair<TLorentzVector, TLorentzVector>& p_l, const pair<TLorentzVector, TLorentzVector>& p_b)
{
    pair<TLorentzVector, TLorentzVector> p_b_match;

    double deltaR[4];
    deltaR[0] = p_l.first.DeltaR(p_b.first);
    deltaR[1] = p_l.first.DeltaR(p_b.second);
    deltaR[2] = p_l.second.DeltaR(p_b.first);
    deltaR[3] = p_l.second.DeltaR(p_b.second);

    int imin = 0;
    for (int i = 1; i < 4; i++)
    {
        if (deltaR[i] < deltaR[imin])
        {
            imin = i;
        }
    }

    if (imin == 0 or imin == 3)
    {
        p_b_match.first = p_b.first;
        p_b_match.second = p_b.second;
    }
    else if (imin == 1 or imin == 2)
    {
        p_b_match.first = p_b.second;
        p_b_match.second = p_b.first;
    }
    return p_b_match;
}
