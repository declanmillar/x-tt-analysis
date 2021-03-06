#include "highest-pt.hpp"

TLorentzVector HighestPt(const std::vector<TLorentzVector>& p_b) {
    // returns the 4-vector with the highest pT from an array of 4-vectors
    int imax = 0;
    for (int i = 1; i < p_b.size(); i++) {
        if (p_b[i].Pt() > p_b[imax].Pt()) {
            imax = i;
        }
    }
    return p_b[imax];
}

std::pair<TLorentzVector, TLorentzVector> TwoHighestPt(const std::vector<TLorentzVector>& p_b) {
    // returns a pair of 4-vectors with the highest pTs from an array of 4-vectors
    std::pair<TLorentzVector, TLorentzVector> p_b_hi;
    std::vector<double> pT_b;
    for (auto p : p_b) {
        pT_b.push_back(p.Pt());
    }
    int imax, imax2;
    if (pT_b[0] > pT_b[1]) {
        imax2 = 1;
        imax = 0;
    } else {
        imax2 = 0;
        imax = 1;
    }
    for (int i = 2; i < pT_b.size(); i++) {
        // use >= n not just > as max and second_max can have same value. E.g. {1, 2, 3, 3}
        if (pT_b[i] >= pT_b[imax]) {
            imax2 = imax;
            imax = i;
        } else if (pT_b[i] > pT_b[imax2]) {
            imax2 = i;
        }
    }
    p_b_hi.first = p_b[imax];
    p_b_hi.second = p_b[imax2];
    return p_b_hi;
}
