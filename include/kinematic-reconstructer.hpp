#ifndef KINEMATIC_RECONSTRUCTER_H
#define KINEMATIC_RECONSTRUCTER_H
#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#include "solve-poly.hpp"
#include "solve-quartic.hpp"

class KinematicReconstructer{

private:
    TLorentzVector m_top;
    TLorentzVector m_tbar;
    TLorentzVector m_ttbar;
    TLorentzVector m_b;
    TLorentzVector m_bbar;
    TLorentzVector m_nu;
    TLorentzVector m_nubar;
    const double m_bMass;
    const double m_WMass;
    const double m_TopMass;
    const bool m_debug = false;
    std::vector<std::pair<TLorentzVector, TLorentzVector> > GetSolutions(const TLorentzVector& p_l1, const TLorentzVector& p_l2, const TLorentzVector& p_b1, const TLorentzVector& p_b2, const TLorentzVector& p_miss);
    void Reset();

public:
    KinematicReconstructer(double, double, double);
    virtual ~KinematicReconstructer();
    bool Reconstruct(const std::pair<TLorentzVector, TLorentzVector>&, const std::vector<TLorentzVector>&, const TLorentzVector&);
    TLorentzVector GetTop(){   return m_top;   };
    TLorentzVector GetTbar(){  return m_tbar;  };
    TLorentzVector GetTtbar(){ return m_ttbar; };
    TLorentzVector GetB(){     return m_b;     };
    TLorentzVector GetBbar(){  return m_bbar;  };
    TLorentzVector GetNu(){    return m_nu;    };
    TLorentzVector GetNubar(){ return m_nubar; };
};
#endif
