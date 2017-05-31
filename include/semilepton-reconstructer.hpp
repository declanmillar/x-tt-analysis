#ifndef SEMILEPTON_RECONSTRUCTER_H
#define SEMILEPTON_RECONSTRUCTER_H
#include <iostream>
#include "TLorentzVector.h"
#include "solve-poly.hpp"

class SemileptonReconstructer{

private:
    TLorentzVector m_top;
    TLorentzVector m_tbar;
    TLorentzVector m_ttbar;
    TLorentzVector m_b;
    TLorentzVector m_bbar;
    TLorentzVector m_nu;
    const double m_WMass;
    const double m_TopMass;
    const bool m_debug = false;
    void Reset();

public:
    SemileptonReconstructer(double, double);
    virtual ~SemileptonReconstructer();
    bool Reconstruct(const std::pair<TLorentzVector, TLorentzVector>&, const std::vector<TLorentzVector>&, const std::vector<TLorentzVector>&, const TLorentzVector&);
    TLorentzVector GetTop(){   return m_top;   };
    TLorentzVector GetTbar(){  return m_tbar;  };
    TLorentzVector GetTtbar(){ return m_ttbar; };
    TLorentzVector GetB(){     return m_b;     };
    TLorentzVector GetBbar(){  return m_bbar;  };
    TLorentzVector GetNu(){    return m_nu;    };
};
#endif
