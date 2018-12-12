#ifndef SEMILEPTON_RECONSTRUCTER_H
#define SEMILEPTON_RECONSTRUCTER_H
#include <iostream>
#include "TLorentzVector.h"
#include "poly34.h"

class SemileptonReconstructer{

private:
    TLorentzVector m_top;
    TLorentzVector m_tbar;
    TLorentzVector m_ttbar;
    TLorentzVector m_b;
    TLorentzVector m_bbar;
    TLorentzVector m_q;
    TLorentzVector m_qbar;
    TLorentzVector m_nu;
    const double m_bTags;
    const double m_WMass;
    const double m_TopMass;
    const bool m_debug = false;
    void Reset();

public:
    SemileptonReconstructer(int, double, double);
    virtual ~SemileptonReconstructer();
    bool Reconstruct(const TLorentzVector&, double charge, const std::vector<TLorentzVector>&, const std::vector<TLorentzVector>&, const TLorentzVector&);
    TLorentzVector GetTop(){   return m_top;   };
    TLorentzVector GetTbar(){  return m_tbar;  };
    TLorentzVector GetTtbar(){ return m_ttbar; };
    TLorentzVector GetB(){     return m_b;     };
    TLorentzVector GetBbar(){  return m_bbar;  };
    TLorentzVector GetNu(){    return m_nu;    };
};
#endif
