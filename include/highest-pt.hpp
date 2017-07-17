#ifndef _HIGHEST_PT_H_
#define _HIGHEST_PT_H_

#include "TLorentzVector.h"

TLorentzVector HighestPt(const std::vector<TLorentzVector>&);
std::pair<TLorentzVector, TLorentzVector> TwoHighestPt(const std::vector<TLorentzVector>&);

#endif
