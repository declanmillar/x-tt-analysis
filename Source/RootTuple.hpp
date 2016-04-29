//  ***************************************************************************
//  *                                                                         *
//  *   This program is free software; you can redistribute it and/or modify  *
//  *   it under the terms of the GNU General Public License as published by  *
//  *   the Free Software Foundation; either version 2 of the License, or     *
//  *   (at your option) any later version.                                   *
//  *                                                                         *
//  *   Author: John Morris (john.morris@cern.ch)                             *
//  *           Queen Mary University of London                               *
//  *   Editor: Declan Millar (declan.millar@cern.ch)                         *
//  *   File Generated on Wed Jun  3 10:31:34 2015                            *
//  *                                                                         *
//  ***************************************************************************/

// This class is for accessing the RootTuple Ntuples
// You should not have to edit this file

#ifndef _NTUPLE_ROOTTUPLE_H_
#define _NTUPLE_ROOTTUPLE_H_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <map>
using std::map;
#include <iostream>
using std::cout;
using std::endl;

class RootTuple{
  public:
    explicit RootTuple(TTree* tree);
    RootTuple(TTree* tree,const bool& isMC,const bool& isAFII);
    virtual ~RootTuple();
    Long64_t totalEvents();
    Long64_t LoadTree(Long64_t entry);

    // public inline member functions -- Use these to get access to the TTree variables
    inline vector<double>*  E() const {b_E->GetEntry(m_currentEvent);return m_E;}
    inline vector<double>*  Px() const {b_Px->GetEntry(m_currentEvent);return m_Px;}
    inline vector<double>*  Py() const {b_Py->GetEntry(m_currentEvent);return m_Py;}
    inline vector<double>*  Pz() const {b_Pz->GetEntry(m_currentEvent);return m_Pz;}
    inline vector<int>*  barcode() const {b_barcode->GetEntry(m_currentEvent);return m_barcode;}
    inline Int_t  iteration() const {b_iteration->GetEntry(m_currentEvent);return m_iteration;}
    inline Double_t  weight() const {b_weight->GetEntry(m_currentEvent);return m_weight;}

    inline Long64_t currentEvent() const {return m_currentEvent;}

  protected:
    Int_t    GetEntry(Long64_t entry);
    void     Init(TTree *tree);

  private:
    RootTuple();
    RootTuple(const RootTuple& rhs);
    void operator=(const RootTuple& rhs);

    bool m_isMC;
    bool m_isAFII;

    TTree          *fChain;
    int             fCurrent;

    Long64_t m_currentEvent;

    vector<double>*  m_E;
    vector<double>*  m_Px;
    vector<double>*  m_Py;
    vector<double>*  m_Pz;
    vector<int>*  m_barcode;
    Int_t  m_iteration;
    Double_t  m_weight;

    TBranch*  b_E;
    TBranch*  b_Px;
    TBranch*  b_Py;
    TBranch*  b_Pz;
    TBranch*  b_barcode;
    TBranch*  b_iteration;
    TBranch*  b_weight;
};
#endif
