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

#include "RootTuple.h" 

RootTuple::RootTuple(TTree* tree) : 
  m_isMC(false), 
  m_isAFII(false) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

RootTuple::RootTuple(TTree* tree,const bool& isMC,const bool& isAFII) : 
  m_isMC(isMC), 
  m_isAFII(isAFII) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

RootTuple::~RootTuple(){ 
  if (!fChain) return; 
  delete m_E; 
  delete m_Px; 
  delete m_Py; 
  delete m_Pz; 
  delete m_barcode; 
} 

int RootTuple::GetEntry(Long64_t entry){ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
} 

Long64_t RootTuple::totalEvents(){
  return fChain->GetEntriesFast();
} 

Long64_t RootTuple::LoadTree(Long64_t entry){ 
  m_currentEvent = entry; 
  if (!fChain) return -5; 
  Long64_t centry = fChain->LoadTree(entry); 
  if (centry < 0) return centry; 
  if (!fChain->InheritsFrom(TChain::Class()))  return centry; 
  TChain *chain = (TChain*)fChain; 
  if (chain->GetTreeNumber() != fCurrent) { 
    fCurrent = chain->GetTreeNumber(); 
  } 
  return centry; 
} 

void RootTuple::Init(TTree* tree){ 
  m_E = 0; 
  m_Px = 0; 
  m_Py = 0; 
  m_Pz = 0; 
  m_barcode = 0; 
  m_iteration = 0; 
  m_weight = 0; 
  m_weightLL = 0; 
  m_weightLR = 0; 
  m_weightRL = 0; 
  m_weightRR = 0; 
  m_weight_ee = 0; 
  m_weight_emu = 0; 
  m_weight_eq = 0; 
  m_weight_qq = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("E", &m_E, &b_E); 
  fChain->SetBranchAddress("Px", &m_Px, &b_Px); 
  fChain->SetBranchAddress("Py", &m_Py, &b_Py); 
  fChain->SetBranchAddress("Pz", &m_Pz, &b_Pz); 
  fChain->SetBranchAddress("barcode", &m_barcode, &b_barcode); 
  fChain->SetBranchAddress("iteration", &m_iteration, &b_iteration); 
  fChain->SetBranchAddress("weight", &m_weight, &b_weight); 
  fChain->SetBranchAddress("weightLL", &m_weightLL, &b_weightLL); 
  fChain->SetBranchAddress("weightLR", &m_weightLR, &b_weightLR); 
  fChain->SetBranchAddress("weightRL", &m_weightRL, &b_weightRL); 
  fChain->SetBranchAddress("weightRR", &m_weightRR, &b_weightRR); 
  fChain->SetBranchAddress("weight_ee", &m_weight_ee, &b_weight_ee); 
  fChain->SetBranchAddress("weight_emu", &m_weight_emu, &b_weight_emu); 
  fChain->SetBranchAddress("weight_eq", &m_weight_eq, &b_weight_eq); 
  fChain->SetBranchAddress("weight_qq", &m_weight_qq, &b_weight_qq); 
} 

