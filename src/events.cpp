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
//  *   File Generated on Mon Dec  5 15:31:51 2016                            * 
//  *                                                                         * 
//  ***************************************************************************/ 

// This class is for accessing the events Ntuples  
// You should not have to edit this file 

#include "events.hpp" 

events::events(TTree* tree) : 
  m_isMC(false), 
  m_isAFII(false) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

events::events(TTree* tree,const bool& isMC,const bool& isAFII) : 
  m_isMC(isMC), 
  m_isAFII(isAFII) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

events::~events(){ 
  if (!fChain) return; 
  delete m_E; 
  delete m_Px; 
  delete m_Py; 
  delete m_Pz; 
  delete m_barcode; 
} 

int events::GetEntry(Long64_t entry){ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
} 

Long64_t events::totalEvents(){ 
  return fChain->GetEntriesFast(); 
} 

Long64_t events::LoadTree(Long64_t entry){ 
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

void events::Init(TTree* tree){ 
  m_E = 0; 
  m_Px = 0; 
  m_Py = 0; 
  m_Pz = 0; 
  m_barcode = 0; 
  m_weight = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("E", &m_E, &b_E); 
  fChain->SetBranchAddress("Px", &m_Px, &b_Px); 
  fChain->SetBranchAddress("Py", &m_Py, &b_Py); 
  fChain->SetBranchAddress("Pz", &m_Pz, &b_Pz); 
  fChain->SetBranchAddress("barcode", &m_barcode, &b_barcode); 
  fChain->SetBranchAddress("weight", &m_weight, &b_weight); 
} 

