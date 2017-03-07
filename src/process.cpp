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

// This class is for accessing the process Ntuples  
// You should not have to edit this file 

#include "process.hpp" 

process::process(TTree* tree) : 
  m_isMC(false), 
  m_isAFII(false) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

process::process(TTree* tree,const bool& isMC,const bool& isAFII) : 
  m_isMC(isMC), 
  m_isAFII(isAFII) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

process::~process(){ 
  if (!fChain) return; 
} 

int process::GetEntry(Long64_t entry){ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
} 

Long64_t process::totalEvents(){ 
  return fChain->GetEntriesFast(); 
} 

Long64_t process::LoadTree(Long64_t entry){ 
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

void process::Init(TTree* tree){ 
  m_cross_section = 0; 
  m_cross_section_uncertainty = 0; 
  m_ebm1 = 0; 
  m_ebm2 = 0; 
  m_idbm1 = 0; 
  m_idbm2 = 0; 
  m_pdfg1 = 0; 
  m_pdfg2 = 0; 
  m_pdfs1 = 0; 
  m_pdfs2 = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("cross_section", &m_cross_section, &b_cross_section); 
  fChain->SetBranchAddress("cross_section_uncertainty", &m_cross_section_uncertainty, &b_cross_section_uncertainty); 
  fChain->SetBranchAddress("ebm1", &m_ebm1, &b_ebm1); 
  fChain->SetBranchAddress("ebm2", &m_ebm2, &b_ebm2); 
  fChain->SetBranchAddress("idbm1", &m_idbm1, &b_idbm1); 
  fChain->SetBranchAddress("idbm2", &m_idbm2, &b_idbm2); 
  fChain->SetBranchAddress("pdfg1", &m_pdfg1, &b_pdfg1); 
  fChain->SetBranchAddress("pdfg2", &m_pdfg2, &b_pdfg2); 
  fChain->SetBranchAddress("pdfs1", &m_pdfs1, &b_pdfs1); 
  fChain->SetBranchAddress("pdfs2", &m_pdfs2, &b_pdfs2); 
} 

