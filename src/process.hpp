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

#ifndef _NTUPLE_PROCESS_H_ 
#define _NTUPLE_PROCESS_H_ 

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

class process{ 
  public: 
    explicit process(TTree* tree); 
    process(TTree* tree,const bool& isMC,const bool& isAFII); 
    virtual ~process(); 
    Long64_t totalEvents(); 
    Long64_t LoadTree(Long64_t entry); 

    // public inline member functions -- Use these to get access to the TTree variables 
    inline Double_t  cross_section() const {b_cross_section->GetEntry(m_currentEvent);return m_cross_section;} 
    inline Double_t  cross_section_uncertainty() const {b_cross_section_uncertainty->GetEntry(m_currentEvent);return m_cross_section_uncertainty;} 
    inline Double_t  ebm1() const {b_ebm1->GetEntry(m_currentEvent);return m_ebm1;} 
    inline Double_t  ebm2() const {b_ebm2->GetEntry(m_currentEvent);return m_ebm2;} 
    inline Double_t  idbm1() const {b_idbm1->GetEntry(m_currentEvent);return m_idbm1;} 
    inline Double_t  idbm2() const {b_idbm2->GetEntry(m_currentEvent);return m_idbm2;} 
    inline Double_t  pdfg1() const {b_pdfg1->GetEntry(m_currentEvent);return m_pdfg1;} 
    inline Double_t  pdfg2() const {b_pdfg2->GetEntry(m_currentEvent);return m_pdfg2;} 
    inline Double_t  pdfs1() const {b_pdfs1->GetEntry(m_currentEvent);return m_pdfs1;} 
    inline Double_t  pdfs2() const {b_pdfs2->GetEntry(m_currentEvent);return m_pdfs2;} 

    inline Long64_t currentEvent() const {return m_currentEvent;} 

  protected: 
    Int_t    GetEntry(Long64_t entry); 
    void     Init(TTree *tree); 

  private: 
    process(); 
    process(const process& rhs);  
    void operator=(const process& rhs); 

    bool m_isMC; 
    bool m_isAFII; 

    TTree          *fChain; 
    int             fCurrent; 

    Long64_t m_currentEvent; 

    Double_t  m_cross_section; 
    Double_t  m_cross_section_uncertainty; 
    Double_t  m_ebm1; 
    Double_t  m_ebm2; 
    Double_t  m_idbm1; 
    Double_t  m_idbm2; 
    Double_t  m_pdfg1; 
    Double_t  m_pdfg2; 
    Double_t  m_pdfs1; 
    Double_t  m_pdfs2; 

    TBranch*  b_cross_section; 
    TBranch*  b_cross_section_uncertainty; 
    TBranch*  b_ebm1; 
    TBranch*  b_ebm2; 
    TBranch*  b_idbm1; 
    TBranch*  b_idbm2; 
    TBranch*  b_pdfg1; 
    TBranch*  b_pdfg2; 
    TBranch*  b_pdfs1; 
    TBranch*  b_pdfs2; 
}; 
#endif 

