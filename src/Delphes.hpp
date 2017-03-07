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
//  *   File Generated on Mon Mar  6 12:14:49 2017                            * 
//  *                                                                         * 
//  ***************************************************************************/ 

// This class is for accessing the Delphes Ntuples  
// You should not have to edit this file 

#ifndef _NTUPLE_DELPHES_H_ 
#define _NTUPLE_DELPHES_H_

#include <vector> 
using std::vector; 
#include <string> 
using std::string; 
#include <map> 
using std::map; 
#include <iostream> 
using std::cout; 
using std::endl; 

#include "TROOT.h"
#include "TApplication.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"
#include <TChain.h> 
#include <TFile.h> 
#include <TLorentzVector.h>
#include <TRef.h> 
#include <TRefArray.h> 

class Delphes{ 
  public: 
    explicit Delphes(TTree* tree); 
    Delphes(TTree* tree,const bool& isMC,const bool& isAFII); 
    virtual ~Delphes(); 
    Long64_t totalEvents(); 
    Long64_t LoadTree(Long64_t entry); 

    // public inline member functions -- Use these to get access to the TTree variables 
    inline TRefArray  EFlowNeutralHadron() const {b_EFlowNeutralHadron->GetEntry(m_currentEvent);return m_EFlowNeutralHadron;} 
    inline Int_t  EFlowNeutralHadron_size() const {b_EFlowNeutralHadron_size->GetEntry(m_currentEvent);return m_EFlowNeutralHadron_size;} 
    inline TRefArray  EFlowPhoton() const {b_EFlowPhoton->GetEntry(m_currentEvent);return m_EFlowPhoton;} 
    inline Int_t  EFlowPhoton_size() const {b_EFlowPhoton_size->GetEntry(m_currentEvent);return m_EFlowPhoton_size;} 
    inline TRef  EFlowTrack() const {b_EFlowTrack->GetEntry(m_currentEvent);return m_EFlowTrack;} 
    inline Int_t  EFlowTrack_size() const {b_EFlowTrack_size->GetEntry(m_currentEvent);return m_EFlowTrack_size;} 
    inline Float_t  Electron() const {b_Electron->GetEntry(m_currentEvent);return m_Electron;} 
    inline Int_t  Electron_size() const {b_Electron_size->GetEntry(m_currentEvent);return m_Electron_size;} 
    inline Float_t  Event() const {b_Event->GetEntry(m_currentEvent);return m_Event;} 
    inline Int_t  Event_size() const {b_Event_size->GetEntry(m_currentEvent);return m_Event_size;} 
    inline TLorentzVector  GenJet() const {b_GenJet->GetEntry(m_currentEvent);return m_GenJet;} 
    inline Int_t  GenJet_size() const {b_GenJet_size->GetEntry(m_currentEvent);return m_GenJet_size;} 
    inline Float_t  GenMissingET() const {b_GenMissingET->GetEntry(m_currentEvent);return m_GenMissingET;} 
    inline Int_t  GenMissingET_size() const {b_GenMissingET_size->GetEntry(m_currentEvent);return m_GenMissingET_size;} 
    inline TLorentzVector  Jet() const {b_Jet->GetEntry(m_currentEvent);return m_Jet;} 
    inline Int_t  Jet_size() const {b_Jet_size->GetEntry(m_currentEvent);return m_Jet_size;} 
    inline Float_t  MissingET() const {b_MissingET->GetEntry(m_currentEvent);return m_MissingET;} 
    inline Int_t  MissingET_size() const {b_MissingET_size->GetEntry(m_currentEvent);return m_MissingET_size;} 
    inline Float_t  Muon() const {b_Muon->GetEntry(m_currentEvent);return m_Muon;} 
    inline Int_t  Muon_size() const {b_Muon_size->GetEntry(m_currentEvent);return m_Muon_size;} 
    inline Float_t  Particle() const {b_Particle->GetEntry(m_currentEvent);return m_Particle;} 
    inline Int_t  Particle_size() const {b_Particle_size->GetEntry(m_currentEvent);return m_Particle_size;} 
    inline Float_t  Photon() const {b_Photon->GetEntry(m_currentEvent);return m_Photon;} 
    inline Int_t  Photon_size() const {b_Photon_size->GetEntry(m_currentEvent);return m_Photon_size;} 
    inline Float_t  ScalarHT() const {b_ScalarHT->GetEntry(m_currentEvent);return m_ScalarHT;} 
    inline Int_t  ScalarHT_size() const {b_ScalarHT_size->GetEntry(m_currentEvent);return m_ScalarHT_size;} 
    inline TRefArray  Tower() const {b_Tower->GetEntry(m_currentEvent);return m_Tower;} 
    inline Int_t  Tower_size() const {b_Tower_size->GetEntry(m_currentEvent);return m_Tower_size;} 
    inline TRef  Track() const {b_Track->GetEntry(m_currentEvent);return m_Track;} 
    inline Int_t  Track_size() const {b_Track_size->GetEntry(m_currentEvent);return m_Track_size;} 
    inline Float_t  Weight() const {b_Weight->GetEntry(m_currentEvent);return m_Weight;} 
    inline Int_t  Weight_size() const {b_Weight_size->GetEntry(m_currentEvent);return m_Weight_size;} 

    inline Long64_t currentEvent() const {return m_currentEvent;} 

  protected: 
    Int_t    GetEntry(Long64_t entry); 
    void     Init(TTree *tree); 

  private: 
    Delphes(); 
    Delphes(const Delphes& rhs);  
    void operator=(const Delphes& rhs); 

    bool m_isMC; 
    bool m_isAFII; 

    TTree          *fChain; 
    int             fCurrent; 

    Long64_t m_currentEvent; 

    TRefArray  m_EFlowNeutralHadron; 
    Int_t  m_EFlowNeutralHadron_size; 
    TRefArray  m_EFlowPhoton; 
    Int_t  m_EFlowPhoton_size; 
    TRef  m_EFlowTrack; 
    Int_t  m_EFlowTrack_size; 
    Float_t  m_Electron; 
    Int_t  m_Electron_size; 
    Float_t  m_Event; 
    Int_t  m_Event_size; 
    TLorentzVector  m_GenJet; 
    Int_t  m_GenJet_size; 
    Float_t  m_GenMissingET; 
    Int_t  m_GenMissingET_size; 
    TLorentzVector  m_Jet; 
    Int_t  m_Jet_size; 
    Float_t  m_MissingET; 
    Int_t  m_MissingET_size; 
    Float_t  m_Muon; 
    Int_t  m_Muon_size; 
    Float_t  m_Particle; 
    Int_t  m_Particle_size; 
    Float_t  m_Photon; 
    Int_t  m_Photon_size; 
    Float_t  m_ScalarHT; 
    Int_t  m_ScalarHT_size; 
    TRefArray  m_Tower; 
    Int_t  m_Tower_size; 
    TRef  m_Track; 
    Int_t  m_Track_size; 
    Float_t  m_Weight; 
    Int_t  m_Weight_size; 

    TBranch*  b_EFlowNeutralHadron; 
    TBranch*  b_EFlowNeutralHadron_size; 
    TBranch*  b_EFlowPhoton; 
    TBranch*  b_EFlowPhoton_size; 
    TBranch*  b_EFlowTrack; 
    TBranch*  b_EFlowTrack_size; 
    TBranch*  b_Electron; 
    TBranch*  b_Electron_size; 
    TBranch*  b_Event; 
    TBranch*  b_Event_size; 
    TBranch*  b_GenJet; 
    TBranch*  b_GenJet_size; 
    TBranch*  b_GenMissingET; 
    TBranch*  b_GenMissingET_size; 
    TBranch*  b_Jet; 
    TBranch*  b_Jet_size; 
    TBranch*  b_MissingET; 
    TBranch*  b_MissingET_size; 
    TBranch*  b_Muon; 
    TBranch*  b_Muon_size; 
    TBranch*  b_Particle; 
    TBranch*  b_Particle_size; 
    TBranch*  b_Photon; 
    TBranch*  b_Photon_size; 
    TBranch*  b_ScalarHT; 
    TBranch*  b_ScalarHT_size; 
    TBranch*  b_Tower; 
    TBranch*  b_Tower_size; 
    TBranch*  b_Track; 
    TBranch*  b_Track_size; 
    TBranch*  b_Weight; 
    TBranch*  b_Weight_size; 
}; 
#endif 

