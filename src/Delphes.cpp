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

#include "Delphes.hpp" 

Delphes::Delphes(TTree* tree) : 
  m_isMC(false), 
  m_isAFII(false) 
{ 
  gSystem->Load("libDelphes");
  m_currentEvent = 0; 
  this->Init(tree); 
} 

Delphes::Delphes(TTree* tree,const bool& isMC,const bool& isAFII) : 
  m_isMC(isMC), 
  m_isAFII(isAFII) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

Delphes::~Delphes(){ 
  if (!fChain) return; 
} 

int Delphes::GetEntry(Long64_t entry){ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
} 

Long64_t Delphes::totalEvents(){ 
  return fChain->GetEntriesFast(); 
} 

Long64_t Delphes::LoadTree(Long64_t entry){ 
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

void Delphes::Init(TTree* tree){ 
  m_EFlowNeutralHadron = 0; 
  m_EFlowNeutralHadron_size = 0; 
  m_EFlowPhoton = 0; 
  m_EFlowPhoton_size = 0; 
  m_EFlowTrack = 0; 
  m_EFlowTrack_size = 0; 
  m_Electron = 0; 
  m_Electron_size = 0; 
  m_Event = 0; 
  m_Event_size = 0; 
  m_GenJet = 0; 
  m_GenJet_size = 0; 
  m_GenMissingET = 0; 
  m_GenMissingET_size = 0; 
  m_Jet = 0; 
  m_Jet_size = 0; 
  m_MissingET = 0; 
  m_MissingET_size = 0; 
  m_Muon = 0; 
  m_Muon_size = 0; 
  m_Particle = 0; 
  m_Particle_size = 0; 
  m_Photon = 0; 
  m_Photon_size = 0; 
  m_ScalarHT = 0; 
  m_ScalarHT_size = 0; 
  m_Tower = 0; 
  m_Tower_size = 0; 
  m_Track = 0; 
  m_Track_size = 0; 
  m_Weight = 0; 
  m_Weight_size = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("EFlowNeutralHadron", &m_EFlowNeutralHadron, &b_EFlowNeutralHadron); 
  fChain->SetBranchAddress("EFlowNeutralHadron_size", &m_EFlowNeutralHadron_size, &b_EFlowNeutralHadron_size); 
  fChain->SetBranchAddress("EFlowPhoton", &m_EFlowPhoton, &b_EFlowPhoton); 
  fChain->SetBranchAddress("EFlowPhoton_size", &m_EFlowPhoton_size, &b_EFlowPhoton_size); 
  fChain->SetBranchAddress("EFlowTrack", &m_EFlowTrack, &b_EFlowTrack); 
  fChain->SetBranchAddress("EFlowTrack_size", &m_EFlowTrack_size, &b_EFlowTrack_size); 
  fChain->SetBranchAddress("Electron", &m_Electron, &b_Electron); 
  fChain->SetBranchAddress("Electron_size", &m_Electron_size, &b_Electron_size); 
  fChain->SetBranchAddress("Event", &m_Event, &b_Event); 
  fChain->SetBranchAddress("Event_size", &m_Event_size, &b_Event_size); 
  fChain->SetBranchAddress("GenJet", &m_GenJet, &b_GenJet); 
  fChain->SetBranchAddress("GenJet_size", &m_GenJet_size, &b_GenJet_size); 
  fChain->SetBranchAddress("GenMissingET", &m_GenMissingET, &b_GenMissingET); 
  fChain->SetBranchAddress("GenMissingET_size", &m_GenMissingET_size, &b_GenMissingET_size); 
  fChain->SetBranchAddress("Jet", &m_Jet, &b_Jet); 
  fChain->SetBranchAddress("Jet_size", &m_Jet_size, &b_Jet_size); 
  fChain->SetBranchAddress("MissingET", &m_MissingET, &b_MissingET); 
  fChain->SetBranchAddress("MissingET_size", &m_MissingET_size, &b_MissingET_size); 
  fChain->SetBranchAddress("Muon", &m_Muon, &b_Muon); 
  fChain->SetBranchAddress("Muon_size", &m_Muon_size, &b_Muon_size); 
  fChain->SetBranchAddress("Particle", &m_Particle, &b_Particle); 
  fChain->SetBranchAddress("Particle_size", &m_Particle_size, &b_Particle_size); 
  fChain->SetBranchAddress("Photon", &m_Photon, &b_Photon); 
  fChain->SetBranchAddress("Photon_size", &m_Photon_size, &b_Photon_size); 
  fChain->SetBranchAddress("ScalarHT", &m_ScalarHT, &b_ScalarHT); 
  fChain->SetBranchAddress("ScalarHT_size", &m_ScalarHT_size, &b_ScalarHT_size); 
  fChain->SetBranchAddress("Tower", &m_Tower, &b_Tower); 
  fChain->SetBranchAddress("Tower_size", &m_Tower_size, &b_Tower_size); 
  fChain->SetBranchAddress("Track", &m_Track, &b_Track); 
  fChain->SetBranchAddress("Track_size", &m_Track_size, &b_Track_size); 
  fChain->SetBranchAddress("Weight", &m_Weight, &b_Weight); 
  fChain->SetBranchAddress("Weight_size", &m_Weight_size, &b_Weight_size); 
} 

