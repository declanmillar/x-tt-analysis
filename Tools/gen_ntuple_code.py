from ROOT import gROOT,TFile,TBranch,TLeaf,TLeafObject,TObjArray,TTree,gDirectory
import os, sys
from datetime import datetime

class genNtupleCode(object):
  """genNtupleCode takes an ntuple (D3PD) and generates code -  Author John Morris <john.morris@cern.ch>  August 2010"""
  def __init__(self,inputFile):
    self.inputNtuple = inputFile
    self.inputFile = TFile(inputFile)
    self.trees = self.determineTrees()
    self.variables = self.allVariables()
    self.alphabeticalVars = self.sortAllVars()
    self.libpwd = os.environ.get("PWD") + '/lib'
 
  def determineTrees(self):
    keys = gDirectory.GetListOfKeys()
    trees=[]
    for i in keys:      
      # if str(i.GetName()).find('physics') != -1  and str(i.GetName()).find('Meta') == -1:
      if str(i.GetName()).find('RootTuple') != -1:
        trees.append(self.inputFile.Get(i.GetName()))
    return trees      
  
  def allVariables(self):
    variables={}
    for i in self.trees:
      variables[i.GetName()]=self.singleTreeVars(i)
    return variables
          
  def singleTreeVars(self,treeName):
    dict={}
    leaves = TObjArray(treeName.GetListOfLeaves())
    for i in leaves:
      branch = i.GetBranch()
      typeName = str(i.GetTypeName())
      if typeName.find('vect') != -1:
        typeName += '*'
      skip = False
      if branch.GetName().partition('.')[2] == 'second':
        typeName = 'map<Int_t,' + typeName + '>*'
        dict[branch.GetName().partition('.')[0]] = typeName
        skip = True
      if skip == False:
        dict[branch.GetName().partition('.')[0]] = typeName

    return dict
   
  def sortAllVars(self):
    sortedVars={}
    for i in self.trees:
      sortedVars[i.GetName()]=self.sortSingleTree(i.GetName())
    return sortedVars
   
  def sortSingleTree(self,treeName):
    dict = self.variables[treeName]
    return self.sortDict(dict)
   
  def sortDict(self,dict):
    keys = dict.keys()
    keys.sort()
    return keys

  def gnugpl(self,file):
    file.write('//  *************************************************************************** \n')
    file.write('//  *                                                                         * \n')
    file.write('//  *   This program is free software; you can redistribute it and/or modify  * \n')
    file.write('//  *   it under the terms of the GNU General Public License as published by  * \n')
    file.write('//  *   the Free Software Foundation; either version 2 of the License, or     * \n')
    file.write('//  *   (at your option) any later version.                                   * \n')
    file.write('//  *                                                                         * \n')
    file.write('//  *   Author: John Morris (john.morris@cern.ch)                             * \n')
    file.write('//  *           Queen Mary University of London                               * \n')
    file.write('//  *   Editor: Declan Millar (declan.millar@cern.ch)                         * \n')
    file.write('//  *   File Generated on ' + datetime.now().ctime() +'                            * \n')              
    file.write('//  *                                                                         * \n')
    file.write('//  ***************************************************************************/ \n')   
    file.write('\n')
  
  def genAllTreeCode(self):
    for i in self.trees:
      print i.GetName()
      if i.GetName().find('L1CaloDB') == -1 and i.GetName().find('TriggerTree') == -1:
        self.genHeader(i.GetName())
        self.genCpp(i.GetName())
      
  def genHeader(self,outputClass):
    varNames = self.alphabeticalVars[outputClass]
    myVars = self.variables[outputClass]
    # if not os.path.exists('D3PDLibs/'):
    #   os.system('mkdir D3PDLibs')  
    # fileName = 'D3PDLibs/' + outputClass + '.h'
    fileName = outputClass + '.h'
    file = open(fileName,'w')
    self.gnugpl(file)
    file.write('// This class is for accessing the ' + outputClass + ' Ntuples  \n')
    file.write('// You should not have to edit this file \n')
    file.write('\n') 
    file.write('#ifndef _NTUPLE_' + outputClass.upper() + '_H_ \n')
    file.write('#define _NTUPLE_' + outputClass.upper() + '_H_ \n')
    file.write('\n')
    file.write('#include <TROOT.h> \n')
    file.write('#include <TChain.h> \n')
    file.write('#include <TFile.h> \n')
    file.write('#include <vector> \n')
    file.write('using std::vector; \n')
    file.write('#include <string> \n')
    file.write('using std::string; \n')
    file.write('#include <map> \n')
    file.write('using std::map; \n')
    file.write('#include <iostream> \n')
    file.write('using std::cout; \n')
    file.write('using std::endl; \n')
    file.write('\n')
    file.write('class ' + outputClass + '{ \n')
    file.write('  public: \n')
    file.write('    explicit ' + outputClass + '(TTree* tree); \n')
    file.write('    ' + outputClass + '(TTree* tree,const bool& isMC,const bool& isAFII); \n')
    file.write('    virtual ~' + outputClass + '(); \n')
    file.write('    Long64_t totalEvents(); \n')
    file.write('    Long64_t LoadTree(Long64_t entry); \n')
    file.write('\n') 
    file.write('    // public inline member functions -- Use these to get access to the TTree variables \n')
    for i in varNames:
      file.write('    inline ' + myVars[i] + '  ' + i + '() const {b_' + i + '->GetEntry(m_currentEvent);return m_' + i + ';} \n')
    file.write('\n')  
    file.write('    inline Long64_t currentEvent() const {return m_currentEvent;} \n')
    file.write('\n')
    file.write('  protected: \n')
    file.write('    Int_t    GetEntry(Long64_t entry); \n')
    file.write('    void     Init(TTree *tree); \n')
    file.write('\n')  
    file.write('  private: \n')
    file.write('    ' + outputClass +'(); \n')
    file.write('    ' + outputClass + '(const ' + outputClass + '& rhs);  \n')              
    file.write('    void operator=(const ' + outputClass + '& rhs); \n')
    file.write('\n')
    file.write('    bool m_isMC; \n')
    file.write('    bool m_isAFII; \n')
    file.write('\n')    
    file.write('    TTree          *fChain; \n')
    file.write('    int             fCurrent; \n')
    file.write('\n')
    file.write('    Long64_t m_currentEvent; \n')
    file.write('\n')  
    for i in varNames:
      file.write('    ' + myVars[i] + '  m_' + i + '; \n')
    file.write('\n')
    for i in varNames:
      file.write('    TBranch*  b_' + i + '; \n')
    file.write('}; \n')
    file.write('#endif \n')
    file.write('\n')          
    file.close()

       
  def genCpp(self,outputClass):
    varNames = self.alphabeticalVars[outputClass]
    myVars = self.variables[outputClass]
    # if not os.path.exists('D3PDLibs/'):
    #   os.system('mkdir D3PDLibs')
    # fileName = 'D3PDLibs/' + outputClass + '.cpp'
    fileName = outputClass + '.cpp'
    file = open(fileName,'w')
    self.gnugpl(file)
    file.write('// This class is for accessing the ' + outputClass + ' Ntuples  \n')
    file.write('// You should not have to edit this file \n')
    file.write('\n') 
    file.write('#include "' + outputClass + '.h" \n')
    file.write('\n')
    file.write(outputClass + '::' + outputClass + '(TTree* tree) : \n')
    file.write('  m_isMC(false), \n')
    file.write('  m_isAFII(false) \n')
    file.write('{ \n')
    file.write('  m_currentEvent = 0; \n')
    file.write('  this->Init(tree); \n')
    file.write('} \n')
    file.write('\n')
    file.write(outputClass + '::' + outputClass + '(TTree* tree,const bool& isMC,const bool& isAFII) : \n')
    file.write('  m_isMC(isMC), \n')
    file.write('  m_isAFII(isAFII) \n')
    file.write('{ \n')
    file.write('  m_currentEvent = 0; \n')
    file.write('  this->Init(tree); \n')
    file.write('} \n')
    file.write('\n')    
    file.write(outputClass + '::~' + outputClass + '(){ \n') 
    file.write('  if (!fChain) return; \n')
    for i in varNames:
      if myVars[i].find('vector') != -1:
        file.write('  delete m_' + i + '; \n')    
    file.write('} \n')
    file.write('\n') 
    file.write('int ' + outputClass + '::GetEntry(Long64_t entry){ \n')
    file.write('  if (!fChain) return 0; \n')
    file.write('  return fChain->GetEntry(entry); \n')
    file.write('} \n')
    file.write('\n')   
    file.write('Long64_t ' + outputClass + '::totalEvents(){ \n')
    file.write('  return fChain->GetEntriesFast(); \n')
    file.write('} \n')
    file.write('\n')       
    file.write('Long64_t ' + outputClass + '::LoadTree(Long64_t entry){ \n')
    file.write('  m_currentEvent = entry; \n')
    file.write('  if (!fChain) return -5; \n')
    file.write('  Long64_t centry = fChain->LoadTree(entry); \n')
    file.write('  if (centry < 0) return centry; \n')
    file.write('  if (!fChain->InheritsFrom(TChain::Class()))  return centry; \n')
    file.write('  TChain *chain = (TChain*)fChain; \n')
    file.write('  if (chain->GetTreeNumber() != fCurrent) { \n')
    file.write('    fCurrent = chain->GetTreeNumber(); \n')
    file.write('  } \n')
    file.write('  return centry; \n')
    file.write('} \n')
    file.write('\n')  
    file.write('void ' + outputClass + '::Init(TTree* tree){ \n')
    for i in varNames:
      file.write('  m_' + i + ' = 0; \n')
    file.write('\n')     
    file.write('  if (!tree) return; \n')
    file.write('  fChain = tree; \n')
    file.write('  fCurrent = -1; \n')
    file.write('  fChain->SetMakeClass(1); \n')    
    file.write('\n') 
    for i in varNames:
      file.write('  fChain->SetBranchAddress("' + i + '", &m_' + i + ', &b_' + i + '); \n')
    file.write('} \n')
    file.write('\n')               
    file.close()

if __name__ == '__main__':
  try:
    print "Generating access class for", sys.argv[1]
  except IndexError:
    sys.exit("Error: no root file specified.")

  g = genNtupleCode(sys.argv[1])
  g.genAllTreeCode()


    

      
      
      
