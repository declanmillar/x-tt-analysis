root
TFile *file = new TFile("dd-AZX-tt-bbllvv.GLR-R-3.13TeV.CT14LL.root", "update")
RootTuple->SetName("events")
RootTuple->SetTitle("events")
gDirectory->Delete("RootTuple;1")
file->Write()
file->Close()
TTree* process = new TTree("process", "process")
double cross_section = 1.1911627439683250E-003
process->Branch("cross_section", &cross_section)
process->Fill()
file->Write()
file->Close()
