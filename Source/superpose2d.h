#ifndef _SUPERPOSE2D_
#define _SUPERPOSE2D_

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH2.h>
#include <TFile.h>
#include "TLegend.h"
#include "TLatex.h"

TCanvas* superpose2d(TString histName1,        TString histTitle1,          TString fileName1, 
                   TString histName2 = "NULL", TString histTitle2 = "NULL", TString fileName2 = "NULL", 
                   TString histName3 = "NULL", TString histTitle3 = "NULL", TString fileName3 = "NULL");

#endif