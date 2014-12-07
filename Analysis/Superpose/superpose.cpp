// ================================= header ====================================
// superpose.cpp
// ------------------
// Author: Declan Millar
//
// Superimposes 2 histograms (TH1Ds) from 2 different root files. The
// histograms must have the same name in all files.

// usage:
// $ superposeHistos <histName> <filename1.root> <filename2.root> --options

// The executable can be added to $PATH and run in the directory
// containing the root files.
// ================================= start =====================================

// c++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

// root includes
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH1.h>
#include <TFile.h>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // copy argc
  const int args = argc;

  if (args<4)
  {
    std::cout << "usage:  " << argv[0] << "<histName> " << " <fileName1.root> "
    << " <fileName2.root> " << std::endl;
    return 0;
  }

  // compulsory arguments
  std::string histName; // name of histogram to compare
  std::string fileName1; // name of first root file
  std::string fileName2; // name of second root file
  std::string epsFileName; // name of epsoutput file

  // optional arguments
  bool epsOutput;
  bool logY;
  bool normalize;
  Double_t rangeMax;
  Double_t rangeMin;
  std::string fileName3; // name of optional 3rd file

  // set optional arguments
  po::options_description desc("Options for my program");
    desc.add_options()
        ("epsOutput,e", po::value<bool>(& epsOutput)->default_value(false),
            "If true save canvas directly to an eps file.")

        ("logY,l", po::value<bool>(& logY)->default_value(false),
            "If true make y-axis logarithmic scale.")

        ("normalize,n", po::value<bool>(& normalize)->default_value(false),
            "If true normalize histos")

        // both rangeMin and rangeMax must be set simultaneously for any effect.

        ("rangeMin,a", po::value<Double_t>(& rangeMin),
            "The minimum range of the x axis.")

        ("rangeMax,b", po::value<Double_t>(& rangeMax),
            "The maximum range of the x axis.")

        ("optionalFile,o", po::value<std::string>(& fileName3)->default_value("NULL"),
            "Optional third input root file.") 
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

  // local variables
  int nFiles;
  std::string name1;
  std::string name2;
  std::string name3;

  // ---------------------------------------------------------------------------
  if (fileName3 != "NULL") nFiles=3;
  else nFiles =2;

  // print 
  printf("Superimposing %i histograms...\n",nFiles);

  histName  = argv[1];
  fileName1 = argv[2];  
  fileName2 = argv[3];

  // check strings
  // printf("Argument strings: %s, %s, %s, %s\n",histName.c_str(),
  // fileName1.c_str(),fileName2.c_str());

  name1 = histName+'@'+fileName1;
  name2 = histName+'@'+fileName2;
  // name3 = histName+'@'+fileName3;

  // ---------------------------------------------------------------------------
  // Find unique info for each histo
  // we have remove 10 for the _<date> before the . and +4 to remove LHC_
  // note this is very inflexible.
  // int firstIndex;
  // int lastIndex;

  // std::string name1; // name of third root file  
  // firstIndex = fileName1.find_first_of("LHC");
  // lastIndex = fileName1.find_last_of(".");
  // name1 = fileName1.substr(firstIndex+4,lastIndex-(firstIndex+3)-10);

  // std::string name2; // name of third root file
  // firstIndex = fileName2.find_first_of("LHC");
  // lastIndex = fileName2.find_last_of(".");
  // name2 = fileName2.substr(firstIndex+4,lastIndex-(firstIndex+3)-10);

  // std::string name3; // name of third root file
  // firstIndex = fileName3.find_first_of("LHC");
  // lastIndex = fileName3.find_last_of(".");
  // name3 = fileName3.substr(firstIndex+4,lastIndex-(firstIndex+3)-10);
  // ---------------------------------------------------------------------------
  
  // TApplication changes argv so the arguments must have already been copied!
  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  // change autogenerated gStyle
  gStyle->SetOptStat(0);
  gStyle->SetLegendFont(132);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetOptTitle(0);
  
  // create canvas
  TCanvas *canvas = new TCanvas( histName.c_str() ,histName.c_str() );
  canvas->cd();

  // 1st root file
  TFile f1( fileName1.c_str() ,"READ" );
  TH1D *h1 = ( TH1D* )f1.Get( histName.c_str() );
  // h1->SetTitle( fileName1.c_str() );
  h1->SetTitle( "SM only" );
  h1->Draw();
  h1->SetLineColor( kBlue );
  h1->GetYaxis()->SetTitleOffset( 1.3 );
  h1->GetXaxis()->SetTitleOffset( 1.2 );

  // 2nd root file
  TFile f2( fileName2.c_str() ,"READ" );
  TH1D *h2 = ( TH1D* )f2.Get( histName.c_str() );
  // h2->SetTitle( fileName2.c_str() );
  h2->SetTitle( "Z' only" );
  h2->Draw( "SAME" );
  h2->SetLineColor( kRed );

  // normalize histograms
  if ( normalize == true )
  { 
    std::string yTitle;
    yTitle=h1->GetYaxis()->GetTitle();
    yTitle="1/#sigma #times " + yTitle;
    h1->GetYaxis()->SetTitle(yTitle.c_str());
    // printf("yTitle=%s\n",yTitle.c_str());
    h1->Scale( 1.0 / h1->Integral() );
    h2->Scale( 1.0 / h2->Integral() );
    // if ( fileName3 != "NULL" )h3->Scale( 1.0 / h3->Integral() );
  }


  // set range user
  // if( rangeMax == rangeMax && rangeMax == rangeMax )
  // {
  //   h1->GetXaxis()->SetRangeUser( rangeMin ,rangeMax ); 
  //   h2->GetXaxis()->SetRangeUser( rangeMin ,rangeMax );
  //   // if ( fileName3 != "NULL" )h3->GetXaxis()->SetRangeUser( rangeMin ,rangeMax );
  // }
  // printf("%f,%f\n", rangeMin,rangeMax);

  // let logarithmic y axis scale
  if ( logY == true ){ canvas->SetLogy(); }

  // Legend
  canvas->BuildLegend( 0.70 ,0.70 ,0.88 ,0.88 ,"" );
  
  // ROOT app
  if ( epsOutput == false )
  {
    printf("Running ROOT app...\n");
    RootApp->Run( kTRUE );
  }
  else if ( epsOutput == true )
  {
    epsFileName=histName+"_2to6_SMMvsSM_EW.eps";
    printf("Saving to %s\n",epsFileName.c_str());
    canvas->SaveAs(epsFileName.c_str());
  }

  // Wrapping up
  printf("Superposition complete. Have a nice day.\n");
  return 0; 
}
// ================================== end ======================================