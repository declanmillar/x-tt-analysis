// =============================================================================
// superpose.cpp
// ------------------
// Author: Declan Millar
//
// Plots 1 histogram (TH1D) from 1 root file.

// usage:
// $ superposeHistos <histName> <filename.root> --options

// The executable can be added to $PATH and run in the directory
// containing the root files.
// -----------------------------------------------------------------------------

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
namespace ma = boost::math;

int main(int argc, char *argv[])
{
  // copy argc
  const int args = argc;

  if (args<3)
  {
    std::cout << "usage:  " << argv[0] << "<histName> " << " <fileName.root> "
     << std::endl;
    return 0;
  }

  // compulsory arguments
  std::string histName; // name of histogram to compare
  std::string fileName; // name of first root file
  std::string epsFileName; // name of epsoutput file

  // optional arguments
  bool epsOutput;
  bool logY;
  bool normalize;
  bool radians;
  Double_t rangeMax;
  Double_t rangeMin;

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

        ("radians,r", po::value<bool>(& radians)->default_value(false),
            "If true x axis is normalised")


        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);  

  histName = argv[1];
  fileName = argv[2];  

  // check strings
  // printf("Argument strings: %s, %s, %s, %s\n",histName.c_str(),
  // fileName1.c_str(),fileName2.c_str());
  std::string name;
  name = histName+'@'+fileName;

  // name3 = histName+'@'+fileName3;
  
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

  // root file
  TFile f( fileName.c_str() ,"READ" );
  TH1D *h = ( TH1D* )f.Get( histName.c_str() );
  h->SetTitle( fileName.c_str() );
  h->Draw();
  h->SetLineColor( kBlue );
  h->GetYaxis()->SetTitleOffset( 1.3 );
  h->GetXaxis()->SetTitleOffset( 1.2 );

  // normalize histograms
  if ( normalize == true )
  { 
    std::string yTitle;
    yTitle=h->GetYaxis()->GetTitle();
    yTitle="1/#sigma  " + yTitle;
    h->GetYaxis()->SetTitle(yTitle.c_str());
    // printf("yTitle=%s\n",yTitle.c_str());
    h->Scale( 1.0 / h->Integral() );
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

  // set logarithmic y axis scale
  if ( logY == true ){ canvas->SetLogy(); }

  // set radians on x axis
  // TAxis* xAxis = h->GetXaxis();  
  // const double pi =3.141592653589793238463;
  // float i = 0;
  // while (i*pi <= xAxis->GetXmax()){
  //     float bin_index = xAxis->FindBin(i*pi);
  //     std::string mark;
  //     std::stringstream ss;
  //     ss<<i;
  //     if(i==0)mark="0";
  //     else if(i==0.5)mark="";
  //     else if(i==1)mark="#pi";
  //     else if(i==1.5)mark="#frac{3#pi}{2}";
  //     else mark=ss.str()+" #pi";
  //     xAxis->SetBinLabel(bin_index,mark.c_str());
  //     i+=0.5;
  //     canvas->Modified();
  //     canvas->Update();
  //     printf("%s",mark.c_str());
  // }

  // set radians on x axis
  if (radians == true ){
    TAxis* xAxis = h->GetXaxis();
    const double pi =3.141592653589793238463;
    xAxis->SetBinLabel(xAxis->FindBin(0.0*pi),"0");
    xAxis->SetBinLabel(xAxis->FindBin(0.5*pi),"#frac{#pi}{2}");
    xAxis->SetBinLabel(xAxis->FindBin(1.0*pi),"#pi");
    xAxis->SetBinLabel(xAxis->FindBin(1.5*pi),"#frac{3#pi}{2}");
    xAxis->SetBinLabel(xAxis->FindBin(2*pi-2*pi/1000),"2#pi");
    xAxis->LabelsOption("h");
    xAxis->SetTicks("+");
  }
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