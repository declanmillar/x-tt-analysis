#!/usr/bin/env python
import AtlasStyle
import ROOT, sys, optparse, os, glob, subprocess
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(False)

usage = "usage: overlay.py hist file [options]"
parser = optparse.OptionParser(usage)

parser.add_option("-n", "--normalise", default = False, action = "store_true" , help = "normalise plots")
parser.add_option("-2","--filename2", default = "", action = "store" , help = "specify second filename")
parser.add_option("-e", "--errors", default = False, action = "store_true" , help = "display errors")
parser.add_option("-s", "--significance", default = False, action = "store_true" , help = "plot significance")

(option, args) = parser.parse_args()

if len(args) < 2:
  sys.exit("%s" % usage)

draw_option = "e1x0p" if option.errors else "histsame"
legend = ROOT.TLegend(0.70, 0.70, 0.88, 0.88, "")

# canvas
canvas = ROOT.TCanvas("canvas","canvas", 1200, 800)
canvas.cd()

# pad
# upper_pad = ROOT.TPad("upper_pad","upper_pad", 0, 0, 1, 1)
#upper_pad.SetFillColor(-1)
# top_margin = 0.1
# right_margin = 0.05
# left_margin = 0.1
# bottom_margin = 0.1
# upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
# upper_pad.Draw()
# upper_pad.cd()

color1 = ROOT.gROOT.GetColor(20);
color2 = ROOT.gROOT.GetColor(20);

histname = str(args[0])
filename = str(args[1])
# hist = ROOT.TH1D()
filename2 = ""
if option.filename2 != "":
    filename2 = option.filename2

if os.path.isfile("%s" % filename) is False:
  sys.exit("%s does not exist" % filename)
file = ROOT.TFile(filename, "read");
if not file.IsOpen():
    print "failed to open %s\n" % filename
# try:
hist = file.Get(histname)
# except ReferenceError:
#     sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename, histnames[0]))
hist.Draw(draw_option)
legend.AddEntry(hist, filename)
hist.SetMarkerColor(1)
hist.SetLineColor(1)

if filename2 != "":
    if os.path.isfile("%s" % filename2) is False:
      sys.exit("%s does not exist" % filename2)
    file = ROOT.TFile(filename2, "read");
    if not file.IsOpen():
        print "failed to open %s\n" % filename2
    try:
        hist2 = file.Get(histname)
        hist2.Draw(draw_option)
        hist2.SetMarkerColor(2)
        hist2.SetLineColor(2)
        legend.AddEntry(hist2, filename2)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename2, histname))
legend.SetBorderSize(0)
legend.Draw()
canvas.SaveAs("%s_%s.eps" % (histname, filename))

#
# // adjust y-axis range
# maxValue *= 1.2;
# if (minValue < 0) minValue *= 1.2;
# else if (minValue > 0) minValue = minValue - minValue*0.2;
# else minValue = 0;
#
# h1->GetYaxis()->SetRangeUser(minValue, maxValue);
#
# if option.significance:
# if (plot2) h2sig = Significance(h2, h1);
# if (plot3) h3sig = Significance(h3, h1);
# if (plot4) h4sig = Significance(h4, h1);
# }
#
# // normalize histograms
if option.normalise:
    ytitle = hist.GetYaxis().GetTitle()
    ytitle = "1/#sigma #times " + ytitle
    hist.GetYaxis().SetTitle(yTitle)
    hist.Scale(1.0/abs(hist.Integral()))
# if (plot2) h2->Scale(std::abs(1.0/h2->Integral()));
# if (plot3) h3->Scale(std::abs(1.0/h3->Integral()));
# if (plot4) h4->Scale(std::abs(1.0/h4->Integral()));

# bool overlap = false;
# double sigPerOverlap = 0;
# // find overlapping area (histograms must have the same user ranges and same number of bins)
# if (overlap == true) {
# if (plot3) {
#   double overlapPerBin = 0;
#   double overlap = 0;
#   for (int i = 0; i < h1->GetNbinsX(); i++) {
#     int bin1 = h1->GetBin(i);
#     int bin2 = h2->GetBin(i);
#     double one = h1->GetBinContent(bin1);
#     double two = h2->GetBinContent(bin2);
#     if (one <= two) overlapPerBin = one;
#     else if (one > two) overlapPerBin = two;
#     if (one <= two) overlap += one;
#     else if (one > two) overlap += two;
#     printf("Overlap per bin = %f\n", overlapPerBin);
#   }
#   printf("Overlapping area = %f\n", overlap);
#   sigPerOverlap = overlap/h2->Integral()*100;
#   printf("Signal in overlapping area: %f%%\n", sigPerOverlap);
# }
# }
#
# // double rangeMin = -999;
# // double rangeMax = -999;
# // rangeMin = 2000;
# // rangeMax = 4000;
# // if(rangeMin != -999 && rangeMax != -999) {
# //   h1->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
# //   if (plot3) h2->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
# //   if (plot3) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
# //   if (plot4) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
# // }
#
# // double yRangeMin = -999;
# // double yRangeMax = -999;
# // yRangeMin = -1;
# // yRangeMax = 1;
# // if(yRangeMin != -999 && yRangeMax != -999) {
# //   h1->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
# //   if (plot3) h2->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
# //   if (plot3) h3->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
# //   if (plot4) h3->GetXaxis()->SetRangeUser(yRangeMin, yRangeMax);
# // }
#
#

# if (overlap) {
# std::ostringstream strs;
# strs << sigPerOverlap;
# std::string str = strs.str();
# TString tstr(str);
# TString SigOverlap = "Signal in overlapping area = " + tstr + "%%";
# TLatex* texBox = new TLatex(0.5,0.5, SigOverlap);
# texBox.Draw();
#
# if (findSignificance) {
# lowerPad->cd();
# if (plot2) h2sig->Draw("HIST");
# if (plot3) h3sig->Draw("HIST SAME");
# if (plot4) h4sig->Draw("HIST SAME");
# h1->GetXaxis()->SetLabelSize(0);
# }




raw_input()
