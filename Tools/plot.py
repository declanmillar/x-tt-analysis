#!/usr/bin/env python
import ROOT, sys, optparse, os, glob, subprocess, math
import atlas_style

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetHatchesSpacing(0.3);
ROOT.gStyle.SetHatchesLineWidth(1);

ROOT.TGaxis.SetMaxDigits(4)

def BinsMatch(hist, hist2):
    if hist.GetNbinsX() != hist2.GetNbinsX():
        return False
    if hist.GetMinimum() != hist2.GetMinimum():
        return False
    if hist.GetMaximum() != hist2.GetMaximum():
        return False
    return True

def PlotSignificance(hist, hist2):
    if not BinsMatch(hist, hist2):
      print "Warning: bins do not match."
    name = hist.GetName() + "_sig"
    hist = hist.Clone(name)
    hist.Add(hist2, -1)
    for i in range(hist.GetNbinsX()):
        error1 = hist.GetBinError(i)
        error2 = hist2.GetBinError(i)
        print error1, error2
        error = math.sqrt(error1*error1 + error2*error2)
        hist.SetBinContent(i, hist.GetBinContent(i)/error)
    hist.GetYaxis().SetTitle("Significance")
    labelSize = hist.GetXaxis().GetLabelSize()
    titleSize = hist.GetXaxis().GetTitleSize()
    titleOffset = hist.GetXaxis().GetTitleOffset()
    print "x label size: %f \n" % xLabelSize
    hist.GetXaxis().SetLabelSize(labelSize*3.2)
    hist.GetXaxis().SetTitleSize(titleSize*3.2)
    hist.GetXaxis().SetTitleOffset(titleOffset/1.5)
    hist.GetYaxis().SetMaxDigits(2)
    hist.GetYaxis().SetLabelSize(labelSize*2)
    hist.GetYaxis().SetTitleSize(titleSize*3.2)
    hist.GetYaxis().SetTitleOffset(titleOffset/3.2)
    # hist.GetYaxis().SetNdivisions(3)
    return hist

usage = "usage: overlay.py hist file [options]"
parser = optparse.OptionParser(usage)

parser.add_option("-n", "--normalise", default = False, action = "store_true" , help = "normalise plots")
parser.add_option("-2","--f2", default = "", action = "store" , help = "specify second f")
parser.add_option("-3","--f3", default = "", action = "store" , help = "specify third f")
parser.add_option("-4","--f4", default = "", action = "store" , help = "specify fourth filename")
parser.add_option("--h2", default = "", action = "store" , help = "specify second histogram")
parser.add_option("--h3", default = "", action = "store" , help = "specify third histogram")
parser.add_option("--h4", default = "", action = "store" , help = "specify fourth histogram")#
parser.add_option("--l1", default = "", action = "store" , help = "set first label")
parser.add_option("--l2", default = "", action = "store" , help = "set second label")
parser.add_option("--l3", default = "", action = "store" , help = "set third label")
parser.add_option("--l4", default = "", action = "store" , help = "set fourth label")
parser.add_option("--c1", default = "", action = "store" , help = "set first color")
parser.add_option("--c2", default = "", action = "store" , help = "set second color")
parser.add_option("--c3", default = "", action = "store" , help = "set third color")
parser.add_option("--c4", default = "", action = "store" , help = "set fourth color")
parser.add_option("--legend_low", default = False, action = "store_true" , help = "put legend at bottom")
parser.add_option("-e", "--errors", default = False, action = "store_true" , help = "display errors")
parser.add_option("-p", "--pause", default = True, action = "store_false" , help = "pause to show graphs")
parser.add_option("-y", "--adjusty", default = False, action = "store_true" , help = "auto adjust range (some issues)")
parser.add_option("--xmin", type="float", default = -99.9, action = "store" , help = "xmin")
parser.add_option("--xmax", type="float", default = -99.9, action = "store" , help = "xmax")
parser.add_option("--ymin", type="float", default = -99.9, action = "store" , help = "ymin")
parser.add_option("--ymax", type="float", default = -99.9, action = "store" , help = "ymax")
parser.add_option("--ytitle", default = "", action = "store" , help = "ytitle")
parser.add_option("-s", "--significance", default = False, action = "store_true" , help = "plot significance")
parser.add_option("-E", "--eps", default = False, action = "store_true" , help = "save plot as eps")
parser.add_option("-o", "--overlap", default = False, action = "store_true" , help = "find overlapping area")
parser.add_option("-P", "--plot_dir", default = "/Users/declan/Code/declans-research-logbook/plots", action = "store_true" , help = "plot directory")

(option, args) = parser.parse_args()

if len(args) < 2:
  sys.exit("%s" % usage)

draw_option = "e2 hist same" if option.errors else "hist same"
if option.legend_low:
    legend = ROOT.TLegend(0.7, 0.20, 0.9, 0.4, "")
else:
    legend = ROOT.TLegend(0.7, 0.70, 0.9, 0.9, "")

# canvas
canvas = ROOT.TCanvas("canvas","canvas", 1920, 1080)
canvas.cd()

# pad

upper_pad = ROOT.TPad("upper_pad","upper_pad", 0, 0, 1, 1)
upper_pad.SetFillColor(-1)
top_margin = 0.07
right_margin = 0.05
left_margin = 0.1
bottom_margin = 0.12
upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
upper_pad.Draw()
upper_pad.cd()

if option.significance:
    lowerPad = ROOT.TPad("lower_pad", "lower_pad", 0, 0.05, 1, 0.3)
    lowerPad.SetFillColor(-1)
    lowerPad.SetTopMargin(0)
    lowerPad.Draw()
    lowerPad.SetBottomMargin(0.3)
    lowerPad.SetTopMargin(0)
    lowerPad.SetRightMargin(right_margin)
    lowerPad.SetLeftMargin(left_margin)

histname = str(args[0])
filename = str(args[1])

if option.f2 == "":
    filename2 = filename
else:
    filename2 = option.f2

if option.f3 == "":
    filename3 = filename
else:
    filename3 = option.f3

if option.f4 == "":
    filename4 = filename
else:
    filename4 = option.f4

if option.f2 == "":
    filename2 = filename
else:
    filename2 = option.f2

if option.f3 == "":
    filename3 = filename
else:
    filename3 = option.f3

if option.f4 == "":
    filename4 = filename
else:
    filename4 = option.f4

if option.h2 == "":
    histname2 = histname
else:
    histname2 = option.h2

if option.h3 == "":
    histname3 = histname
else:
    histname3 = option.h3

if option.h4 == "":
    histname4 = histname
else:
    histname4 = option.h4

if os.path.isfile("%s" % filename) is False:
  sys.exit("%s does not exist" % filename)
file1 = ROOT.TFile(filename, "read")
if not file1.IsOpen():
    print "failed to open %s\n" % filename
try:
    hist = file1.Get(histname)
    xmin = option.xmin
    xmax = option.xmax
    if xmin != -99.9 and xmax != -99.9:
        hist.GetXaxis().SetRangeUser(xmin, xmax)
    ymin = option.ymin
    ymax = option.ymax
    if ymin != -99.9 and ymax != -99.9:
        hist.GetYaxis().SetRangeUser(ymin, ymax)

    if option.l1 != "":
        labelname1 = option.l1
    else:
        if option.f2 == "" and option.h2 != "":
            labelname1 = hist.GetTitle()
        else:
            labelname1 = filename

    if option.normalise:
        ytitle = hist.GetYaxis().GetTitle()
        # print hist.Integral()
        if ytitle == "Events" or "AFB" in histname:
            ytitle = "Normalised " + ytitle
        else:
            ytitle = "1/#sigma #times " + ytitle
        hist.GetYaxis().SetTitle(ytitle)
        hist.Scale(1.0/abs(hist.Integral()))
    hist.SetLineColor(ROOT.kAzure-7)
    hist.SetMarkerColor(ROOT.kAzure-7)
    hist.SetMarkerStyle(0)
    hist.DrawCopy("h hist")
    hist.SetFillColor(ROOT.kAzure-7)
    hist.SetFillStyle(3354)
    hist.DrawCopy("e2 same")
    legend.AddEntry(hist, labelname1)

except ReferenceError:
    sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename, histname))

if option.f2 != "" or option.h2 != "":
    if os.path.isfile("%s" % filename2) is False:
      sys.exit("%s does not exist" % filename2)
    file2 = ROOT.TFile(filename2, "read")
    if not file2.IsOpen():
        print "failed to open %s\n" % filename2
    try:
        hist2 = file2.Get(histname2)
        if xmin != -99.9 and xmax != -99.9:
            hist2.GetXaxis().SetRangeUser(xmin, xmax)
        if option.l2 != "":
            labelname2 = option.l2
        else:
            if option.f2 == "" and option.h2 != "":
                labelname2 = hist2.GetTitle()
            else:
                labelname2 = filename2
            # print hist2.Integral()
        if option.normalise:
            ytitle2 = hist2.GetYaxis().GetTitle()
            if ytitle2 == "Events" or "AFB" in histname2:
                ytitle2 = "Normalised " + ytitle2
            else:
                ytitle2 = "1/#sigma #times " + ytitle2
            hist2.GetYaxis().SetTitle(ytitle2)
            hist2.Scale(1.0/abs(hist2.Integral()))
        hist2.SetMarkerColor(ROOT.kRed-7)
        hist2.SetLineColor(ROOT.kRed-7)
        hist2.SetMarkerStyle(0)
        hist2.DrawCopy("h hist same")
        hist2.SetFillColor(ROOT.kRed-7)
        hist2.SetFillStyle(3354)
        hist2.DrawCopy("e2 same")
        legend.AddEntry(hist2, labelname2)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename2, histname))

if option.f3 != "" or option.h3 != "":
    if os.path.isfile("%s" % filename3) is False:
      sys.exit("%s does not exist" % filename3)
    file3 = ROOT.TFile(filename3, "read")
    if not file3.IsOpen():
        print "failed to open %s\n" % filename3
    try:
        hist3 = file3.Get(histname3)
        if xmin != -99.9 and xmax != -99.9:
            hist3.GetXaxis().SetRangeUser(xmin, xmax)
        if option.l3 != "":
            labelname3 = option.l3
        else:
            if option.f3 == "" and option.h3 != "":
                labelname3 = hist3.GetTitle()
            else:
                labelname3 = filename3
            # print hist3.Integral()
        if option.normalise:
            ytitle3 = hist3.GetYaxis().GetTitle()
            if ytitle3 == "Events" or "AFB" in histname3:
                ytitle3 = "Normalised " + ytitle3
            else:
                ytitle3 = "1/#sigma #times " + ytitle3
            hist3.GetYaxis().SetTitle(ytitle3)
            hist3.Scale(1.0/abs(hist3.Integral()))
        hist3.SetMarkerColor(ROOT.kSpring-7)
        hist3.SetLineColor(ROOT.kSpring-7)
        hist3.SetMarkerStyle(0)
        hist3.DrawCopy("h hist same")
        hist3.SetFillColor(ROOT.kSpring-7)
        hist3.SetFillStyle(3354)
        hist3.DrawCopy("e2 same")
        legend.AddEntry(hist3, labelname3)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename3, histname))

if option.f4 != "" or option.h4 != "":
    if os.path.isfile("%s" % filename4) is False:
      sys.exit("%s does not exist" % filename4)
    file4 = ROOT.TFile(filename4, "read")
    if not file4.IsOpen():
        print "failed to open %s\n" % filename4
    try:
        hist4 = file2.Get(histname4)
        hist4.Draw(draw_option)
        hist4.SetMarkerColor(ROOT.kViolet-7)
        hist4.SetLineColor(ROOT.kViolet-7)
        hist4.SetMarkerStyle(0)
        # hist.SetLineStyle(2)
        if option.l4 != "":
            labelname4 = option.l4
        else:
            if option.f4 == "" and option.h4 != "":
                labelname4 = hist4.GetTitle()
            else:
                labelname4 = filename4
        legend.AddEntry(hist4, labelname4)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename4, histname))

if option.ytitle != "":
    hist.GetYaxis().SetTitle(option.ytitle)

legend.SetBorderSize(0)
legend.Draw()

# adjust y-axis range
if option.adjusty:
    min_value = hist.GetBinContent(hist.GetMinimumBin())
    max_value = hist.GetBinContent(hist.GetMaximumBin())
    if option.h2 != "":
        if hist2.GetBinContent(hist2.GetMinimumBin()) < min_value:
            min_value = hist2.GetBinContent(hist2.GetMinimumBin())
        if hist2.GetBinContent(hist2.GetMaximumBin()) > max_value:
            max_value = hist2.GetBinContent(hist2.GetMaximumBin())
    if option.h3 != "":
        if hist3.GetBinContent(hist3.GetMinimumBin()) < min_value:
            min_value = hist3.GetBinContent(hist3.GetMinimumBin())
        if hist3.GetBinContent(hist3.GetMaximumBin()) > max_value:
            max_value = hist3.GetBinContent(hist3.GetMaximumBin())
    if option.h4 != "":
        if hist4.GetBinContent(hist4.GetMinimumBin()) < min_value:
            min_value = hist4.GetBinContent(hist4.GetMinimumBin())
        if hist4.GetBinContent(hist4.GetMaximumBin()) > max_value:
            max_value = hist4.GetBinContent(hist4.GetMaximumBin())

    max_value *= 1.2
    if min_value < 0:
        min_value *= 1.2
    elif min_value > 0:
        min_value = min_value - min_value*0.2
    else:
        min_value = 0
    hist.GetYaxis().SetRangeUser(min_value, max_value);

if option.significance:
    if filename2 != "":
        sighist2 = PlotSignificance(hist2, hist)
    if filename3 != "":
        sighist3 = PlotSignificance(hist3, hist)
    if filename4 != "":
        sighist4 = PlotSignificance(hist4, hist)

sigPerOverlap = 0
# find overlapping area (histograms must have the same user ranges and same number of bins)
if option.overlap:
    overlapPerBin = 0
    overlap = 0
    for i in range(hist.GetNbinsX()):
        bin1 = hist.GetBin(i)
        bin2 = hist2.GetBin(i)
        one = hist.GetBinContent(bin1)
        two = hist2.GetBinContent(bin2)
        if one <= two:
            overlapPerBin = one
        elif one > two:
            overlapPerBin = two
        if one <= two:
            overlap += one
        elif one > two:
            overlap += two
        print "Overlap per bin = %f\n" % overlapPerBin

    printf("Overlapping area = %f\n", overlap)
    sigPerOverlap = overlap/hist2.Integral()*100
    printf("Signal in overlapping area: %f%%\n", sigPerOverlap)

ymin = option.ymin
ymax = option.ymax
if ymin != -99.9 and ymax != -99.9:
    hist.GetYaxis().SetRangeUser(ymin, ymax)
    # if filename2 != "":
    #     hist2.GetYaxis().SetRangeUser(ymin, ymax)
    # if filename3 != "":
    #     hist3.GetYaxis().SetRangeUser(ymin, ymax)
    # if filename4 != "":
    #     hist4.GetYaxis().SetRangeUser(ymin, ymax)

if option.overlap:
    sigOverlap = "Signal in overlapping area = " + str(sigPerOverlap) + "%%"
    texBox = ROOT.TLatex(0.5,0.5, SigOverlap)
    texBox.Draw()

xmin = option.xmin
xmax = option.xmax
if xmin != -99.9 and xmax != -99.9:
    hist.GetXaxis().SetRangeUser(xmin, xmax)
    if option.f2 != "":
        hist2.GetXaxis().SetRangeUser(xmin, xmax)
    if option.f3 != "":
        hist3.GetXaxis().SetRangeUser(xmin, xmax)
    if option.f4 != "":
        hist4.GetXaxis().SetRangeUser(xmin, xmax)

if option.significance:
    lower_pad.cd()
    if filename2 != "":
        hist2sig.Draw("HIST")
    if filename3 != "":
        h3sig.Draw("HIST SAME")
    if filename4 != "":
        h4sig.Draw("HIST SAME")
    hist.GetXaxis().SetLabelSize(0)

if option.pause:
    raw_input()

if option.eps:
    canvas.SaveAs("%s/%s_%s.eps" % (option.plot_dir,histname, filename))
