#!/usr/bin/env python
import AtlasStyle
import ROOT, sys, optparse, os, glob, subprocess, math
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(False)

def BinsMatch(hist1, hist2):
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        return False
    if hist1.GetMinimum() != hist2.GetMinimum():
        return False
    if hist1.GetMaximum() != hist2.GetMaximum():
        return False
    return True

def PlotSignificance(hist1, hist2):
    if not BinsMatch(hist1, hist2):
      print "Warning: bins do not match."
    name = hist1.GetName() + "_sig"
    hist = hist1.Clone(name)
    hist.Add(hist2, -1)
    for i in range(hist.GetNbinsX()):
        error1 = hist1.GetBinError(i)
        error2 = hist2.GetBinError(i)
        print error1, error2
        error = math.sqrt(error1*error1 + error2*error2)
        hist.SetBinContent(i, hist.GetBinContent(i)/error)
    hist.GetYaxis().SetTitle("Significance")
    labelSize = hist1.GetXaxis().GetLabelSize()
    titleSize = hist1.GetXaxis().GetTitleSize()
    titleOffset = hist1.GetXaxis().GetTitleOffset()
    print "x label size: %f \n" % xLabelSize
    hist.GetXaxis().SetLabelSize(labelSize*3.2)
    hist.GetXaxis().SetTitleSize(titleSize*3.2)
    hist.GetXaxis().SetTitleOffset(titleOffset/1.5)

    hist.GetYaxis().SetLabelSize(labelSize*2)
    hist.GetYaxis().SetTitleSize(titleSize*3.2)
    hist.GetYaxis().SetTitleOffset(titleOffset/3.2)
    # hist.GetYaxis().SetNdivisions(3)
    return hist

usage = "usage: overlay.py hist file [options]"
parser = optparse.OptionParser(usage)

parser.add_option("-n", "--normalise", default = False, action = "store_true" , help = "normalise plots")
parser.add_option("-2","--filename2", default = "", action = "store" , help = "specify second filename")
parser.add_option("-3","--filename3", default = "", action = "store" , help = "specify third filename")
parser.add_option("-4","--filename4", default = "", action = "store" , help = "specify fourth filename")
parser.add_option("-e", "--errors", default = False, action = "store_true" , help = "display errors")
parser.add_option("-r", "--adjust_range", default = False, action = "store" , help = "adjust range")
parser.add_option("-l", "--xmin", type="float", default = -99.9, action = "store" , help = "xmin")
parser.add_option("-u", "--xmax", type="float", default = -99.9, action = "store" , help = "xmax")
parser.add_option("-s", "--significance", default = False, action = "store_true" , help = "plot significance")
parser.add_option("-E", "--eps", default = False, action = "store_true" , help = "save plot as eps")
parser.add_option("-o", "--overlap", default = False, action = "store_true" , help = "find overlapping area")
parser.add_option("--label1", default = "", action = "store" , help = "set first label")
parser.add_option("--label2", default = "", action = "store" , help = "set second label")
parser.add_option("--label3", default = "", action = "store" , help = "set third label")
parser.add_option("--label4", default = "", action = "store" , help = "set fourth label")


(option, args) = parser.parse_args()

if len(args) < 2:
  sys.exit("%s" % usage)

draw_option = "e1x0p" if option.errors else "histsame"
legend = ROOT.TLegend(0.70, 0.70, 0.88, 0.88, "")

# canvas
canvas = ROOT.TCanvas("canvas","canvas", 1200, 800)
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

color1 = ROOT.gROOT.GetColor(20)
color2 = ROOT.gROOT.GetColor(20)

histname = str(args[0])
filename = str(args[1])

filename2 = ""
if option.filename2 != "":
    filename2 = option.filename2
filename3 = ""
if option.filename3 != "":
    filename3 = option.filename3
filename4 = ""
if option.filename4 != "":
    filename4 = option.filename4

if os.path.isfile("%s" % filename) is False:
  sys.exit("%s does not exist" % filename)
file1 = ROOT.TFile(filename, "read")
if not file1.IsOpen():
    print "failed to open %s\n" % filename
try:
    hist = file1.Get(histname)
    hist.Draw(draw_option)
    legend.AddEntry(hist, filename)
    # hist.SetMarkerColor(ROOT.kRed-3)
    hist.SetLineColor(ROOT.kRed-3)
    hist.SetLineStyle(1)

except ReferenceError:
    sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename, histname))

if filename2 != "":
    if os.path.isfile("%s" % filename2) is False:
      sys.exit("%s does not exist" % filename2)
    file2 = ROOT.TFile(filename2, "read")
    if not file2.IsOpen():
        print "failed to open %s\n" % filename2
    try:
        hist2 = file2.Get(histname)
        hist2.Draw(draw_option)
        # hist2.SetMarkerColor(ROOT.kAzure+7)
        hist2.SetLineColor(ROOT.kAzure+7)
        hist.SetLineStyle(2)
        legend.AddEntry(hist2, filename2)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename2, histname))

if filename3 != "":
    if os.path.isfile("%s" % filename3) is False:
      sys.exit("%s does not exist" % filename3)
    file3 = ROOT.TFile(filename3, "read")
    if not file3.IsOpen():
        print "failed to open %s\n" % filename3
    try:
        hist3 = file2.Get(histname)
        hist3.Draw(draw_option)
        # hist3.SetMarkerColor(ROOT.kAzure+7)
        hist3.SetLineColor(ROOT.kAzure+7)
        hist.SetLineStyle(2)
        legend.AddEntry(hist3, filename3)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename3, histname))

if filename4 != "":
    if os.path.isfile("%s" % filename4) is False:
      sys.exit("%s does not exist" % filename4)
    file4 = ROOT.TFile(filename4, "read")
    if not file4.IsOpen():
        print "failed to open %s\n" % filename4
    try:
        hist4 = file2.Get(histname)
        hist4.Draw(draw_option)
        hist4.SetMarkerColor(ROOT.kAzure+7)
        hist4.SetLineColor(ROOT.kAzure+7)
        hist.SetLineStyle(2)
        legend.AddEntry(hist4, filename4)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename4, histname))

legend.SetBorderSize(0)
legend.Draw()

if option.eps:
    canvas.SaveAs("%s_%s.eps" % (histname, filename))

# adjust y-axis range
if option.adjust_range:
    min_value = hist1.GetBinContent(hist1.GetMinimumBin())
    max_value = hist1.GetBinContent(hist1.GetMaximumBin())
    max_value *= 1.2
    if min_value < 0:
        min_value = 1.2
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
#
# normalize histograms
if option.normalise:
    ytitle = hist.GetYaxis().GetTitle()
    ytitle = "1/#sigma #times " + ytitle
    hist.GetYaxis().SetTitle(yTitle)
    hist.Scale(1.0/abs(hist.Integral()))
    if filename2 != "":
        hist2.Scale(1.0/abs(hist2.Integral()))
    if filename2 != "":
        hist3.Scale(1.0/abs(hist3.Integral()))
    if filename2 != "":
        hist4.Scale(1.0/abs(hist4.Integral()))

sigPerOverlap = 0
# find overlapping area (histograms must have the same user ranges and same number of bins)
if option.overlap:
    overlapPerBin = 0
    overlap = 0
    for i in range(hist1.GetNbinsX()):
        bin1 = hist1.GetBin(i)
        bin2 = hist2.GetBin(i)
        one = hist1.GetBinContent(bin1)
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

xmin = option.xmin
xmax = option.xmax
if xmin != -99.9 and xmax != -99.9:
    hist.GetXaxis().SetRangeUser(xmin, xmax)
    if filename2 != "":
        hist2.GetXaxis().SetRangeUser(xmin, xmax)
    if filename3 != "":
        hist3.GetXaxis().SetRangeUser(xmin, xmax)
    if filename4 != "":
        hist4.GetXaxis().SetRangeUser(xmin, xmax)

if option.overlap:
    sigOverlap = "Signal in overlapping area = " + str(sigPerOverlap) + "%%"
    texBox = ROOT.TLatex(0.5,0.5, SigOverlap)
    texBox.Draw()

if option.significance:
    lower_pad.cd()
    if filename2 != "":
        hist2sig.Draw("HIST")
    if filename3 != "":
        h3sig.Draw("HIST SAME")
    if filename4 != "":
        h4sig.Draw("HIST SAME")
    hist1.GetXaxis().SetLabelSize(0)

raw_input()
