#!/usr/bin/env python
import ROOT, sys, optparse, os, glob, subprocess, math
import atlas_style

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetHatchesSpacing(0.3)
ROOT.gStyle.SetHatchesLineWidth(1)
ROOT.TGaxis.SetMaxDigits(4)

usage = "usage: overlay.py hist file [options]"
parser = optparse.OptionParser(usage)

parser.add_option("-n", "--normalise", default = False, action = "store_true" , help = "normalise plots")
parser.add_option("-2","--f2", default = "", action = "store" , help = "specify second filename")
parser.add_option("-3","--f3", default = "", action = "store" , help = "specify third filename")
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
parser.add_option("--z1", default = False, action = "store_true" , help = "set first errors to zero")
parser.add_option("--z2", default = False, action = "store_true" , help = "set second errors to zero")
parser.add_option("--z3", default = False, action = "store_true" , help = "set third errors to zero")
parser.add_option("--z4", default = False, action = "store_true" , help = "set fourth errors to zero")
parser.add_option("--legend_centre", default = False, action = "store_true" , help = "put legend at centre")
parser.add_option("--legend_bottom", default = False, action = "store_true" , help = "put legend at bottom")
parser.add_option("--legend_left", default = False, action = "store_true" , help = "put legend on left")
parser.add_option("--legend_bottom_left", default = False, action = "store_true" , help = "put legend on bottom left")
parser.add_option("-e", "--errors", default = False, action = "store_true" , help = "display errors")
parser.add_option("-p", "--pause", default = True, action = "store_false" , help = "pause to show graphs")
parser.add_option("-y", "--adjusty", default = False, action = "store_true" , help = "auto adjust range (some issues)")
parser.add_option("--xmin", type="float", default = -99.9, action = "store" , help = "xmin")
parser.add_option("--xmax", type="float", default = -99.9, action = "store" , help = "xmax")
parser.add_option("--ymin", type="float", default = -99.9, action = "store" , help = "ymin")
parser.add_option("--ymax", type="float", default = -99.9, action = "store" , help = "ymax")
parser.add_option("--ytitle", default = "", action = "store" , help = "ytitle")
parser.add_option("--xtitle", default = "", action = "store" , help = "xtitle")
parser.add_option("-d", "--distribution", default = True, action = "store_false" , help = "plot distribution")
parser.add_option("-O", "--draw_option", default = "h hist", action = "store" , help = "draw option")
parser.add_option("-s", "--significance", default = False, action = "store_true" , help = "plot significance")
parser.add_option("-C", "--combined", default = False, action = "store_true" , help = "plot combined significance")
parser.add_option("-S", "--significance2", default = False, action = "store_true" , help = "plot significance for h1/h2 & h3/h4")
parser.add_option("-l", "--logy", default = False, action = "store_true" , help = "log y axis")
parser.add_option("-E", "--eps", default = False, action = "store_true" , help = "save plot as eps")
parser.add_option("-D", "--pdf", default = False, action = "store_true" , help = "save plot as pdf")
parser.add_option("-o", "--overlap", default = False, action = "store_true" , help = "find overlapping area")
parser.add_option("-t", "--tag", default = "", action = "store" , help = "add tag to output")
parser.add_option("-P", "--plot_dir", default = "/Users/declan/Code/declans-research-logbook/plots", action = "store" , help = "plot directory")
parser.add_option("-T", "--texbox", default = "", action = "store" , help = "add tag texbox")

(option, args) = parser.parse_args()

def BinsMatch(hist, hist2):
    if hist.GetNbinsX() != hist2.GetNbinsX():
        return False
    if hist.GetMinimum() != hist2.GetMinimum():
        return False
    if hist.GetMaximum() != hist2.GetMaximum():
        return False
    return True

def ZeroErrors(hist):
    for i in range(1,sighist.GetNbinsX()+1):
        hist.SetBinError(i, 0)
    return

def PlotSignificance(hist, hist2):
    if not BinsMatch(hist, hist2):
      print "Warning: bins do not match."
    name = hist.GetName() + "_sig"
    sighist = hist.Clone(name)
    for i in range(1,sighist.GetNbinsX()+1):
        if i == 0:
            continue
        n = hist.GetBinContent(i)
        n2 = hist2.GetBinContent(i)
        error1 = hist.GetBinError(i)
        error2 = hist2.GetBinError(i)
        error = math.sqrt(error1*error1 + error2*error2)
        sighist.SetBinContent(i, abs(n-n2)/error)
    sighist.GetYaxis().SetTitle("Significance")
    labelSize = sighist.GetXaxis().GetLabelSize()
    titleSize = sighist.GetXaxis().GetTitleSize()
    titleOffset = sighist.GetXaxis().GetTitleOffset()
    sighist.GetXaxis().SetLabelSize(labelSize*3.2)
    sighist.GetXaxis().SetTitleSize(titleSize*3.2)
    sighist.GetXaxis().SetTitleOffset(titleOffset/1.5)
    # sighist.GetYaxis().SetMaxDigits(2)
    sighist.GetYaxis().SetLabelSize(labelSize*2)
    sighist.GetYaxis().SetTitleSize(titleSize*3.2)
    sighist.GetYaxis().SetTitleOffset(titleOffset/3.2)
    # hist.GetYaxis().SetNdivisions(3)
    return sighist

fillstyle = 1001

if len(args) < 2:
    sys.exit("%s" % usage)

if (option.significance2):
    option.significance = True
    option.distribution = False
    option.z2 = True
    option.z4 = True

ROOT.gStyle.SetGridStyle(1)

# colors
color4 = ROOT.TColor.GetColor(250.0/255.0, 0.0/255.0, 0.0/255.0)
color2 = 46 #ROOT.TColor.GetColor(0.0/255.0, 90.0/255.0, 130.0/255.0)
color3 = ROOT.TColor.GetColor(0.0/255.0, 130.0/255.0, 90.0/255.0)
color1 = ROOT.TColor.GetColor(64.0/255.0, 64.0/255.0, 64.0/255.0)

draw_option = option.draw_option

if option.errors: draw_option = "e2 " + draw_option
# Create Legend
if option.legend_bottom:
    legend = ROOT.TLegend(0.7, 0.2, 0.9, 0.4, "")
elif option.legend_left:
    legend = ROOT.TLegend(0.15, 0.5, 0.35, 0.65, "")
elif option.legend_bottom_left:
    legend = ROOT.TLegend(0.15, 0.2, 0.35, 0.4, "")
elif option.legend_centre:
    legend = ROOT.TLegend(0.45, 0.7, 0.65, 0.9, "")
else:
    legend = ROOT.TLegend(0.6, 0.65, 0.9, 0.9, "")

# canvas
canvasx = 1920
canvasy = 1080
canvas = ROOT.TCanvas("canvas","canvas", canvasx, canvasy)
canvas.cd()

ROOT.gStyle.SetPalette(57)
ROOT.gStyle.SetNumberContours(999)

# pad
if option.significance and option.distribution:
    upper_pad = ROOT.TPad("upper_pad","upper_pad", 0, 0.16, 1, 1)
    upper_pad.Draw()
    lower_pad = ROOT.TPad("lower_pad", "lower_pad", 0, 0, 1, 0.25)
    lower_pad.Draw()
else:
    upper_pad = ROOT.TPad("upper_pad","upper_pad", 0, 0, 1, 1)
upper_pad.SetFillColor(-1)
if (option.logy):
    upper_pad.SetLogz()
if draw_option == "COLZ": upper_pad.SetLogz()

top_margin = 0.025
right_margin = 0.01
left_margin = 0.1
bottom_margin = 0.13
upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
upper_pad.Draw()
upper_pad.cd()

upper_pad.SetTheta(20) # default: 30
upper_pad.SetPhi(-30) # default: 30
upper_pad.Update()


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
    # hist.Scale(0.2)
    if option.significance:
        hist.GetXaxis().SetLabelSize(0)
        # hist.GetXaxis().SetTickLength(0)
        hist.GetXaxis().SetLabelOffset(999)
        hist.GetXaxis().SetTitleOffset(999)
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
    if option.z1:
        for i in range(1,hist.GetNbinsX()+1):
            hist.SetBinError(i,0)
    if option.normalise:
        ytitle = hist.GetYaxis().GetTitle()
        # print hist.Integral()
        if ytitle == "Events" or "AFB" in histname:
            ytitle = "Normalised " + ytitle
        else:
            ytitle = "1/#sigma #times " + ytitle
        hist.GetYaxis().SetTitle(ytitle)
        hist.Scale(1.0/abs(hist.Integral()))
    # hist.SetError(0)
    # hist.SetFillColor(color1)
    hist.SetLineColor(color1)
    hist.SetMarkerColor(color1)
    hist.SetMarkerStyle(0)
    hist.GetXaxis().SetLabelColor(color1)
    hist.GetYaxis().SetLabelColor(color1)
    hist.GetXaxis().SetLabelSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetTitleSize(0.06)
    # hist.GetXaxis().SetNdivisions(8)
    # hist.GetYaxis().SetNdivisions(2)
    # hist.GetYaxis().SetTitle('cos#theta* (reco)')
    # hist.GetXaxis().CenterTitle()
    # hist.GetYaxis().CenterTitle()
    # hist.GetXaxis().SetTitleOffset(1.4)
    hist.GetYaxis().SetTitleOffset(0.8)
    # hist.GetZaxis().SetTitleOffset(1.3)

    dummy = hist.Clone()
    dummy.SetLineColor(0)
    dummy.SetFillStyle(0)
    dummy.SetFillColor(0)
    dummy.SetLineStyle(0)
    dummy.SetMarkerColor(0)
    dummy.SetMarkerStyle(0)
    dummy.SetTitle("")
    # legend.AddEntry(dummy,"")


    # temporary!
    if option.ytitle != "":
        hist.GetYaxis().SetTitle(option.ytitle)
        # hist.GetYaxis().SetTitle("Expected events with L = 300fb^{-1}")

    if option.xtitle != "":
        hist.GetXaxis().SetTitle(option.xtitle)
        # hist.GetYaxis().SetTitle("Expected events with L = 300fb^{-1}")


    hist.DrawCopy(draw_option)
    if not option.z1:
        #  hist.SetFillColor(color1)
         hist.SetFillColorAlpha(color1, 0.2)
         hist.SetFillStyle(fillstyle)
    if not "COLZ" in  draw_option:
        if not "LEGO" in draw_option:
            hist.DrawCopy("e2 same")
    if not (option.significance2):
        legend.AddEntry(hist, labelname1,"lf")


except:
    sys.exit("Error: check %s contains histogram '%s'" % (filename, histname))

if option.f2 != "" or option.h2 != "":
    if os.path.isfile("%s" % filename2) is False:
      sys.exit("%s does not exist" % filename2)
    file2 = ROOT.TFile(filename2, "read")
    if not file2.IsOpen():
        print "failed to open %s\n" % filename2
    try:
        hist2 = file2.Get(histname2)
        # hist2.Scale(0.2)
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
        if option.z2:
            for i in range(1,hist2.GetNbinsX()+1):
                hist2.SetBinError(i,0)
        if option.normalise:
            ytitle2 = hist2.GetYaxis().GetTitle()
            if ytitle2 == "Events" or "AFB" in histname2:
                ytitle2 = "Normalised " + ytitle2
            else:
                ytitle2 = "1/#sigma #times " + ytitle2
            hist2.GetYaxis().SetTitle(ytitle2)
            hist2.Scale(1.0/abs(hist2.Integral()))
        hist2.SetMarkerColor(color2)
        hist2.SetLineColor(color2)
        hist2.SetMarkerStyle(0)
        hist2.DrawCopy(draw_option + " same")
        if not option.z1:
            hist2.SetFillColorAlpha(color2, 0.2)
            # hist2.SetFillColor(color2)
            # hist2.SetFillStyle(3354)
            hist2.SetFillStyle(fillstyle)
        hist2.DrawCopy("e2 same")
        if not (option.significance2):
            legend.AddEntry(hist2, labelname2,"lf")
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
        # hist3.Scale(0.2)
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
        if option.z3:
            for i in range(1,hist3.GetNbinsX()+1):
                hist3.SetBinError(i,0)
        if option.normalise:
            ytitle3 = hist3.GetYaxis().GetTitle()
            if ytitle3 == "Events" or "AFB" in histname3:
                ytitle3 = "Normalised " + ytitle3
            else:
                ytitle3 = "1/#sigma #times " + ytitle3
            hist3.GetYaxis().SetTitle(ytitle3)
            hist3.Scale(1.0/abs(hist3.Integral()))
        hist3.SetMarkerColor(color3)
        hist3.SetLineColor(color3)
        hist3.SetMarkerStyle(0)
        hist3.DrawCopy(draw_option + " same")
        if not option.z3:
            # hist3.SetFillColor(color3)
            hist3.SetFillColorAlpha(color3, 0.2)
            hist3.SetFillStyle(fillstyle)
        else:
            hist3.SetFillStyle(1)
        hist3.DrawCopy("e2 same")
        if not (option.significance2):
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
        hist4 = file4.Get(histname4)
        if xmin != -99.9 and xmax != -99.9:
            hist4.GetXaxis().SetRangeUser(xmin, xmax)
        if option.l4 != "":
            labelname4 = option.l4
        else:
            if option.f4 == "" and option.h4 != "":
                labelname4 = hist4.GetTitle()
            else:
                labelname4 = filename4
            # print hist4.Integral()
        if option.z4:
            for i in range(1,hist4.GetNbinsX()+1):
                hist4.SetBinError(i,0)
        if option.normalise:
            ytitle4 = hist4.GetYaxis().GetTitle()
            if ytitle4 == "Events" or "AFB" in histname4:
                ytitle4 = "Normalised " + ytitle4
            else:
                ytitle4 = "1/#sigma #times " + ytitle4
            hist4.GetYaxis().SetTitle(ytitle4)
            hist4.Scale(1.0/abs(hist4.Integral()))
        hist4.SetMarkerColor(color4)
        hist4.SetLineColor(color4)
        hist4.SetMarkerStyle(0)
        hist4.DrawCopy(draw_option + " same")
        if not option.z4:
            hist4.SetFillColorAlpha(color4, 0.2)
            # hist4.SetFillColor(color4)
            hist4.SetFillStyle(fillstyle)
        else:
            hist4.SetFillStyle(1)
        hist4.DrawCopy("e2 same")
        if not (option.significance2):
            legend.AddEntry(hist4, labelname4)
    except ReferenceError:
        sys.exit("ReferenceError: check %s contains histogram '%s'" % (filename4, histname))

if option.ytitle != "":
    hist.GetYaxis().SetTitle(option.ytitle)

legend.SetBorderSize(0)
if not "COLZ" in  draw_option:
    if not "LEGO" in draw_option:
        # legend.AddEntry(dummy, "#bf{Model: GLR-R}")
        # legend.AddEntry(dummy, "#bf{#sqrt{s} = 13 TeV}")
        # legend.AddEntry(dummy, "#bf{#int L dt = 100 fb^{-1}}")
        t1 = ROOT.TLatex(0.15, 0.9, "#bf{Model: GLR-R}")
        t1.SetNDC(True)
        t1.Draw()
        t2 = ROOT.TLatex(0.15, 0.8, "#bf{#sqrt{s} = 13 TeV}")
        t2.SetNDC(True)
        t2.Draw()
        t3 = ROOT.TLatex(0.15, 0.7, "#bf{#int L dt = 100 fb^{-1}}")
        t3.SetNDC(True)
        t3.Draw()
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
    hist.GetYaxis().SetRangeUser(min_value, max_value)

if option.significance:
    if filename2 != "":
        sighist2 = PlotSignificance(hist, hist2)
    if filename3 != "":
        sighist3 = PlotSignificance(hist3, hist4)
    # if filename4 != "":
    #     sighist4 = PlotSignificance(hist4, hist)
    if option.distribution:
        sighist2.GetXaxis().SetLabelSize(0.15)
        sighist2.GetXaxis().SetTickLength(0.1)
        sighist2.GetXaxis().SetLabelOffset(0.01)
        sighist2.GetXaxis().SetTitleOffset(0.9)
        sighist2.GetYaxis().SetLabelSize(0.12)
        sighist2.GetYaxis().SetLabelOffset(0.01)
        sighist2.GetYaxis().SetTitleOffset(0.3)
    else:
        sighist2.GetXaxis().SetLabelSize(0.05)
        sighist2.GetXaxis().SetTickLength(0.1)
        sighist2.GetXaxis().SetLabelOffset(0.01)
        sighist2.GetXaxis().SetTitleOffset(1.0)
        sighist2.GetXaxis().SetTitleSize(0.05)
        sighist2.GetYaxis().SetLabelSize(0.05)
        sighist2.GetYaxis().SetLabelOffset(0.01)
        sighist2.GetYaxis().SetTitleOffset(0.8)
        sighist2.GetYaxis().SetTitleSize(0.05)
sigPerOverlap = 0
# find overlapping area (histograms must have the same user ranges and same number of bins)
if option.overlap:
    overlapPerBin = 0
    overlap = 0
    for i in range(1,hist.GetNbinsX()+1):
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



canvas.Draw()


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
    if option.distribution:
        # lower_pad = ROOT.TPad("lower_pad", "lower_pad", 0, 0.05, 1, 0.3)
        lower_pad.SetFillColor(-1)
        lower_pad.SetTopMargin(0)
        lower_pad.SetBottomMargin(0.3)
        lower_pad.SetTopMargin(0)
        lower_pad.SetRightMargin(right_margin)
        lower_pad.SetLeftMargin(left_margin)
        lower_pad.cd()
    combhist = hist.Clone()
    for i in range(1,sighist2.GetNbinsX()+1):
        if i == 0:
            continue
        s2 = sighist2.GetBinContent(i)
        s3 = sighist3.GetBinContent(i)
        comb = math.sqrt(s3*s3 + s2*s2)
        combhist.SetBinContent(i, comb)
    FuckYouTLatexYouCunt = hist.Clone()
    FuckYouTLatexYouCunt.SetLineColor(ROOT.kWhite)
    FuckYouTLatexYouCunt.SetFillStyle(0)
    FuckYouTLatexYouCunt.SetFillColor(0)
    FuckYouTLatexYouCunt.SetLineStyle(0)
    FuckYouTLatexYouCunt.SetMarkerColor(0)
    FuckYouTLatexYouCunt.SetMarkerStyle(0)
    legend.AddEntry(FuckYouTLatexYouCunt,"Model: GLR-R")
    legend.AddEntry(FuckYouTLatexYouCunt," ")
    if option.combined:
        combhist.SetLineColor(ROOT.kOrange+7)
        combhist.SetLineStyle(2)
        combhist.SetMarkerColor(0)
        combhist.SetMarkerStyle(0)
        # combhist.SetFillColor(color1)
        combhist.SetFillStyle(0)
        combhist.SetFillColor(0)
        combhist.GetXaxis().SetLabelSize(0.05)
        combhist.GetXaxis().SetTickLength(0.1)
        combhist.GetXaxis().SetLabelOffset(0.01)
        combhist.GetXaxis().SetTitleOffset(1.0)
        combhist.GetXaxis().SetTitleSize(0.05)
        combhist.GetYaxis().SetLabelSize(0.05)
        combhist.GetYaxis().SetLabelOffset(0.01)
        combhist.GetYaxis().SetTitleOffset(0.8)
        combhist.GetYaxis().SetTitleSize(0.05)
        combhist.SetTitle("Combined")
        # legend.AddEntry(combhist, "Combined")
        combhist.Draw("HIST")
    if filename2 != "":
        sighist2.SetFillStyle(1001)
        sighist2.Draw("HIST")
        legend.AddEntry(sighist2, labelname1)
        if option.ytitle != "":
            sighist2.GetYaxis().SetTitle(option.ytitle)
    if filename3 != "":
        # sighist3.SetFillStyle(3354)
        sighist3.Draw("HIST SAME")
        legend.AddEntry(sighist3, labelname2)

    if option.combined:
        legend.AddEntry(combhist, "Combined")

if not "COLZ" in draw_option:
    if not "LEGO" in draw_option:
        legend.Draw()

if option.pause: raw_input()

siglabel = ""
if option.distribution is False and option.significance is True:
    siglabel += "_sig"

outfilename = filename
outfilename = option.plot_dir + "/" + histname + siglabel + "_" + os.path.splitext(filename)[0] + option.tag
print outfilename

if option.eps:
    canvas.SaveAs("%s.eps" % outfilename)

if option.pdf:
    canvas.SaveAs("%s.pdf" % outfilename)
