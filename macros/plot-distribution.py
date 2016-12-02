#!/usr/bin/env python
# Declan Millar (declan.millar@cern.ch)
import ROOT, sys, optparse, os, glob, subprocess, math

class HistPainter():
    "Makes a canvas. Has members to add histograms and text boxes."
    canvas = ROOT.TCanvas("canvas", "canvas")
    grey = ROOT.TColor.GetColor(64.0/255.0, 64.0/255.0, 64.0/255.0)


    def __init__(self, xwidth, ywidth):

        self.members = []
        self.histograms = []
        self.setlogy = False
        self.setlogz = False
        self.is_th2d = False
        self.ymin = -99.9
        self.ymax = -99.9
        self.xmin = -99.9
        self.xmax = -99.9
        self.ztitle = ""
        self.ytitle = ""
        self.xtitle = ""

        HistPainter.canvas.SetCanvasSize(xwidth, ywidth)
        HistPainter.canvas.Draw()
        HistPainter.canvas.cd()

        return None


    def AddHistogram(self, histname, filename, label, color, style = 1):

        # check file exists
        if os.path.isfile("%s" % filename) is False:
          sys.exit("%s does not exist." % filename)

        # open file
        f = ROOT.TFile(filename, "read")
        if not f.IsOpen():
            sys.exit("failed to open %s\n" % filename)

        # get histogram
        try:
            hist = f.Get(histname)
            hist.SetDirectory(0) # "detach" the histogram from the file
            self.histograms.append(hist)
        except:
            sys.exit("Error: check %s contains %s" % (filename, histname))

        # set range
        if self.ymin != -99.9 and self.ymax != -99.9:
            hist.GetYaxis().SetRangeUser(self.ymin, self.ymax)

        if self.xmin != -99.9 and self.xmax != -99.9:
            hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)

        # axis style
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetTitleSize(0.06)
        hist.GetXaxis().SetTitleOffset(1.05)
        hist.GetYaxis().SetTitleOffset(0.9)

        # hist style
        if not label == "":
            hist.SetTitle(label)
        if not self.ytitle == "":
            hist.GetYaxis().SetTitle(self.ytitle)
        if not self.xtitle == "":
            hist.GetXaxis().SetTitle(self.xtitle)
        if not self.ztitle == "":
            hist.GetZaxis().SetTitle(self.ztitle)
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(0)
        hist.SetLineStyle(style)

        # 2D
        if self.is_th2d:
            hist.SetLineStyle(30)
            hist.SetLineWidth(1)
            hist.GetXaxis().CenterTitle()
            hist.GetYaxis().CenterTitle()
            hist.GetZaxis().CenterTitle()
            hist.GetYaxis().SetNdivisions(4)
            hist.GetXaxis().SetTitleSize(0.06)
            hist.GetYaxis().SetTitleSize(0.06)
            hist.GetZaxis().SetTitleSize(0.06)
            hist.GetXaxis().SetTitleOffset(1.2)
            # hist.GetYaxis().SetTitleOffset(1.2)
            hist.GetYaxis().SetTitleOffset(0.7)
            hist.GetZaxis().SetTitleOffset(0.8)
            hist.SetMinimum(hist.GetBinContent(hist.GetMinimumBin()))
            # hist.Draw("lego2 same fb")
            hist.Draw("colz fb")
        else:
            hist.Draw("hist same")

            # clone to draw errors
            clone = hist.Clone("")
            self.members.append(clone)

            # Uncertainty style
            clone.SetFillColorAlpha(color, 0.2)
            clone.SetFillStyle(1001)
            clone.DrawCopy("e2 same")


        f.Close()


    def GetHistType(self, histname, filename):

        # check file exists
        if os.path.isfile("%s" % filename) is False:
          sys.exit("%s does not exist." % filename)

        # open file
        f = ROOT.TFile(filename, "read")
        if not f.IsOpen():
            sys.exit("failed to open %s\n" % filename)

        # get histogram
        try:
            hist = f.Get(histname)
            if "TH2D" in type(hist).__name__:
                self.is_th2d = True
            hist.SetDirectory(0) # "detach" the histogram from the file
        except:
            sys.exit("Error: check %s contains %s" % (filename, histname))

        if "TH2D" in type(hist).__name__:
            self.is_th2d = True

        f.Close()


    def AddInfoBox(self, model, x, y):

        t1 = ROOT.TLatex(x + 0.02, y - 0.075, "#bf{#it{m_{Z'}} = 3 TeV}")
        t1.SetTextAlign(11)
        t1.SetTextSize(0.04)
        t1.SetNDC(True)
        t1.Draw("same")
        self.members.append(t1)

        # t2 = ROOT.TLatex(x, y, "#bf{#int #it{L dt} = 100 fb^{-1}}, #bf{#it{#sqrt{s}} = 13 TeV}")
        t2 = ROOT.TLatex(x + 0.02, y, "#bf{#it{#sqrt{s}} = 13 TeV}")
        t2.SetNDC(True)
        t2.SetTextSize(0.04)
        t2.Draw("same")
        self.members.append(t2)

        t3 = ROOT.TLatex(x + 0.02, y - 0.15, "#bf{Model: %s}" % model)
        t3.SetNDC(True)
        t3.SetTextSize(0.04)
        t3.Draw("same")
        self.members.append(t3)


    def AddLegend(self, xlow, ylow, xup, yup):
        legend = ROOT.TLegend(xlow, ylow, xup, yup, "")
        legend.SetBorderSize(0)
        for h in self.histograms:
            legend.AddEntry(h, h.GetTitle(), "lf")
        legend.Draw()
        self.members.append(legend)


    def AddPads(self):
        upper_pad = ROOT.TPad("upper_pad", "upper_pad", 0, 0, 1, 1)

        if (self.is_th2d):
            top_margin = 0.05
            left_margin = 0.09
            right_margin = 0.15
            bottom_margin = 0.15
        else:
            top_margin = 0.05
            right_margin = 0.05
            left_margin = 0.12
            bottom_margin = 0.15

        upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
        upper_pad.Draw()
        upper_pad.cd()

        # print "hist type = ", type(h).__name__

        # if "TH2D" in type(h).__name__:
        #     upper_pad.SetTheta(20) # default: 30
        #     # upper_pad.SetPhi(-30) # default: 30
        #     upper_pad.SetPhi(-150)
        #     upper_pad.Update()

        upper_pad.SetGridx()
        upper_pad.SetGridy()


        if (self.setlogy):
            upper_pad.SetLogy()
            print "Setting log y"
        if (self.setlogz):
            upper_pad.SetLogz()
            print "Setting log z"
        self.members.append(upper_pad)


    def SetRange(self, ymin, ymax):
        self.ymin = ymin
        self.ymax = ymax


    def SetXtitle(self, xtitle):
        self.xtitle = xtitle


    def SetYtitle(self, ytitle):
        self.ytitle = ytitle


    def SetZtitle(self, ztitle):
        self.ztitle = ztitle


    def SetHistTitle(self, i, title):
        self.histograms[i].SetTitle("#bf{%s}" % title)


    def SetDomain(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax


    def SetLogy(self):
        self.setlogy = True
    def SetLogz(self):
        self.setlogz = True


    def Save(self, name):
        HistPainter.canvas.SaveAs(name)


    def SetStyle(self):
        import atlas_style


red = ROOT.TColor.GetColor(250.0 / 255.0, 0.0 / 255.0, 0.0 / 255.0)
blue = ROOT.TColor.GetColor(0.0 / 255.0, 90.0 / 255.0, 130.0 / 255.0)
green = ROOT.TColor.GetColor(0.0 / 255.0, 130.0 / 255.0, 90.0 / 255.0)
grey = ROOT.TColor.GetColor(64.0 / 255.0, 64.0 / 255.0, 64.0 / 255.0)
black = ROOT.TColor.GetColor(1.0, 1.0, 1.0)

hs = ["m_t", "m_tbar", "m_tt", "m_tt_perf", "pT_t", "pT_tbar", "pT_t_perf", "pT_tbar_perf", "E_t", "E_tbar", "costheta_tt", "costheta_tt_perf", "m_tt_pT_t_perf", "m_tt_costheta_tt_perf"]#, "deltaR_tt"]

f = "GLR-R-3_dd-AZX-tt-bbllvv.a.root"

i = 0
for h in hs:
    art = HistPainter(1920, 1080)
    # if (h == "m_tt" or h == "E_t"):

    art.GetHistType(h, f)
    
    if ("m_tt_costheta_tt" in h or "m_tt_pT_t" in h):
        art.SetLogz()
    else:
        art.SetLogy()
    art.SetStyle()
    art.AddPads()
    # art.SetXtitle("#it{m_{tt}} [TeV]")
    # art.SetYtitle("d#it{#sigma} / d#it{m_{tt}} [pb/TeV]")
    # art.SetYtitle("A^{*}_{FB}")
    # art.SetZtitle("#it{cos#theta_{l}}")

    art.AddHistogram(h, f, "#bf{SM}", blue)
    if not ("perf" in h): art.AddHistogram(h + "_R", f, "#bf{SM}", red, 2)

    if ("m_tt_" in h): 
        pass
    elif ("perf" in h):
        art.AddInfoBox("SM", 0.2, 0.85)
    else:
        if ("costheta" in h): 
            art.AddInfoBox("SM", 0.45, 0.85)
        else:
            art.AddInfoBox("SM", 0.75, 0.85)

    art.SetHistTitle(0, "generated")
    if not ("perf" in h): 
        art.SetHistTitle(1, "reconstructed")
        if ("costheta" in h): 
            art.AddLegend(0.46, 0.45, 0.62, 0.65)
        elif (h == "m_t" or h == "m_tbar"):
            art.AddLegend(0.5, 0.7, 0.7, 0.9)
        else:
            art.AddLegend(0.76, 0.45, 0.92, 0.65)

    art.Save("~/Website/figures/dilepton/" + h + "-dd-AZX.pdf")
    i += 1
