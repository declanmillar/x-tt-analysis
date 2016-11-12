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
        # hist.SetLineStyle(30)
        # hist.SetLineWidth(1)
        # hist.GetXaxis().CenterTitle()
        # hist.GetYaxis().CenterTitle()
        # hist.GetYaxis().SetNdivisions(4)
        # hist.GetXaxis().SetTitleSize(0.06)
        # hist.GetYaxis().SetTitleSize(0.06)
        # hist.GetZaxis().SetTitleSize(0.06)
        # hist.GetXaxis().SetTitleOffset(1.2)
        # hist.GetYaxis().SetTitleOffset(1.2)
        # hist.GetZaxis().SetTitleOffset(0.8)
        # hist.SetMinimum(hist.GetBinContent(hist.GetMinimumBin()))
        # hist.Draw("lego2 same fb")

        hist.Draw("hist same")

        # clone to draw errors
        clone = hist.Clone("")
        self.members.append(clone)

        # Uncertainty style
        clone.SetFillColorAlpha(color, 0.2)
        clone.SetFillStyle(1001)
        clone.DrawCopy("e2 same")

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
        top_margin = 0.05
        right_margin = 0.05
        left_margin = 0.12
        bottom_margin = 0.15

        # 2D
        # right_margin = 0.03
        # bottom_margin = 0.08

        upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
        upper_pad.Draw()
        upper_pad.cd()

        # if True:#"TH2D" in type(hist).__name__:
            # upper_pad.SetTheta(20) # default: 30
            # # upper_pad.SetPhi(-30) # default: 30
            # upper_pad.SetPhi(-150)
            # upper_pad.Update()

        if (self.setlogy):
          upper_pad.SetLogy()
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


    def Save(self, name):
        HistPainter.canvas.SaveAs(name)


    def SetStyle(self):
        import atlas_style


red = ROOT.TColor.GetColor(250.0 / 255.0, 0.0 / 255.0, 0.0 / 255.0)
blue = ROOT.TColor.GetColor(0.0 / 255.0, 90.0 / 255.0, 130.0 / 255.0)
green = ROOT.TColor.GetColor(0.0 / 255.0, 130.0 / 255.0, 90.0 / 255.0)
grey = ROOT.TColor.GetColor(64.0 / 255.0, 64.0 / 255.0, 64.0 / 255.0)
black = ROOT.TColor.GetColor(1.0, 1.0, 1.0)

hs = ["m_t", "mt_bar", "m_tt", "pT_t", "pT_tbar", "E_t", "Et_bar"]#, "deltaR_tt"]
hs2 = ["mt", "mtbar", "mtt", "pTt", "pTtbar", "Et", "Etbar"]#, "deltaRtt"]

i = 0
for h in hs:
    art = HistPainter(1920, 1080)
    art.SetStyle()
    art.AddPads()
    # art.SetXtitle("#it{m_{tt}} [TeV]")
    # art.SetYtitle("d#it{#sigma} / d#it{m_{tt}} [pb/TeV]")
    # art.SetYtitle("A^{*}_{FB}")
    # art.SetZtitle("#it{cos#theta_{l}}")

    art.AddHistogram(h, "GLR-R-3_qq-AZX-tt-bbllvv_0-4_5x10M.a.root", "#bf{SM}", blue)
    art.AddHistogram(hs2[i] + "_R", "GLR-R-3_qq-AZX-tt-bbllvv_0-4_5x10M.a.root", "#bf{SM}", red, 2)

    art.AddInfoBox("GLR-R", 0.75, 0.85)

    art.SetHistTitle(0, "generated")
    art.SetHistTitle(1, "reconstructed")
    # art.AddLegend(0.13, 0.7, 0.35, 0.9)
    art.AddLegend(0.76, 0.45, 0.92, 0.65)

    art.Save("~/Website/dilepton-plots/" + h +"_GLR-R-3_uudd_0-4.pdf")
    i += 1
