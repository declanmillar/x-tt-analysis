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

        # hist.SetMinimum(0.3)
        # hist.Scale(1 / hist.Integral())


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
            hist.GetXaxis().SetTitleOffset(1.4)
            hist.GetYaxis().SetTitleOffset(1.3)
            hist.GetZaxis().SetTitleOffset(0.7)
            hist.SetMinimum(hist.GetBinContent(hist.GetMinimumBin()))
            hist.Draw("lego same fb")
            # hist.Draw("lego2 same fb")
            # hist.Draw("colz fb")
        else:
            hist.Draw("hist same")

            # clone to draw errors
            clone = hist.Clone("")
            self.members.append(clone)

            # Uncertainty style
            clone.SetFillColorAlpha(color, 0.3)
            # clone.SetFillStyle(1001)
            clone.SetFillStyle(3554)
            # if not ".SM." in filename: clone.DrawCopy("e2 same")


        f.Close()


    def SetHistType(self, histname, filename):

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


    def AddInfoBox(self, model, xmass, xratio, energy, luminosity, efficiency, x, y):

        # if model == "SM":
        #     t = ROOT.TLatex(0.85, y, "#bf{Model: %s} " % model)
        #     t.SetNDC(True)
        #     t.SetTextSize(0.04)
        #     t.Draw("same")
        #     self.members.append(t)
        # else:
        #     t = ROOT.TLatex(x, y, "#bf{Model: %s, #it{m_{Z'}} = %i TeV, #it{#Gamma_{Z'} / m_{Z'}} = %.1f %%}" % (model, xmass, xratio))
        #     t.SetNDC(True)
        #     t.SetTextSize(0.04)
        #     t.Draw("same")
        #     self.members.append(t)

        #bf{#it{m_{Z'}}  =  %i TeV, #it{#Gamma_{Z'} / m_{Z'}}  =  %.1f %%}

        # t = ROOT.TLatex(x, y - 0.075, "#bf{#it{m_{Z'}}  =  %i TeV, #it{#Gamma_{Z'} / m_{Z'}}  =  %.1f %%}" % (xmass, xratio)) #bf{#it{#Gamma_{Z'} / m_{Z'}}  =  %.1f %%}" % xratio
        # t.SetTextAlign(11)
        # t.SetTextSize(0.04)
        # t.SetNDC(True)
        # t.Draw("same")
        # self.members.append(t)

        t = ROOT.TLatex(x - 0.057, y - 0.21, "#bf{#it{#Gamma_{Z'} / m_{Z'}}  =  %.1f %%}" % xratio)
        t.SetTextAlign(11)
        t.SetTextSize(0.04)
        t.SetNDC(True)
        t.Draw("same")
        self.members.append(t)

        # t = ROOT.TLatex(x - 0.5, y, "#bf{#it{#sqrt{s}}  =  %i TeV, #int #it{L dt}  =  %i fb^{-1}}" % (energy, luminosity))
        # t.SetNDC(True)
        # t.SetTextSize(0.04)
        # t.Draw("same")
        # self.members.append(t)

        t = ROOT.TLatex(x - 0.0395, y, "#bf{#int #it{L dt}  =  %i fb^{-1}}" % luminosity)
        t.SetNDC(True)
        t.SetTextSize(0.04)
        t.Draw("same")
        self.members.append(t)

        t = ROOT.TLatex(x - 0.002, y - 0.07, "#bf{#it{#sqrt{s}}  =  %i TeV}" % (energy))
        t.SetNDC(True)
        t.SetTextSize(0.04)
        t.Draw("same")
        self.members.append(t)

        t = ROOT.TLatex(x - 0.01, y - 0.14, "#bf{#it{m_{Z'}}  =  %i TeV}" % 4)
        t.SetNDC(True)
        t.SetTextSize(0.04)
        t.Draw("same")
        self.members.append(t)

        # t = ROOT.TLatex(x, y - 0.375, "#bf{\it{#epsilon}  =  %.1f}" % efficiency)
        # t.SetNDC(True)
        # t.SetTextSize(0.04)
        # t.Draw("same")
        # self.members.append(t)


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
            left_margin = 0.1
            right_margin = 0.03
            bottom_margin = 0.15
        else:
            top_margin = 0.05
            right_margin = 0.05
            left_margin = 0.12
            bottom_margin = 0.15

        upper_pad.SetMargin(left_margin, right_margin, bottom_margin, top_margin)
        upper_pad.Draw()
        upper_pad.cd()

        if self.is_th2d:
            upper_pad.SetTheta(20) # default: 30
            upper_pad.SetPhi(-30) # default: 30
            # upper_pad.SetPhi(-150)
            upper_pad.Update()

        # upper_pad.SetGridx()
        # upper_pad.SetGridy()


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


#############################################################################


red = ROOT.TColor.GetColor(250.0 / 255.0, 0.0 / 255.0, 0.0 / 255.0)
orange = ROOT.TColor.GetColor(255.0 / 255.0, 165.0 / 255.0, 0.0 / 255.0)
blue = ROOT.TColor.GetColor(0.0 / 255.0, 90.0 / 255.0, 130.0 / 255.0)
green = ROOT.TColor.GetColor(0.0 / 255.0, 130.0 / 255.0, 90.0 / 255.0)
grey = ROOT.TColor.GetColor(64.0 / 255.0, 64.0 / 255.0, 64.0 / 255.0)
black = ROOT.TColor.GetColor(0.0, 0.0, 0.0)

# hs = []
hs = ["AL1", "AL_R"]

models = ["GLR-R-4"]
# models = ["GLR-R-4", "GLR-Y-4", "GLR-LR-4"]
# models = ["GSM-Q-4", "GLR-BL-4"]
# models = ["GSM-T3L-4", "GSM-SM-4"]
widths = [2.5, 2.4, 2.1, 4.7]

fs = []

for model in models:
    fs.append("ggqqdduu-AZX-tt-bbllvv." + model + ".13TeV.CT14LL.2-5.weighted.r2.L300.root")
    fs.append("ggqqdduu-AZX-tt-bbllvv." + model + ".13TeV.CT14LL.2-5.weighted.r2.L300.fid.root")

# models.append("SM")
# fs.append("ggqqdduu-AZ-tt-bbllvv.SM.13TeV.CT14LL.3-5.20x2M.weighted.r2.L300.fid.root")

# j = 0
# for f in fs:
i = 0
# for h in hs:
art = HistPainter(1920, 1080)

art.SetHistType(hs[0], fs[0])

# art.SetLogy()
# art.SetLogz()
art.SetStyle()
art.AddPads()

# art.SetDomain(3.025, 4.975)
if "AtFB" in hs[0]: art.SetRange(-0.5, 0.5)
if "AL" in hs[0]: art.SetRange(-1.0, 1.0)

art.SetXtitle("#it{m_{tt}} [TeV]")
if "AtFB" in hs[0]:
    art.SetYtitle("#it{A^{t}_{FB^{*}}}")
elif "AL" in hs[0]:
    art.SetYtitle("#it{A_{L}}")
elif "mtt_costhetastar_R" in hs[0]:
    art.SetYtitle("#it{cos#theta*}")
    art.SetZtitle("Expected events")
elif "mtt_costhetal_R" in hs[0]:
    art.SetYtitle("#it{cos#theta_{l}}")
    art.SetZtitle("Expected events")
else: 
    # art.SetYtitle("d#it{#sigma} / d#it{m_{tt}} [pb/TeV]")
    art.SetYtitle("Expected events")

art.AddHistogram(hs[0], fs[0], " #bf{truth}", black, 1)
art.AddHistogram(hs[1], fs[0], " #bf{recon}", black, 2)
art.AddHistogram(hs[0], fs[1], " #bf{truth, p_{T} > 25 GeV, |#eta| < 2.5}", red, 1)
art.AddHistogram(hs[1], fs[1], " #bf{recon, p_{T} > 25 GeV, |#eta| < 2.5}", red, 2)

art.AddInfoBox(models[0].rsplit('-',1)[0], 4, 2.5, 13, 300, 0.6, 0.77, 0.87)
art.AddLegend(0.15, 0.52, 0.45, 0.92)
# art.AddLegend(0.15, 0.16, 0.45, 0.52)

art.Save("~/Desktop/" + hs[0] + "-" + "truth-vs-reco-vs-fid.pdf")
# i += 1
# j += 1
