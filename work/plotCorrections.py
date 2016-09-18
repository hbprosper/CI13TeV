#!/usr/bin/env python
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
INPFILENAME = '../data/ElectroWeakCorrection.root'
#-----------------------------------------------------------------------------
def main():

    gSystem.Load('libCI.so')

    def funNP(x, p):
        # from Paolo Gunnellini (2015-12-17)
        A = p[0]
        B = p[1]
        n = p[2]
        pT= x[0]
        return A + B / pT**n
    fNP = TF1('NPcor', funNP, 500, 2500, 3)
    fNP.SetParameters(1.00139, 433.922, 1.67615)
    fNP.SetLineColor(kRed)
    fNP.SetLineWidth(2)

    setStyle()        
    hfile = TFile(INPFILENAME)
    if not hfile.IsOpen():
        print "** can't open file %s" % INPFILENAME
        sys.exit(0)

    htmp1 = mkhist1('htmp1', 'Jet p_{T} (GeV)',
                    'Correction Factor', 50, 500, 2500)
    htmp1.SetMinimum(0.8)
    htmp1.SetMaximum(1.2)
    htmp1.SetNdivisions(505, "Y")
    
    h = hfile.Get("histo")
    h.GetYaxis().SetTitle('Corrections')
    h.GetYaxis().SetTitleOffset(1.6)
    h.SetNdivisions(505, "Y")
                
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.GetXaxis().SetTitleOffset(1.2)
    h.SetNdivisions(505, "X")
    h.SetMinimum(0)
    h.SetMaximum(1.5)
    h.SetLineColor(kBlue)
    h.SetLineWidth(2)
    h.SetAxisRange(500, 2500)

    legend = mklegend(0.24, 0.75, 0.62, 0.18)
    legend.AddEntry(h,   'Electroweak correction', 'l')    
    legend.AddEntry(fNP, 'Non-pertubative correction', 'l')

    
    c1 = TCanvas('figures/fig_corrections', "Corrections", 10, 10, 500, 500)
    c1.cd()
    gPad.SetGrid()
    htmp1.Draw()
    fNP.Draw('c same')
    h.Draw('c same')
    legend.Draw()
    scribe = addTitle('CMS Preliminary  #surds=13TeV', 0.035)
    c1.Update()
    c1.SaveAs('.png')

    sleep(10)
#-----------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
