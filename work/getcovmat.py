#!/usr/bin/env python
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#------------------------------------------------------------------------------
def main():
    hfile = TFile('../data/CovMatrix-CDEFGH-ak7.root')
    covmat= hfile.Get('CovMatrix_1bin')
    
    setStyle() 

    covmat.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    covmat.GetXaxis().SetRangeUser(1032, 3832)
    covmat.GetXaxis().SetTitleOffset(1.4)
    covmat.SetNdivisions(505, "X")

    covmat.GetYaxis().SetTitle('Jet p_{T} (GeV)')
    covmat.GetYaxis().SetRangeUser(1032, 3832)
    covmat.GetYaxis().SetTitleOffset(1.4)
    covmat.SetNdivisions(505, "Y")

    hcov = covmat.Clone("covmat")
    zmin = hcov.GetMinimum()
    zmax = hcov.GetMaximum()
    print("zmin: %10.2e\tzmax: %10.2e" % (zmin, zmax))

    c = TCanvas('fig_covmat', "covmat", 10, 10, 800, 800)
    c.cd()
    c.SetLogz()
    hcov.Draw("lego2")
    
    c.Update()
    gSystem.ProcessEvents()
    c.SaveAs('.png')
    sleep(5)
#------------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
