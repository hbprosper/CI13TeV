#!/usr/bin/env python
#-------------------------------------------------------------------------------
import os, sys, re
from time import sleep
from histutil import *
from ROOT import *
#-------------------------------------------------------------------------------
HISTNAMES = '''
_0.500_0.500
_0.500_1.000
_1.000_0.500
_1.000_1.000
_1.000_2.000
_2.000_1.000
_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
#-------------------------------------------------------------------------------
def main():
    gSystem.Load('libCI.so')
    
    cinlo= []
    print '\t...loading NLO CI histograms...'
    member = 0
    histdir = '../CT14/%3.3d' % member
    for postfix in HISTNAMES:
        hname = 'nlo%s' % postfix
        cinlo.append( CISpectrum(hname, hname, histdir, hname) )
        print hname

    Lambda = 30.0
    l = 1.0/Lambda**2
    kappa = vector('double')(6, 0)

    setStyle()
    c = TCanvas('figs/CT14/fig_CI_Lambda_%2.2d' % int(Lambda), '',
                10, 10, 1000, 500)
    c.Divide(2, 1)

    color = [kRed, kOrange, kYellow+3, kGreen+2, kBlue, kMagenta, kCyan+1]
    h = []
    for i, (k, ymin, ymax) in enumerate([(-1, -0.001, 0.001),
                                         ( 1, -0.001, 0.001)]):
        c.cd(i+1)
        gPad.SetGrid()        
        kappa[0] = k
        scribe = Scribe(0.50,0.80)
        option = 'l'
        for j, CI in enumerate(cinlo):
            h.append( CI(l, kappa) )
            h[-1].SetMinimum(ymin)
            h[-1].SetMaximum(ymax)
            h[-1].SetLineColor(color[j])
            h[-1].SetLineStyle(j+1)
            h[-1].SetLineWidth(3)
	    h[-1].Draw(option)
            option = 'l same'
        scribe.write('#sqrt{s} = 13TeV')
        scribe.write('#Lambda  = %dTeV' % int(Lambda))
        scribe.write('#kappa  = (%d,0,0,0,0,0)' % k)
        c.Update()
    c.SaveAs('.png')
    gSystem.ProcessEvents()
    sleep(2)
    
#------------------------------------------------------------------------------ 
main()

