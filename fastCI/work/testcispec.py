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
        pT = cinlo[-1].pT()
        print hname, histdir, pT[0], pT[-1]

    Lambda = 20.0
    l = 1.0/Lambda**2
    kappa = vector('double')(6, 0)

    setStyle()
    c = TCanvas('figs/CT14/fig_CI_Lambda_%2.2d' % int(Lambda), '',
                10, 10, 1000, 500)
    c.Divide(2, 1)

    hdummy = mkhist1('hdummy',
                         'Jet p_{T} (GeV)',
                         'd^{2}#sigma/dp_{T}dy',
                         50, 0, 4000)
    hdummy.SetMinimum(-3.e-3)
    hdummy.SetMaximum( 3.e-3)
    
    color = [kRed, kOrange, kYellow+2, kGreen+1, kBlue, kMagenta, kCyan]
    h = []
    for i, k in enumerate([-1, 1]):
        kolor = 0
        c.cd(i+1)
        gPad.SetGrid()
        hdummy.Draw()
        kappa[0] = k
        scribe = Scribe(0.50,0.80)
        option = 'c same'
        for j, CI in enumerate(cinlo):
            h.append( CI(l, kappa) )
            h[-1].SetLineWidth(2)
            h[-1].SetLineColor(color[kolor])
            #h[-1].SetLineStyle(j+1)
            h[-1].Draw(option)
            kolor += 1
            if kolor >= len(color):
                kolor = 0
        scribe.write('#sqrt{s} = 13 TeV')
        scribe.write('#Lambda  = %d TeV' % int(Lambda))
        scribe.write('#kappa   = (%d,0,0,0,0,0)' % k)
        c.Update()
        gSystem.ProcessEvents()
    c.SaveAs('.png')
    sleep(5)
    
#------------------------------------------------------------------------------ 
main()

