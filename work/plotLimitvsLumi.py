#!/usr/bin/env python
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
DATA = '''
100			16.7	16.2	11.2	11.1
200			19.3	18.4	12.1	12.1
400			22.3	20.8	13.0	13.1
1000		27.2	24.6	14.7	14.6
1600		29.8	26.4	15.6	15.1
3200		34.3	29.4	19.7	15.9
10000		43.2	34.6	37.1	17.2
30000		53.6	38.8	49.4	18.5
'''
DATA = map(lambda x: map(atof, x),
           map(split, split(strip(DATA), '\n')))
#-----------------------------------------------------------------------------
def main():
    setStyle()
    
    x  = array('d')
    y1 = array('d')
    y2 = array('d')
    y3 = array('d')
    y4 = array('d')
    for lumi, l1, l2, l3, l4 in DATA:
        lumi /= 1000
        x.append(lumi)
        y1.append(l1)
        y2.append(l2)
        y3.append(l3)
        y4.append(l4)

    g1 = mkgraph(x, y1, 'integrated luminosity (1/fb)', '#Lambda_{95}',
                 0, 20, ymin=0.0, ymax=60,
                 color=kBlue)

    g2 = mkgraph(x, y2, '#int L(t)dt (1/fb^)', '#Lambda_{95}',
                 0, 20, ymin=0.0, ymax=60,
                 color=kRed)

    legend = mklegend(0.24, 0.75, 0.20, 0.20)
    legend.AddEntry(g1, 'No syst.', 'l')    
    legend.AddEntry(g2, 'syst.', 'l')

    
    c1 = TCanvas('figures/fig_limitsvslumi_c', "constructive", 10, 10, 500, 500)
    c1.cd()
    g1.Draw('ac')
    g2.Draw('c same')
    legend.Draw()
    scribe = addTitle('CMS Preliminary  #surds=13TeV, CT14, LL^{-}', 0.035)
    c1.Update()
    c1.SaveAs('.png')


    g3 = mkgraph(x, y3, 'integrated luminosity (1/fb)', '#Lambda_{95}',
                 0, 20, ymin=0.0, ymax=60,
                 color=kBlue)

    g4 = mkgraph(x, y4, '#int L(t)dt (1/fb^)', '#Lambda_{95}',
                 0, 20, ymin=0.0, ymax=60,
                 color=kRed)

    legend = mklegend(0.24, 0.75, 0.20, 0.20)
    legend.AddEntry(g1, 'No syst.', 'l')    
    legend.AddEntry(g2, 'syst.', 'l')

    
    c2 = TCanvas('figures/fig_limitsvslumi_d', "destructive", 10, 520,
                 500, 500)
    c2.cd()
    g3.Draw('ac')
    g4.Draw('c same')
    legend.Draw()
    scribe = addTitle('CMS Preliminary  #surds=13TeV, CT14, LL^{+}', 0.035)
    c2.Update()
    c2.SaveAs('.png')
    

    sleep(10)
#-----------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
