#!/usr/bin/env python
import os, sys
import ROOT
from string import *

LUMI = 35.1 # 1/fb

def main():
    hfile = ROOT.TFile("forCI.root")
    hdata = hfile.Get("ak7/y_0.0-0.5/hptRun2016").Clone("hdata")

    if hdata:
        print "total number of jets: %d" % hdata.Integral()
    else:
        sys.exit("can't get histogram")
        

    
    out = open('xsec_13TeV_L035.1ifb.txt', 'w')
    record = ' %9s %9s %9s %9s' % ('lower', 'upper', 'count', 'x-sect.')
    print record
    out.write('%s\n' % record)

    M = 0
    for ii in xrange(hdata.GetNbinsX()):
        ibin  = ii + 1
        xlow  = hdata.GetBinLowEdge(ibin)
        if xlow < 700: continue

        M += 1
        width = hdata.GetBinWidth(ibin)
        xupp  = xlow + width
        xsec  = hdata.GetBinContent(ibin)
        count = xsec * width * LUMI
        record = ' %9.0f %9.0f %9.0f %9.3e' % (xlow, xupp, count, xsec)
        print record
        out.write('%s\n' % record)

    print("number of bins: %d" % M)        
main()

