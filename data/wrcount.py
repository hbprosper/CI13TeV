#!/usr/bin/env python
import os, sys
import ROOT
from string import *

LUMI = 35.1 # 1/fb

def main():
    hfile = ROOT.TFile("data_13TeV_L035.1ifb_plus_corrections.root")
    hdata = hfile.Get("ak7/y_0.0-0.5/Njet").Clone("hdata")

    if hdata:
        print "total number of jets: %d" % hdata.Integral()
    else:
        sys.exit("can't get histogram")
        

    
    out = open('data_13TeV_L035.1ifb.txt', 'w')
    record = ' %9s %9s %9s' % ('lower', 'upper', 'count')
    print record
    out.write('%s\n' % record)

    total = 0
    M = 0
    for ii in xrange(hdata.GetNbinsX()):
        ibin  = ii + 1
        xlow  = hdata.GetBinLowEdge(ibin)
        if xlow < 700: continue

        M += 1
        width = hdata.GetBinWidth(ibin)
        xupp  = xlow + width
        count = hdata.GetBinContent(ibin)
        record = ' %9.0f %9.0f %9.0f' % (xlow, xupp, count)
        print record
        out.write('%s\n' % record)
        total += count
        
    print("number of bins: %10d" % M)
    print("total countL    %10d" % total)
main()

