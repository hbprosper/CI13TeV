#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Make a histogram of the jet count per bin and simultaneously
# create a table data.tex containing the observed counts
# Created 18-Dec-2015 @ CERN HBP
# Updated 25-Oct-2017 @ CERN HBP
#-------------------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from array import array
from time import sleep
from ROOT import *
#-------------------------------------------------------------------------------
NAME = 'data_13TeV_L035.1ifb'
INPFILENAME = '../data/%s.txt'  % NAME
OUTFILENAME = '../data/%s.root' % NAME
#-------------------------------------------------------------------------------
def divideByBinWidth(h):
    for ii in xrange(h.GetNbinsX()):
        count = h.GetBinContent(ii+1)
        error = h.GetBinError(ii+1)
        width = h.GetBinWidth(ii+1)
        h.SetBinContent(ii+1, count/width)
        h.SetBinError(ii+1, error/width)
        
def main():
    inp = open(INPFILENAME)
    table = '''
\\begin{table}[htp]
  \\caption{Jet yield for each $p_\\textrm{T}$ bin and $|y| < 0.5$.}
  \\label{tab:yield}
  \\medskip
  \\centering
  \\begin{tabular}{|r|rr|r|}
  \\hline
  bin\t& $p_\\textrm{T,min}$\t& $p_\\textrm{T,max}$\t& jet yield \\\\ \\hline
 '''
    pt = array('d')
    y  = []
    ptmin  = 0.0
    ptmax  = 0.0
    header = split(inp.readline())
    ii = 0
    while 1:
        t = split(inp.readline())
        if len(t) == 0: break
        ptmin, ptmax, count = map(atof, t)

        pt.append(ptmin)
        y.append(count)
        
        ii += 1
        record = '%4d\t&%8.0f\t&%8.0f\t&%10.0f\t\\\\ \\hline' % \
          (ii, ptmin, ptmax, count)
        print record
        table += '%s\n' % record
    pt.append(ptmax)
    
    table += '\\end{tabular}\n\\end{table}\n'
    open('data.tex','w').write(table)

    print

    # ---------------------------------------------------------------
    # Make and plot histogram
    # ---------------------------------------------------------------
    L = 35.1 # 1/fb

    setStyle() 
    hfile = TFile(OUTFILENAME, 'recreate')

    hdata = TH1D('hdata', '', len(y), pt)
    hdata.SetMinimum(2)
    hdata.SetMaximum(2e+6)    
    hdata.SetAxisRange(700, 4000)
    hdata.GetYaxis().SetTitle('count / bin')
    hdata.GetYaxis().SetTitleOffset(1.6)
    hdata.SetNdivisions(505, "Y")
                
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetXaxis().SetTitleOffset(1.2)
    hdata.SetNdivisions(505, "X")
    hdata.SetMarkerSize(0.8)
    
    for ii in xrange(len(y)):
        hdata.SetBinContent(ii+1, y[ii])
        hdata.SetBinError(ii+1, sqrt(y[ii]))

    print 'plot...'
    cdata = TCanvas('figures/fig_%s' % NAME, "observed", 10, 10, 500, 500)
    cdata.cd()
    gPad.SetLogy()
    hdata.Draw('ep')

    scribe = addTitle('CMS Preliminary  '\
                      '#surds=13TeV L=%5.2fpb^{-1}' % L ,
                      0.035)
    cdata.Update()
    cdata.SaveAs('.png')


    # now divide by bin width
    hfile.cd()
    hdata2 = hdata.Clone('hdata2')
    divideByBinWidth(hdata2)
    hdata2.Scale(1.0/L)
    hdata2.SetMinimum(1e-3)
    hdata2.SetMaximum(1e+3)
    hdata2.GetYaxis().SetTitle('d^{2}#sigma/dp_{T}dy (fb/GeV)')


    cdata2 = TCanvas('figures/fig_%s2' % NAME, "observed", 520, 10, 500, 500)
    cdata2.cd()
    gPad.SetLogy()
    hdata2.Draw('ep')

    scribe2 = addTitle('CMS Preliminary  '\
                       '#surds=13TeV L=%5.2ffb^{-1}' % L,
                      0.035)
    cdata2.Update()
    gSystem.ProcessEvents()
    cdata2.SaveAs('.png')
    
    sleep(5)
    hfile.Write()
    hfile.Close()
#------------------------------------------------------------------------------
try:    
    main()
except KeyboardInterrupt:
    print "ciao!"
    
