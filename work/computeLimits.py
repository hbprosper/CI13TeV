#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        computeLimits.py
# Description: read CI workspace and compute limits
# Created:     02-Mar-2014 HBP
#              22-May-2014 HBP slighty modified version of program
#                          used to compute 7 TeV limits
#              09-Jun-2014 HBP add check function
#              15-Nov-2014 HBP use new pdf
#              18-Sep-2016 HBP
#              18-May-2018 HBP permit model-dependent lambda range
#-----------------------------------------------------------------
import os,sys,re
from array import array
from time import sleep
from histutil import *
from string import *
from ROOT import *
#-----------------------------------------------------------------
EXPECTED = True    # if true compute expected limits

ENERGY   = 13      # TeV
LUMI     = 35100.0 # 1/fb
BINMIN   =  1      # first bin to use (ROOT bin number convention)
BINMAX   = 28      # last bin to use
LMIN     = 0.0     # lower limit of lambda = 1/Lambda^2
LMAX     = {'LL': {1: 0.0015, -1: 0.004},
            'RR': {1: 0.0015, -1: 0.004},
            'VV': {1: 0.0015, -1: 0.004},
            'AA': {1: 0.0015, -1: 0.004},
            'V-A':{1: 0.0015, -1: 0.004}} # upper limit on lambda

YMAX     = 0.060   # maximum Y limit of posterior density plot
WSPACE   = 'CI'    # name of workspace

RELTOL   = 1.e-4   # minimum ratio of min(like)/max(like)
# Models
#                        kappa
#                   0   1   2   3   4   5
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']
#-----------------------------------------------------------------
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
#-----------------------------------------------------------------        
def makePlot(ws, likelihood, q,
             hname,
             color=kBlack,
             lstyle=1,
             lwidth=2,
             ymin=0.0,
             ymax=YMAX):

    xmin   = q.getMin()
    xmax   = q.getMax()
    xbins  = q.getBins()
    xtitle = q.GetTitle()
    ytitle = likelihood.GetTitle()
    xstep  = (xmax-xmin)/xbins

    h = TH1D(hname, "", xbins, xmin, xmax)

    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.GetXaxis().SetTitle('#lambda (TeV^{-2})')
    h.GetXaxis().SetNdivisions(505)

    h.GetYaxis().SetTitle('p(#lambda | D)')
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetTitleOffset(1.7)
    h.SetLineWidth(lwidth)

    ymax = 0.0
    reltol = RELTOL
    scan = True
    for ii in xrange(xbins):
        x = xmin + (ii+0.5)*xstep
        q.setVal(x)
        y = likelihood.getVal()
        h.SetBinContent(ii+1, y)
        h.SetBinError(ii+1, 0)
        if y > ymax: ymax = y
        if scan:
            if ymax > 0:
                if abs(y / ymax) < reltol:
                    scan = False
                    xmax = x
                    
    total = h.Integral()
    if total == 0:
        print "total is ZERO!"
        sys.exit(0)
        
    h.Scale(1.0/total)
    ymax = h.GetMaximum()
    ii   = int(ymax/0.04)
    ymax = (ii+1)*0.04
    h.SetMinimum(ymin)
    h.SetMaximum(ymax)

    #print "\tx-range: (%-6.4f ... %6.4f)" % (xmin, xmax)
    #print "\ty-range: (%-6.4f ... %6.4f)" % (ymin, ymax)
    return (h, xmin, xmax, ymin, ymax)
#-----------------------------------------------------------------
def setValues(ws, d, scale=1.0):
    if type(d) == type(""): d = ws.set(d)
    iterator = d.iterator()
    ii = 0
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        ws.var(name).setVal(d[name].getVal()/scale)
        print "%5d\t%10.2f" % (ii, ws.var(name).getVal())
        ii += 1
#-----------------------------------------------------------------
def toVector(ws, setname):
    data = ws.set(setname)
    if data == None: return None
    iterator = data.iterator()
    vdata = vector('double')()
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        vdata.push_back( ws.var(name).getVal() )
    return vdata
#-----------------------------------------------------------------
def main():
    # --------------------------------------
    # load various codes needed for the
    # calculatioms
    # --------------------------------------
    gSystem.Load('libCI')

    # set up some standard graphics style
    setStyle()

    # --------------------------------------
    # load workspace into memory
    # --------------------------------------
    argv = sys.argv[1:]
    if len(argv) == 0:
        print '''
    Usage:
         ./computeLimits.py <workspace-file> [LL RR etc.]
        '''
        sys.exit(0)
        
    filename = argv[0]
    wfile = TFile(filename)
    if not wfile.IsOpen():
        print "** can't open %s" % filename
        sys.exit(0)
    ws = wfile.Get(WSPACE)
    check(ws, "can't access workspace %s" % WSPACE)

    if len(argv) > 1:
        models = argv[1:]
    else:
        models = MODEL

    energy ='%d' % ENERGY
    lumi   = '%5.2f' % LUMI
        
    prefix = nameonly(filename)
    dname  = split(prefix, '_')[0]
    os.system('mkdir -p figures/%s' % dname)
    
    # --------------------------------------
    # get model etc.
    # --------------------------------------
    model = ws.pdf('model')
    Nset  = ws.set('Nset')
    poi   = ws.var('lambda')
    data  = toVector(ws, 'Nset')
    nbins = model.numberOfBins()

    # if EXPECTED = True then use internally
    # created Asimov data set
    model.setAsimov(EXPECTED, LUMI)
    if EXPECTED:
        print "\n\t== Use Asimov data set (Lumi = %10.1f / pb)" % LUMI
        Asimov = model.Asimov()
        for ii in xrange(Asimov.size()):
            print "\t%4d\t%10.1f" % (ii+1, Asimov[ii])

    # --------------------------------------
    # set range of bins to use.
    # NOTE: use ROOT bin labeling convention
    # first bin is 1 and last is nbins
    # --------------------------------------
    binmin = BINMIN
    binmax = BINMAX
    model.setBinRange(binmin, binmax)
         
    print "="*80
    print "input filename:  %s" % filename
    print "workspace:       %s" % WSPACE
    print "models:          %s" % models
    print "bin range:       [%d ... %d]" % (binmin, binmax)
    print "integrated lumi: %8.2f / pb" % LUMI
    print "="*80

    # --------------------------------------
    # create a graph of likelihood
    # --------------------------------------

    # create wrapper for model
    pdf = PDFWrapper(model, Nset, poi)
    CL  = 0.95  # confidence level

    # loop over models for which limits are to be calculated
    lmin = LMIN
    for key in models:

        for sign in [1, -1]:

            lmax = LMAX[key][sign]
            
            # set range depending on whether we have constructive or destructive
            # interference
            poi.setRange(lmin, lmax)
            print "lambda-range:    (%-6.3f ... %-6.3f) 1/TeV^2" % (lmin, lmax)
                
            # delete histograms to avoid memory leaks
            try:
                del hnom
            except:
                pass
            try:
                del havg
            except:
                pass
            try:
                del hclone
            except:
                pass

            Limit1 = -1
            Limit2 = -1
            
            if sign > 0:
                name = '%s_constructive' % key
                fname = 'figures/%s/%s_limit_%s_%5.5d' % \
                  (dname, prefix, name, int(LUMI+0.5))
                clike1 = TCanvas(fname, fname, 10, 10, 500, 500)
                clike = clike1
            else:
                name = '%s_destructive' % key
                fname = 'figures/%s/%s_limit_%s_%5.5d' % \
                  (dname, prefix, name, int(LUMI+0.5))
                clike2 = TCanvas(fname, fname, 515, 10, 500, 500)
                clike = clike2

            kappa = map(lambda x: sign*x, KAPPA[key])
            print "\n\tmodel: %5s\t%s" % (key, kappa)
            for ii in xrange(len(kappa)):
                vname = 'kappa%d' % ii
                ws.var(vname).setVal(kappa[ii])
                
            # --------------------------------------
            # compute nominal limits
            # --------------------------------------
            model.initialize(0) # use nominal cross section
                      
            hnom, xmin, xmax, ymin, ymax = makePlot(ws, model, poi,
                                                    "hnom",
                                                    color=kBlue, lstyle=2)
            swatch = TStopwatch()
            swatch.Start()

            bayes  = Bayes(pdf, data, xmin, xmax)
            limit1 = bayes.quantile(CL)
            if limit1 > 0: Limit1= 1.0/sqrt(limit1)
            print "\tLambda > %8.1f TeV @ %3.1f%s CL (no-syst.)" % \
              (Limit1, 100*CL, '%')
            clike.cd()
            hnom.Draw('l')
            clike.Update()

            # --------------------------------------
            # compute limits with systematic
            # uncertainties
            # --------------------------------------
            model.initialize(-1) # integrate over systematic uncertainties
                        
            havg, xmin, xmax, ymin, ymax  = makePlot(ws, model, poi,
                                                     "havg",
                                                     color=kRed,
                                                     lstyle=1,
                                                     lwidth=2)
            clike.cd()
            havg.Draw('l same')            
            clike.Update()

            bayes  = Bayes(pdf, data, xmin, xmax)        
            limit2 = bayes.quantile(CL)
            if limit2 > 0: Limit2 = 1.0/sqrt(limit2)
            print "\tLambda > %8.1f TeV @ %3.1f%s CL (syst.)" % \
              (Limit2, 100*CL, '%')
            print "\t\t\t==> real time: %8.3f s " % swatch.RealTime()

            # shade 95% region
            hclone = havg.Clone()
            hclone.SetAxisRange(lmin, limit2)
            hclone.SetFillStyle(3001)
            hclone.SetFillColor(30)
            clike.cd()
            hclone.Draw('l same')
                
            # --------------------------------------
            # plot posterior density and limits
            # --------------------------------------
            clike.cd()
            scribe = addTitle('CMS Preliminary  '\
                              '#surds=%sTeV  L=%spb^{-1}' % \
                              (energy, lumi),
                              0.04)
            scribe.vspace()
            scribe.write("model: %s, #kappa = %s" % (key, kappa), 0.06)
            scribe.write("#color[4]{#Lambda} > %4.1fTeV @ 95%s CL" % \
                         (Limit1, '%'), 0.06)
            scribe.write("#color[2]{#Lambda} > %4.1fTeV @ 95%s CL" % \
                         (Limit2, '%'), 0.06)
            clike.Update()
            clike.SaveAs('.png')
            
        sleep(5)
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
