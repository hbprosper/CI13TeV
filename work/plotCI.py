#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotWorkspace.py
# read spectra and plot them
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re, optparse
from math import *
from histutil import *
from string import *
from array import array
from time import sleep
from ROOT import *
#-----------------------------------------------------------------------------
gSystem.Load("../CI/lib/libCI.so")
#from ROOT import hutil
    
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-8  # minimum and maximum y-limits for spectrum
YMAX    = 2.0
#                        kappa
#                   1   2   3   4   5   6
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']

LUMI    = 35100.0  # 1/pb

class Context:
    pass    
#-----------------------------------------------------------------------------
getpdfdir   = re.compile('(?<=fastNLO/).+(?=/[0-9])|(?<=fastCI/).+(?=/[0-9])')
getsmearing = re.compile('(?<=/[0-9][0-9][0-9]/).+')
PDFS = {'ALL'  : 'CT14nlo, MMHT2014nlo, NNPDF30_nlo',
        'CT14' : 'CT14nlo',
        'MMHT' : 'MMHT2014nlo',
        'NNPDF': 'NNPDF30_nlo'}

def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
def waitForKey(message):
    print message
    sys.stdout.flush()
    raw_input("")
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '24-Jun-2016'
    USAGE = '''
    python plotCI.py [options] <workspace-file>

    options
       -m<model>      LL, RR, VV, AA, V-A     [ALL]
       -L<Lambda>                             [20 (TeV)]
       -r<reference>  median, nominal (QCD)   [nominal]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)

    parser.add_option('-L', '--Lambda',
                      action='store',
                      dest='Lambda',
                      type='string',
                      default='20',
                      help='CI mass scale')

    parser.add_option('-m', '--model',
                      action='store',
                      dest='model',
                      type='string',
                      default='LL RR VV AA V-A',
                      help='model')
        
    options, args = parser.parse_args()
    if len(args) == 0:
        sys.exit(USAGE)

    filename = args[0]

    # Set up model
    Lambda = float(options.Lambda)
    kappa = vector('double')(6, 0)
    kappa[0] =-1
    model = upper(options.model)
    if model[0] == '-':
        model = model[1:]
        key   = model
        sign  =-1
        interf= 'd'
    else:
        key   = model
        sign  = 1
        interf= 'c'
    if KAPPA.has_key(key):
        for ii in xrange(kappa.size()):
            kappa[ii] = sign*KAPPA[key][ii]

    return (filename, Lambda, kappa, model, interf)    
#-----------------------------------------------------------------------------
def makePlot(context):
    model  = context.model
    Lambda = context.Lambda
    kappa  = context.kappa
    interf = context.interf
    hqcd   = context.hqcd
    pT     = context.pT
    cispectrum  = context.cispectrum
    
    print "\n\t<=== %s ===> %s" % (context.model, context.interf)
    
    lam   = 1.0/Lambda**2
    pTmin = pT[0]
    pTmax = pT[-1]

    # --------------------------------------------------------
    # number of bins and range
    # --------------------------------------------------------    
    nbins  = hqcd.GetNbinsX()
    pTLOW  = hqcd.GetBinLowEdge(1)
    pTHIGH = hqcd.GetBinLowEdge(nbins)+hqcd.GetBinWidth(nbins)
        
    # --------------------------------------------------------
    # determine nature of smearing 
    # --------------------------------------------------------    
    dirname = cispectrum[0].dirname()
    smearing= getsmearing.findall(dirname)[0]
    
    if smearing == 'JEC':
        label  = 'JEC uncertainties only'
    elif smearing == 'PDF':
        label  = 'PDF uncertainties only'
    else:
        label  = 'JEC + PDF uncertainties'
    postfix = '_%s' % smearing

    # get PDF (directory) name
    prefix = nameonly(dirname)
    pdfdir = split(dirname, '/')[2]
    
    # --------------------------------------------------------
    print "="*80
    postfix += '_%s_%3.1f_%s' % (model, Lambda, interf[0])
    pcCI = PercentileCurve(nbins)
    hCI  = []    
    
    for index, CI in enumerate(cispectrum[1:]):
        hci = CI(lam, kappa)
        hutil.divideByWidth(hci) # convert to a density
        pcCI.add(hci)
        hCI.append(hci)

        if index % 1000 == 0:
            print "%5d" % index, hci.GetName()
    
    # --------------------------------------------------------
    # plot spectrum
    # --------------------------------------------------------
    os.system('mkdir -p figures/%s' % pdfdir)
    
    name = 'figures/%s/%s_xsection%s' % (pdfdir, prefix, postfix)
    print
    print name
    cspect = TCanvas(name, name, 10, 10, 500, 500)
    x = map(lambda i: (pT[i+1]+pT[i])/2, range(nbins))

    # decide what spectrum is to be plotted
    pc   = pcCI
    hist = hCI 
        
    n = hqcd.GetNbinsX()
    print "QCD0: bins %d, pT-range = (%-6.1f ... %-6.1f) GeV" % \
      (n, hqcd.GetBinLowEdge(1),
       hqcd.GetBinLowEdge(n)+hqcd.GetBinWidth(n))

    ymin =-2.e-5
    ymax = 2.e-5
    hist[0].SetMinimum(ymin)
    hist[0].SetMaximum(ymax)
        
    # compute percentile spectra
    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95 = mkpline(x, curve[0], curve[-1], hist[0], cspect, color=kGreen)
    p68 = mkpline(x, curve[1], curve[-2], hist[0], cspect, color=kYellow)
    p50 = mkgraph(x, curve[2],
                  "Jet p_{T} (GeV)",
                  "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                  pTmin, pTmax,
                  ymin=ymin, ymax=ymax,
                  color=kRed,
                  lwidth=1)
    p50.SetLineWidth(2)
    
    cspect.cd()
    cspect.SetGrid()
    p50.Draw('ac')
    p95.Draw('f same')
    p68.Draw('f same')
    p50.Draw('c same')
    hqcd.Draw('c same')

    bigtab = ' '*6
    energy = context.energy
    lumi   = context.lumi
    scribe = addTitle('%sCMS Preliminary  #surds=%sTeV L=%s/fb' % \
                      (' '*10, energy, lumi),
                      0.035)
    scribe.vspace()
    scribe.write("%s(#Lambda = %3.1f TeV) %s" % (model,
                                                     Lambda,
                                                     interf), 0.05)
    scribe.write(PDFS[pdfdir], 0.07)
    scribe.write(label, 0.11)

    xp = 0.70
    yp = 0.68
    xw = 0.18
    yw = 0.24
    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hqcd, 'QCD', 'l')
    lg.AddEntry(p50,  'median', 'l')
    lg.AddEntry(p68,  '68%s' % '%', 'f')
    lg.AddEntry(p95,  '95%s' % '%', 'f')
    lg.Draw('same')
    
    cspect.Update()
    gSystem.ProcessEvents()
    cspect.SaveAs('.png')
    
    # --------------------------------------------------------
    # now plot the ratio
    # --------------------------------------------------------
    cratio = TCanvas('figures/%s/%s_xsection%s_ratio' % (pdfdir,
                                                         prefix, postfix),
                     '%s - xsection-ratio' % dirname,
                     515, 10, 500, 500)


    title = 'CI / QCD_{NLO}#otimesNP#otimesEWK'
    color = [kRed, kOrange, kGreen, kBlue, kMagenta]
    pc  = PercentileCurve(nbins)
    jj = 0
    for ii, h in enumerate(hist):
        h.Divide(hqcd)
        pc.add(h)

    curve = []
    for p in PERCENT: curve.append( pc(p) )

    ymin = 0.0
    ymax = 3.0
    hist[0].SetMinimum(ymin)
    hist[0].SetMaximum(ymax)
    p95r = mkpline(x, curve[0], curve[-1], hist[0], cratio, color=kGreen)
    p68r = mkpline(x, curve[1], curve[-2], hist[0], cratio, color=kYellow)
    p50r = mkgraph(x, curve[2],
                   "Jet p_{T} (GeV)", title,
                   pTmin, pTmax,
                   ymin=ymin, ymax=ymax,
                   color=kRed,
                   lwidth=1)
        
    cratio.cd()
    cratio.SetGrid()    
    p50r.Draw('ac')
    p95r.Draw('f same')
    p68r.Draw('f same')
    p50r.Draw('c same')
    
    scriber = addTitle('%sCMS Preliminary  #surds=%sTeV L=%s/fb' % \
                       ('', energy, lumi),
                       0.035)
    scriber.vspace()
    scriber.write("%s(#Lambda=%3.1f TeV) %s" % (model,
                                                    Lambda,
                                                    interf), 0.04)
    scriber.write(PDFS[pdfdir], 0.04)
    scriber.write(label, 0.04)

    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(p50r, 'median', 'l')
    lg.AddEntry(p68r, '68%s' % '%', 'f')
    lg.AddEntry(p95r, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cratio.Update()
    gSystem.ProcessEvents()
    cratio.SaveAs('.png')
    sleep(5)    
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotWorkspace.py ===>"
    setStyle()

    filename, Lambda, kappa, model, interf = decodeCommandLine()

    context = Context()
    context.filename = filename
    context.Lambda = Lambda
    context.kappa  = kappa
    context.energy = '13'
    lumi = LUMI/1000
    context.lumi   = '%5.2f' % lumi
        
    # --------------------------------------------------------        
    print "\nloading workspace..."
    # --------------------------------------------------------    
    hfile = TFile(filename)
    ws = hfile.Get('CI')
    check(ws, "can't get workspace CI")

    # --------------------------------------------------------      
    # get model and CI spectra
    # --------------------------------------------------------    
    pdf = ws.pdf('model')
    nspectra = pdf.size()
    print "number of spectra %d" % nspectra

    cispectrum = []
    for ii in xrange(nspectra):
        ci = pdf.CI(ii)
        if ci == None: break
        cispectrum.append(ci)

        if ii % 1000 == 0:
            print "%4d %s %s" % (ii, ci.dirname(), ci.histname())

    context.cispectrum  = cispectrum

    # --------------------------------------------------------      
    # get nominal QCD spectrum
    # --------------------------------------------------------
    hqcd = pdf.QCD(0)()
    check(hqcd, "can't get hqcd")
    nbins   = hqcd.GetNbinsX()
    hqcd.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hqcd.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hqcd.SetMarkerSize(0.8)
    hqcd.SetMarkerStyle(20)
    hutil.divideByWidth(hqcd) # convert to density
    ymin = hqcd.SetMinimum(YMIN)
    ymax = hqcd.SetMaximum(YMAX)
    pT = hutil.binlowedges(hqcd)
    pT.push_back(pT[-1] + hqcd.GetBinWidth(nbins))

    context.hqcd = hqcd
    context.pT = pT
    
    # --------------------------------------------------------
    # now plot
    # --------------------------------------------------------
    models = split(model)
    print "\nplotting..."
    for model in models:
        context.model = model
        
        context.interf= 'constructive'
        print '\=> interference: %s' % context.interf
        for ii in xrange(kappa.size()):
            kappa[ii] = KAPPA[model][ii]
        context.kappa = kappa
        makePlot(context)

        context.interf= 'destructive'                        
        print '\=> interference: %s' % context.interf        
        for ii in xrange(kappa.size()):
            kappa[ii] =-KAPPA[model][ii]
        context.kappa = kappa
        makePlot(context)
            
    sleep(5)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
