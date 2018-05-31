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
    
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-8  # minimum and maximum y-limits for spectrum
YMAX    = 4.0
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
    python plotWorkspace.py [options] <workspace-file>

    options
       -m<model>      LL, RR, VV, AA, V-A, QCD     [ALL]
       -L<Lambda>                                  [20 (TeV)]
       -r<reference>  median, nominal (QCD)        [median]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)

    parser.add_option('-L', '--Lambda',
                      action='store',
                      dest='Lambda',
                      type='string',
                      default='25',
                      help='CI mass scale')

    parser.add_option('-m', '--model',
                      action='store',
                      dest='model',
                      type='string',
                      default='QCD LL RR VV AA V-A',
                      help='model')

    parser.add_option('-s', '--scale',
                      action='store_true',
                      dest='scale',
                      default=False,
                      help='normalize median spectrum to data')    
        
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
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

    return (filename, Lambda, kappa, model, interf, options.scale)    
#-----------------------------------------------------------------------------
def makePlot(context):
    model    = context.model
    Lambda   = context.Lambda
    kappa    = context.kappa
    interf   = context.interf
    hdata    = context.hdata
    pT       = context.pT
    scaleToData   = context.scaleToData
    qcdspectrum = context.qcdspectrum
    cispectrum  = context.cispectrum
    
    print "\n\t<=== %s ===> %s" % (context.model, context.interf)
    
    lam   = 1.0/Lambda**2
    pTmin = pT[0]
    pTmax = pT[-1]

    # --------------------------------------------------------
    # number of bins and range
    # --------------------------------------------------------    
    nbins  = hdata.GetNbinsX()
    pTLOW  = hdata.GetBinLowEdge(1)
    pTHIGH = hdata.GetBinLowEdge(nbins)+hdata.GetBinWidth(nbins)
        
    # --------------------------------------------------------
    # determine nature of smearing 
    # --------------------------------------------------------    
    dirname = qcdspectrum[0].dirname()
    smearing= getsmearing.findall(dirname)[0]
    
    if smearing == 'JEC':
        label  = 'JEC uncertainties only'
    elif smearing == 'PDF':
        label  = 'PDF uncertainties only'
    else:
        label  = 'JEC + PDF uncertainties'
    postfix = '_%s' % smearing

    # get PDF (directory) name
    prefix  = nameonly(context.filename)
    dirname = split(prefix, '_')[0]

    
    # --------------------------------------------------------
    print "="*80
    addCI = model != 'QCD'
    if addCI:
        postfix += '_%s_%3.1f_%s' % (model, Lambda, interf[0])
    else:
        postfix += '_QCD'

    # --------------------------------------------------------
    # get reference QCD spectrum
    # --------------------------------------------------------    
    pc = PercentileCurve(nbins)
    for index, QCD in enumerate(qcdspectrum[1:]):
        hqcd = QCD()
        hutil.divideByWidth(hqcd) # convert to a density
        pc.add(hqcd)

    hqcdnom = qcdspectrum[0]()
    hqcdnom.SetLineColor(kBlue)                
    hutil.divideByWidth(hqcdnom) # convert to a density
    
    if scaleToData:
        print "scale median prediction to data"
        c50 = pc(0.50)
        hqcdref = hdata.Clone('hqcdref')
        for ii in range(len(c50)):
            hqcdref.SetBinContent(ii+1, c50[ii])
            hqcdref.SetBinError(ii+1, 0)
    else:
        hqcdref = hqcdnom.Clone('hqcdref')

    xdata = hdata.Integral()
    xsect = hqcdref.Integral()
    scale = xdata / xsect
    print "\n==> data/theory: %10.2f\n" % scale
    hqcdnom.Scale(scale)

    # --------------------------------------------------------
    # compute percentile spectra and scale to data
    # --------------------------------------------------------    
    pcQCD   = PercentileCurve(nbins)
    pcQCDCI = PercentileCurve(nbins)
    
    hQCD    = []
    hQCDCI  = []
    
    for index, QCD in enumerate(qcdspectrum[1:]):
        hqcd = QCD()
        hqcd.Scale(scale)
        hutil.divideByWidth(hqcd) # convert to a density
        pcQCD.add(hqcd)
        hQCD.append(hqcd)

        if addCI:
            CI  = cispectrum[index]
            hci = CI(lam, kappa)
            hci.Scale(scale)            
            hutil.divideByWidth(hci) # convert to a density
            
            hci.Add(hqcd)
            pcQCDCI.add(hci)
            hQCDCI.append(hci)

        if index % 1000 == 0:
            print "%5d" % index, hqcd.GetName()
    
    # --------------------------------------------------------
    # plot spectra
    # --------------------------------------------------------
    os.system('mkdir -p figures/%s' % dirname)
    
    fname = 'figures/%s/%s_xsection%s' % (dirname, prefix, postfix)
    if scaleToData:
        fname = fname + '_scaled'
        
    print
    print fname
    cspect = TCanvas(fname, fname, 10, 10, 500, 500)
    x = map(lambda i: (pT[i+1]+pT[i])/2, range(nbins))

    # decide what spectrum is to be plotted
    if addCI:
        pc   = pcQCDCI
        hist = hQCDCI
    else:
        pc   = pcQCD
        hist = hQCD
        
    print
    n = hdata.GetNbinsX()
    print "data: bins %d, pT-range = (%-6.1f ... %-6.1f) GeV" % \
      (n, hdata.GetBinLowEdge(1), hdata.GetBinLowEdge(n)+hdata.GetBinWidth(n))
    
    n = hqcdref.GetNbinsX()
    print "QCD(reference): bins %d, pT-range = (%-6.1f ... %-6.1f) GeV" % \
      (n, hqcdref.GetBinLowEdge(1),
       hqcdref.GetBinLowEdge(n)+hqcdref.GetBinWidth(n))
       

    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95 = mkpline(x, curve[0], curve[-1], hdata, cspect, color=kGreen)
    p68 = mkpline(x, curve[1], curve[-2], hdata, cspect, color=kYellow)
    p50 = mkgraph(x, curve[2],
                  "Jet p_{T} (GeV)",
                  "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                  pTmin, pTmax,
                  ymin=YMIN,
                  ymax=YMAX,
                  color=kRed,
                  lwidth=1)
    p50.SetLineWidth(2)
    #h50 = p50.GetHistogram()
    
    cspect.cd()
    gPad.SetLogy()
    hdata.Draw("ep")
    p95.Draw('f same')
    p68.Draw('f same')
    p50.Draw('c same')
    hqcdnom.Draw('c same')
    hdata.Draw("ep same")

    energy = context.energy
    lumi   = context.lumi
    scribe = addTitle('CMS Preliminary  #surds=%sTeV L=%s/fb' % \
                      (energy, lumi),
                      0.035)
    scribe.vspace()
    if addCI:
        scribe.write("%s(#Lambda = %3.1f TeV) %s" % (model,
                                                     Lambda,
                                                     interf), 0.05)
    scribe.write(PDFS[dirname], 0.07)
    scribe.write(label, 0.11)

    xp = 0.70
    yp = 0.68
    xw = 0.18
    yw = 0.24
    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hdata, 'data', 'p')
    #lg.AddEntry(hqcdnom, 'QCD', 'l')
    lg.AddEntry(p50, 'median', 'l')
    lg.AddEntry(p68, '68%s' % '%', 'f')
    lg.AddEntry(p95, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cspect.Update()
    gSystem.ProcessEvents()
    cspect.SaveAs('.png')
    
    # --------------------------------------------------------
    # now plot the ratio
    # --------------------------------------------------------
    cratio = TCanvas('%s_ratio' % fname,
                     '%s - xsection-ratio' % dirname,
                     515, 10, 500, 500)

    hratio = hdata.Clone('hratio')
    hratio.Divide(hqcdnom)
    hratio.SetMinimum(YMIN)
    hratio.SetMaximum(YMAX)

    title = 'data / QCD_{NLO}#otimesNP#otimesEWK'
    hratio.GetYaxis().SetTitle(title)

    cratio.cd()
    hratio.Draw('ep')
    cratio.Update()
    gSystem.ProcessEvents()
    
    color = [kRed, kOrange, kGreen, kBlue, kMagenta]
    pc  = PercentileCurve(nbins)
    jj = 0
    for ii, h in enumerate(hist):
        h.Divide(hqcdnom)
        pc.add(h)
    hqcdnomr = hqcdnom.Clone('hqcdnomr')
    hqcdnomr.Divide(hqcdnom)
    hqcdnomr.SetLineWidth(2)
    hqcdnomr.SetLineColor(kBlue)
            
    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95r = mkpline(x, curve[0], curve[-1], hratio, cratio, color=kGreen)
    p68r = mkpline(x, curve[1], curve[-2], hratio, cratio, color=kYellow)
    p50r = mkgraph(x, curve[2],
                   "Jet p_{T} (GeV)", title,
                   pTmin, pTmax,
                   ymin=YMIN,
                   ymax=YMAX,
                   color=kRed,
                   lwidth=3)


    cratio.cd()
    gPad.SetLogy(kFALSE)
    hratio.Draw("ep")
    p95r.Draw('f same')
    p68r.Draw('f same')
    p50r.Draw('c same')
    hqcdnomr.Draw('c same')
    hratio.Draw("ep same")
    
    scriber = addTitle('CMS Preliminary  #surds=%sTeV L=%s/fb' % \
                       (energy, lumi),
                       0.035)
    scriber.vspace()
    if addCI:
        scriber.write("%s(#Lambda=%3.1f TeV) %s" % (model,
                                                    Lambda,
                                                    interf), 0.04)
    scriber.write(PDFS[dirname], 0.04)
    scriber.write(label, 0.04)
    scriber.write('data/theory = %3.2f' % scale, 0.24)

    xp = 0.24
    yp = 0.52
    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hratio, 'data', 'p')
    #lg.AddEntry(hqcdnomr, 'QCD', 'l')
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

    filename, Lambda, kappa, model, interf, scaleToData = decodeCommandLine()
    
    context = Context()
    context.filename = replace(filename, '_workspace', '')
    context.Lambda   = Lambda
    context.kappa    = kappa
    context.energy   = '13'
    lumi = LUMI/1000
    context.lumi     = '%5.2f' % lumi
    context.scaleToData= scaleToData
    
    # --------------------------------------------------------        
    print "\nloading workspace..."
    # --------------------------------------------------------    
    addCI = model != "QCD"
    qcdspectrum = []
    cispectrum  = []
    hfile = TFile(filename)
    ws = hfile.Get('CI')
    check(ws, "can't get workspace CI")

    # --------------------------------------------------------      
    # get data
    # --------------------------------------------------------
    hdata = ws.obj('hdata')
    check(hdata, "can't get hdata")
    nbins = hdata.GetNbinsX()
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hdata.SetMarkerSize(0.8)
    hdata.SetMarkerStyle(20)
    hutil.divideByWidth(hdata)
    hdata.Scale(1.0/LUMI)

    ymin = hdata.SetMinimum(YMIN)
    ymax = hdata.SetMaximum(YMAX)

    # get bin boundaries
    pt = hutil.binlowedges(hdata)
    pt.push_back(pt[-1]+hdata.GetBinWidth(nbins))
    pT = array('d')
    for i in range(len(pt)): pT.append(pt[i])
    context.pT = pT
    
    # --------------------------------------------------------      
    # get model
    # --------------------------------------------------------    
    pdf = ws.pdf('model')
    nspectra = pdf.size()
    print "number of spectra %d" % nspectra

    for ii in xrange(nspectra):
        qcd = pdf.QCD(ii)
        if qcd == None: break
        qcdspectrum.append(qcd)

        if addCI:
            ci = pdf.CI(ii)
            if ci == None: break
            cispectrum.append(ci)

        if ii % 1000 == 0:
            print "%4d %s %s" % (ii, qcd.dirname(), qcd.histname())

    # --------------------------------------------------------
    # now plot
    # --------------------------------------------------------
    context.hdata  = hdata
    context.interf = ''
    context.qcdspectrum = qcdspectrum
    context.cispectrum  = cispectrum

    models = split(model)
    print "\nplotting..."
    for model in models:
        context.model = model
         
        if context.model == 'QCD':
            makePlot(context)
        else:
            for ii in xrange(kappa.size()):
                kappa[ii] = KAPPA[model][ii]
            context.kappa = kappa
            context.interf= 'constructive'
            makePlot(context)
            
            for ii in xrange(kappa.size()):
                kappa[ii] =-KAPPA[model][ii]
            context.kappa = kappa
            context.interf= 'destructive'                
            makePlot(context)
            
    sleep(5)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
