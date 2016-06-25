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
gSystem.Load("libCI.so")
from ROOT import hutil
    
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#                        kappa
#                   1   2   3   4   5   6
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']

LUMI    = 71.5  # 1/pb

class Context:
    pass    
#-----------------------------------------------------------------------------
getpdfdir   = re.compile('(?<=fastNLO/).+(?=/[0-9])|(?<=fastCI/).+(?=/[0-9])')
getsmearing = re.compile('(?<=/[0-9][0-9][0-9]/).+')
PDFS = {'ALL'  : 'CT14, MMHT2014, NNPDF30',
        'CT14' : 'CT14',
        'MMHT' : 'MMHT2014',
        'NNPDF': 'NNPDF30'}
#-----------------------------------------------------------------    
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
       -r<reference>  median, nominal (QCD)        [nominal]
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
                      default='QCD LL RR VV AA V-A',
                      help='model')
        
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

    return (filename, Lambda, kappa, model, interf)    
#-----------------------------------------------------------------------------
def makePlot(context):
    model    = context.model
    Lambda   = context.Lambda
    kappa    = context.kappa
    interf   = context.interf
    hdata    = context.hdata
    pT       = context.pT
    
    qcdspectrum = context.qcdspectrum
    cispectrum  = context.cispectrum
    
    print "\n\t<=== %s ===> %s" % (context.model, context.interf)
    
    lam   = 1.0/Lambda**2
    pTmin = pT[0]
    pTmax = pT[-1]
    nbins = hdata.GetNbinsX()
    # --------------------------------------------------------
    smearing = getsmearing.findall(qcdspectrum[0].dirname())[0]
    dirname  = split(context.filename, '_')[0]

    if   smearing == 'NONE':
        label = 'no smearing'
        postfix = '_NONE'
    elif smearing == 'JESJER':
        label  = 'JES+JER uncert.'
        postfix= '_JESJER'         
    elif smearing == 'PDF':
        label  = 'PDF uncert.'
        postfix= '_PDF'
    else:
        label  = 'JES+JER+PDF uncert.'
        postfix= '_JESJERPDF'
    
    # --------------------------------------------------------
    print "="*80
    addCI = model != 'QCD'
    if addCI:
        postfix += '_%s_%3.1f_%s' % (model, Lambda, interf[0])
    else:
        postfix += '_QCD'
    
    pcQCD   = PercentileCurve(nbins)
    pcQCDCI = PercentileCurve(nbins)
    hQCD    = []
    hQCDCI  = []
    
    for index, QCD in enumerate(qcdspectrum[1:]):
        hqcd = QCD()
        hutil.divideByWidth(hqcd) # convert to a density
        pcQCD.add(hqcd)
        hQCD.append(hqcd)

        if addCI:
            CI  = cispectrum[index]
            hci = CI(lam, kappa)
            hutil.divideByWidth(hci) # convert to a density
            hci.Add(hqcd)
            pcQCDCI.add(hci)
            hQCDCI.append(hci)

        if index % 1000 == 0:
            print "%5d" % index, hqcd.GetName()
    
    # --------------------------------------------------------
    # plot spectrum
    # --------------------------------------------------------
    os.system('mkdir -p figures/%s' % dirname)
    
    name = 'figures/%s/%s_d2sigma_dpTdy%s' % (dirname,
                                              dirname, postfix)
    cspect = TCanvas(name, name, 10, 10, 500, 500)
    x = map(lambda i: (pT[i+1]+pT[i])/2, range(nbins))

    # decide what spectrum is to be plotted
    if addCI:
        pc   = pcQCDCI
        hist = hQCDCI
    else:
        pc   = pcQCD
        hist = hQCD
        
    # get reference QCD curve.
    hqcdnom = qcdspectrum[0]()
    hutil.divideByWidth(hqcdnom) # convert to a density
    hqcdnom.SetLineColor(kBlue)            
                
    xdata = hdata.Integral()
    xsect = hqcdnom.Integral()
    scale = xdata / xsect
    print "\n==> data/theory: %10.3f\n" % scale
    
    # compute percentile spectra
    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95 = mkpline(x, curve[0], curve[-1], hdata, color=kGreen)
    p68 = mkpline(x, curve[1], curve[-2], hdata, color=kYellow)
    p50 = mkgraph(x, curve[2],
                  "Jet p_{T} (GeV)",
                  "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                  pTmin, pTmax,
                  ymin=YMIN,
                  ymax=YMAX,
                  color=kRed,
                  lwidth=1)
    p50.SetLineWidth(2)
    h50 = p50.GetHistogram()
    
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
    scribe = addTitle('CMS Preliminary  #surds=%sTeV L=%s/pb' % \
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
    lg.AddEntry(hqcdnom, 'QCD', 'l')
    lg.AddEntry(p50, 'median', 'l')
    lg.AddEntry(p68, '68%s' % '%', 'f')
    lg.AddEntry(p95, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cspect.Update()
    cspect.SaveAs('.png')

    # --------------------------------------------------------
    # now plot the ratio
    # --------------------------------------------------------
    cratio = TCanvas('figures/%s/%s_data_over_theory%s' % (dirname,
                                                           dirname, postfix),
                     '%s - spectrum-ratio' % dirname,
                     515, 10, 500, 500)
   
    hqcdnom.Scale(scale)
    hratio = hdata.Clone('hratio')
    hratio.Divide(hqcdnom)
    hratio.SetMinimum(0.0)
    hratio.SetMaximum(3.0)

    title = 'data / QCD_{NLO}#otimesNP#otimesEWK'
    hratio.GetYaxis().SetTitle(title)

    cratio.cd()
    hratio.Draw('ep')
    cratio.Update()

    color = [kRed, kOrange, kGreen, kBlue, kMagenta]
    pc  = PercentileCurve(nbins)
    jj = 0
    for ii, h in enumerate(hist):
        h.Divide(hqcdnom)
        pc.add(h)

    curve = []
    for p in PERCENT: curve.append( pc(p) )
    p95r = mkpline(x, curve[0], curve[-1], hratio, color=kGreen)
    p68r = mkpline(x, curve[1], curve[-2], hratio, color=kYellow)
    p50r = mkgraph(x, curve[2],
                   "Jet p_{T} (GeV)", title,
                   pTmin, pTmax,
                   ymin=0,
                   ymax=3,
                   color=kRed,
                   lwidth=1)

    xx = array('d'); xx.append(pTmin); xx.append(pTmax)
    yy = array('d'); yy.append(1); yy.append(1)
    qcdr = TGraph(2, xx, yy)
    qcdr.SetLineWidth(2)
    qcdr.SetLineColor(kBlue)
        
    cratio.cd()
    gPad.SetLogy(kFALSE)
    hratio.Draw("ep")
    p95r.Draw('f same')
    p68r.Draw('f same')
    p50r.Draw('c same')
    qcdr.Draw('c same')
    hratio.Draw("ep same")
    
    scriber = addTitle('CMS Preliminary  #surds=%sTeV L=%s/pb' % \
                       (energy, lumi),
                       0.035)
    scriber.vspace()
    if addCI:
        scriber.write("%s(#Lambda=%3.1f TeV) %s" % (model,
                                                    Lambda,
                                                    interf), 0.04)
    scriber.write(PDFS[dirname], 0.04)
    scriber.write(label, 0.04)

    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hratio, 'data', 'p')
    lg.AddEntry(qcdr, 'QCD', 'l')
    lg.AddEntry(p50r, 'median', 'l')
    lg.AddEntry(p68r, '68%s' % '%', 'f')
    lg.AddEntry(p95r, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cratio.Update()
    cratio.SaveAs('.png')    

    sleep(10)           
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
    context.lumi   = '%5.2f' % LUMI
        
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
    nbins   = hdata.GetNbinsX()
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hdata.SetMarkerSize(0.8)
    hdata.SetMarkerStyle(20)
    hutil.divideByWidth(hdata)
    hdata.Scale(1.0/LUMI)
    pT = hutil.binlowedges(hdata)
    pT.push_back(pT[-1]+hdata.GetBinWidth(nbins))

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
    context.hdata = hdata
    context.interf= ''    
    context.pT = pT
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

            waitForKey("\n**hit any key to continue ")
            
            for ii in xrange(kappa.size()):
                kappa[ii] =-KAPPA[model][ii]
            context.kappa = kappa
            context.interf= 'destructive'                
            makePlot(context)
            sleep(4)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'