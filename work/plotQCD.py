#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotQCD.py
# read QCD spectra and plot them
# created 16-Mar-2018 Harrison B. Prosper & Bipen Kotwal
#-----------------------------------------------------------------------------
import os, sys, re, argparse
from math import *
from histutil import *
from string import *
from array import array
from time import sleep
from glob import glob
from random import randint
from ROOT import *
#-----------------------------------------------------------------------------
gSystem.Load("../CI/lib/libCI.so")
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-8  # minimum and maximum y-limits for spectrum
YMAX    = 2.0
LUMI    = 35100.0  # 1/pb
DATAFILENAME  = '../data/data_13TeV_L035.1ifb.root'
class Context:
    pass    
#-----------------------------------------------------------------------------
PDFS = {'ALL'  : 'CT14nlo, MMHT2014nlo, NNPDF30_nlo',
        'CT14' : 'CT14nlo',
        'MMHT' : 'MMHT2014nlo',
        'NNPDF': 'NNPDF30_nlo'}

def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '16-Mar-2018'
    USAGE = '''
    python plotQCD.py [options] PDFdir (e.g., CT14)

    options
       -s<smearing>   pdf, jec, jecpdf [jecpdf]
       -m<member>     PDF member       [*]
       -r             rescale QCD prediction to data
    '''
    parser = argparse.ArgumentParser(usage=USAGE,
                                    version=VERSION)
    group  = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--smearing',
                        action='store',
                        default='jecpdf',
                        help='level of smearing to apply')
    
    group.add_argument('-m', '--member',
                        action='store',
                        default='*',
                        help='PDF member')

    group.add_argument('-r', '--rescale',
                        action='store_true',
                        help='scale QCD prediction to data')
        
    parser.add_argument('pdfdir')
    
    options = vars(parser.parse_args())
    member  = options['member']
    pdfdir  = options['pdfdir']
    rescale = options['rescale']
    smearing= upper(options['smearing'])
    dirname = '../fastNLO/%s' % pdfdir
    
    if member != '*': member = '%3.3d' % atoi(member)
    return (dirname, pdfdir, smearing, member, rescale)
#-----------------------------------------------------------------------------
def makePlot(context):
    model    = context.model
    smearing = context.smearing
    hdata    = context.hdata
    dirname  = context.dirname
    pdfdir   = context.pdfdir
    prefix   = pdfdir
    rescale  = context.rescale
    qcdspectrum = context.qcdspectrum
    pT       = context.pT
    
    print "\n\t<=== %s ===>" % context.model

 
    # --------------------------------------------------------
    # number of bins and range
    # --------------------------------------------------------    
    nbins = hdata.GetNbinsX()
    pTmin = pT[0]
    pTmax = pT[-1]
        
    print
    print "\t=> bins %d, pT-range = (%-6.1f ... %-6.1f) GeV" % \
      (nbins, pTmin, pTmax)

    # get reference QCD curve.
    hqcdnom = qcdspectrum[0]()
    hutil.divideByWidth(hqcdnom) # convert to a density
    hqcdnom.SetLineColor(kBlue)                
       
    xdata = hdata.Integral()
    xsect = hqcdnom.Integral()
    scale = xdata / xsect
    print xdata, xsect
    print "\n==> data/theory: %10.2f\n" % scale

    # scale QCD prediction to data if requested
    if rescale:
        print "\n\t<= rescaling QCD prediction to data =>\n"
        hqcdnom.Scale(scale)
        
    # --------------------------------------------------------
    # determine nature of smearing 
    # --------------------------------------------------------    
    if smearing == 'JEC':
        label  = 'JEC uncertainties only'
    elif smearing == 'PDF':
        label  = 'PDF uncertainties only'
    else:
        label  = 'JEC + PDF uncertainties'
    postfix = '_%s_QCD' % smearing

    # --------------------------------------------------------
    # accumulate histograms in order to compute percentile
    # curves
    # --------------------------------------------------------
    print "="*80
    pc   = PercentileCurve(nbins)
    hist = []
    for index, QCD in enumerate(qcdspectrum[1:]):
        hqcd = QCD()
        if rescale: hqcd.Scale(scale)
        hutil.divideByWidth(hqcd) # convert to a density
        pc.add(hqcd)
        hist.append(hqcd)
    
    # --------------------------------------------------------
    # plot spectra
    # --------------------------------------------------------
    os.system('mkdir -p figures/%s' % pdfdir)
    
    name = 'figures/%s/%s_xsection%s' % (pdfdir, prefix, postfix)
    cspect = TCanvas(name, name, 10, 10, 500, 500)

    # compute percentile spectra
    p50, p68, p95 = pc.plines(hdata, cspect,
                                "Jet p_{T} (GeV)",
                                "d^{2}#sigma /dp_{T}dy (pb/GeV)",
                                xmin=pTmin, xmax=pTmax)
    
    p50.SetLineColor(kRed)
    p50.SetLineWidth(2)
    
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
    scribe.write(PDFS[pdfdir], 0.07)
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
    gSystem.ProcessEvents()
    cspect.SaveAs('.png')
    # --------------------------------------------------------
    # now plot the ratio
    # --------------------------------------------------------
    cratio = TCanvas('figures/%s/%s_xsection%s_ratio' % (pdfdir,
                                                         prefix,
                                                             postfix),
                        '%s - xsection-ratio' % dirname,
                        515, 10, 500, 500)

    hratio = hdata.Clone('hratio')
    hratio.Divide(hqcdnom)
    hratio.SetMinimum(0.0)
    hratio.SetMaximum(3.0)

    ytitle = 'data / QCD_{NLO}#otimesNP#otimesEWK'
    hratio.GetYaxis().SetTitle(ytitle)

    cratio.cd()
    hratio.Draw('ep')
    cratio.Update()

    p50r, p68r, p95r = pc.plines(hratio, cratio,
                                     "Jet p_{T} (GeV)",
                                     ytitle,
                                     pTmin, pTmax,
                                     denom=hqcdnom)
    p50r.SetLineColor(kRed)
    p50r.SetLineWidth(2)
            
    qcdr = qcdspectrum[0]().Clone('qcdr')
    qcdr.Divide(hqcdnom)
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
    
    scriber = addTitle('CMS Preliminary  #surds=%sTeV L=%s/fb' % \
                       (energy, lumi),
                       0.035)
    scriber.vspace()
    scriber.write(PDFS[pdfdir], 0.04)
    scriber.write(label, 0.04)
    scriber.write('data / theory = %3.2f' % scale, 0.04)

    lg = mklegend(xp, yp, xw, yw)
    lg.AddEntry(hratio, 'data', 'p')
    lg.AddEntry(qcdr, 'QCD', 'l')
    lg.AddEntry(p50r, 'median', 'l')
    lg.AddEntry(p68r, '68%s' % '%', 'f')
    lg.AddEntry(p95r, '95%s' % '%', 'f')
    lg.Draw('same')
    
    cratio.Update()
    gSystem.ProcessEvents()
    cratio.SaveAs('.png')
    sleep(5)    
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotQCD.py ===>"
    setStyle()

    dirname, pdfdir, smearing, member, rescale = decodeCommandLine()
    
    context = Context()
    context.dirname  = dirname
    context.pdfdir   = pdfdir
    context.smearing = smearing
    context.member   = member
    context.rescale  = rescale
    context.model    = 'QCD'
    context.energy   = '13'
    lumi = LUMI/1000
    context.lumi     = '%5.2f' % lumi
    qcdspectrum      = []

    # --------------------------------------------------------      
    # get data
    # --------------------------------------------------------
    hfile = TFile(DATAFILENAME)
    hdata = hfile.Get('hdata')
    check(hdata, "can't get hdata")

    context.hdata = hdata

    print "\t=> number of jets: %d" % hdata.Integral()
    
    nbins   = hdata.GetNbinsX()
    hdata.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    hdata.GetYaxis().SetTitle('d^{2}#sigma /dp_{T}dy (pb/GeV)')
    hdata.SetMarkerSize(0.8)
    hdata.SetMarkerStyle(20)
    hutil.divideByWidth(hdata)
    hdata.Scale(1.0/LUMI)
    
    print "\t=> observed cross section: %5.2f pb" % hdata.Integral()
    
    ymin = hdata.SetMinimum(YMIN)
    ymax = hdata.SetMaximum(YMAX)
    
    pT = hutil.binlowedges(hdata)
    pT.push_back(pT[-1]+hdata.GetBinWidth(nbins))
    context.pT = pT    
    pTmin = pT[0]
    pTmax = pT[-1]
    print '\t=> pT-range: %10.0f %10.0f  GeV' % (pTmin, pTmax)
    
    # --------------------------------------------------------      
    # get QCD spectra
    # --------------------------------------------------------
    cmd = '%s/%s/%s/qcd.root' % (dirname, member, smearing)
    qcdfiles = glob(cmd)
    qcdfiles.sort()
    nspectra = len(qcdfiles)
    print "\t=> number of spectra %d" % nspectra
    print "\t=> first file: %s" % qcdfiles[0]
    print "\t=> last file:  %s" % qcdfiles[-1]

    # get directory names
    qcddirs = [os.path.dirname(x) for x in qcdfiles]

    # get histogram names
    hf = TFile(qcdfiles[0])
    histnames = hutil.histogramNames(hf)
    print
    for x in histnames:
        print x
    print
    hf.Close()
    
    # place nominal QCD spectrum first
    jj = 0
    name = "QCD%5.5d" % jj
    qcddir = qcddirs[jj]
    histname = 'nlo_1.000_1.000_000'
    qcdnom = QCDSpectrum(name, name, qcddir, histname, nbins, pTmin)        
    qcdspectrum.append(qcdnom)
    print
    print "\t=> %4d %s %s" % (jj, qcdnom.dirname(), qcdnom.histname())

    for ii in xrange(nspectra):
        # randomly pick a PDF member from the PDF set
        member = randint(1, nspectra-1)
        qcddir = qcddirs[member]

        # randomly pick a histogram (each corresponding
        # to different choices of renormalization and
        # factorization scales)
        jj = randint(0, len(histnames)-1)
        histname = histnames[jj]

        # create a QCD spectrum
        name = "QCD%5.5d" % jj
        qcd  = QCDSpectrum(name, name, qcddir, histname, nbins, pTmin)        
        qcdspectrum.append(qcd)
        if ii % 50 == 0:
            print "%4d %s %s" % (ii, qcd.dirname(), qcd.histname())

    # --------------------------------------------------------
    # now plot
    # --------------------------------------------------------
    context.qcdspectrum = qcdspectrum
    makePlot(context)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
