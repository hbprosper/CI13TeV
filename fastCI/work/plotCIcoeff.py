#!/usr/bin/env python
#-----------------------------------------------------------------------------
# plotCICoeff.py
# plot CI coefficient spectra
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#-----------------------------------------------------------------------------
import os, sys, re
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from random import shuffle
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#-----------------------------------------------------------------------------
gSystem.Load("libCI.so")
getname = re.compile('(?<=[0-9]/)[ab].+(?=[._])')
SQRTS = 13
COEFF = 'bi42'
#-----------------------------------------------------------------------------
def makePlot(hname, spectrum, pt, COLOR=kBlue):
    x = array('d')
    nbins = pt.size()-1
    for ii in xrange(nbins+1): x.append(pt[ii])
    h = TH1D(hname, '', nbins, x)
    h.GetXaxis().SetTitle('p_{T} (GeV)')
    h.GetYaxis().SetTitle('EWK #otimes NP #otimes '\
                          '#d^{2}#sigma /dp_{T}dy (pb/GeV)')
    h.SetLineColor(COLOR)
    for ii in xrange(nbins):
        pT = (pt[ii+1]+pt[ii])/2
        h.SetBinContent(ii+1, spectrum(pT))
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== plotCIcoeff.py ===>"
    setStyle()

    # --------------------------------------------------------
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        sys.exit('''
    ./plotCIcoeff.py PDFset [PDFmember=all] [smearing-dir=none]
        ''')
    PDFset = argv[0]
    os.system('mkdir -p figs/%s' % PDFset)

    if argc > 1:
        PDFmin = atoi(argv[1])
        PDFmax = PDFmin
    else:
        PDFmin = 0
        PDFmax = 200
    
    if argc > 2:
        smearDir = '%s/' % argv[2]
    else:
        smearDir = ''


    # --------------------------------------------------------
    print "="*80
    print "\tPDFset:        %s" % PDFset
    print "\tPDFmember:     %d" % PDFmin
    print "\tsmearingDir:   %s" % smearDir
    print "="*80
    # --------------------------------------------------------
    
    tdir  = TDirectory('CI', 'C')
    hist  = []
    kolor = [kRed, kOrange, kYellow+1, kGreen+1, kBlue, kMagenta]
    jcolor=0       
    for index in xrange(PDFmin, PDFmax+1):
        dirname = '../%s/%3.3d/%s' % (PDFset, index, smearDir)
        print dirname

        rootfiles = glob('%s/*.root' % dirname)
        # ----------------------------------------------------
        # sort files
        # ----------------------------------------------------
        filelist = []
        for kk, rootfile in enumerate(rootfiles):
            if not os.path.exists(rootfile):
                hutil.error("plotCIcoeff.py",
                            "can't find rootfile %s" % rootfile)
            hfile = TFile(rootfile)
            if not hfile.IsOpen():
                hutil.error("plotCIcoeff.py",
                            "can't open rootfile %s" % rootfile)
            name = getname.findall(rootfile)[0]
            filelist.append(('%4s' % name, rootfile))
        filelist.sort()

        # ----------------------------------------------------
        # plot
        # ----------------------------------------------------             
        figname = 'figs/%s/coeff_%3.3d_%dTeV' %  (PDFset, index, SQRTS)
        canvas = TCanvas(figname, figname, 10, 10, 1000, 800)
        canvas.Divide(6,10)

        coeff = TCanvas('coeff', 'coeff', 1100, 10, 500, 500)
        
        lastname = ''
        for kk, (name, rootfile) in enumerate(filelist):
            if lastname in ['ai5',  'ai43', 'aij8',
                            'aig5', 'ai4g3','aijg8',
                            'bi5',  'bi43', 'bij8']:
                jcolor += 1
                if jcolor >= len(kolor): jcolor = 0
            lastname = strip(name)
            color = kolor[jcolor]
            # get histograms from current file
            if not os.path.exists(rootfile):
                hutil.error("plotCIcoeff.py",
                            "can't find rootfile %s" % rootfile)
            hfile = TFile(rootfile)
            if not hfile.IsOpen():
                hutil.error("plotCIcoeff.py",
                            "can't open rootfile %s" % rootfile)
            
            histnames = filter(lambda x: x[:3] == 'nlo',
                               hutil.histogramNames(hfile))
            canvas.cd(kk+1)
            option = 'l'
            for histname in histnames:
                h = hfile.Get(histname)
                fixhist2(h)
                tdir.cd()
                h = h.Clone('%s%d' % (histname, kk))
                h.SetLineColor(color)
                scribe = Scribe(0.5, 0.6, 0.24)
                h.Draw(option)
                scribe.write(name)
                hist.append(h)
                option = 'l same'
            canvas.Update()

            coeff.cd()
            option = 'l'
            for histname in histnames:
                h = hfile.Get(histname)
                fixhist2(h)
                tdir.cd()
                h = h.Clone('%s%d-bi42' % (histname, kk))
                h.SetLineColor(color)
                scribe = Scribe(0.6, 0.7, 0.04)
                h.Draw(option)
                scribe.write(name)
                scribe.write('#sqrt{s} = %dTeV' % SQRTS)
                scribe.write('PDF = %s, %d' % (PDFset, PDFmin))
                hist.append(h)
                option = 'l same'
            coeff.Update()
            gSystem.ProcessEvents()
            
            figname = 'figs/%s/coeff_%3.3d_%dTeV_%s.png' %  (PDFset, index,
                                                             SQRTS, name)
            coeff.SaveAs(figname)
            sleep(1)
            hfile.Close()
        canvas.SaveAs('.png')
          
    sleep(2)       
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
