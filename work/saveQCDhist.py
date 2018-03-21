#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File: saveQCDhist.py
# created 12-Oct-2014 Harrison B. Prosper - a rewrite
#         08-Feb-2015 HBP - add nominal QCD spectrum, that is, the
#                     spectrum computed with PDF member zero
#         16-May-2016 HBP - update to use random sampling implemented and
#                     tested by Roberto
#-----------------------------------------------------------------------------
import os, sys, re, optparse, histutil
from math import *
from string import *
from glob import glob
from array import array
from time import sleep
from random import shuffle, randint
from ROOT import gSystem, TFile, kFALSE, kTRUE, vector
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
CIPATH  = os.environ['CIPATH']
DATAFILE= '../data/data_13TeV_L000.07152ifb.root'
PDFDIR  = 'CT14'
NMEMBERS= 200 # number of PDF members/PDF set to use
NSAMPLE = 500
LUMI    = 72.0 # /pb
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '17-Sep-2016'
    USAGE = '''
    python saveQCDhists.py [options] [PDFset default=CT14 NNPDF MMHT]

    options
       -s<smearing>  pdf, jec, jecpdf          [jecpdf]
       -n<number of PDF members>               [%d]
       -b                                      [no boostrap] 
    ''' % NMEMBERS
    if len(sys.argv) < 2:
        print USAGE
        sys.exit()
        
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)
    parser.add_option('-s', '--smearing',
                      action='store',
                      dest='smearing',
                      type='string',
                      default='jecpdf',
                      help='level of smearing to use')

    parser.add_option('-n', '--nmembers',
                      action='store',
                      dest='nmembers',
                      type='int',
                      default=NMEMBERS,
                      help='number of PDF members to use/PDF set')
    
    parser.add_option('-b', '--bootstrap',
                      action='store_true',
                      dest='bootstrap',
                      default=True,
                      help='use a bootstrap sample of spectra')    
                                      
    options, PDFset = parser.parse_args()
    if len(PDFset) == 0:
        PDFset = 'CT14'
    else:
        PDFset = PDFset[0]

    # create directory name to contain smeared spectra
    directory = upper(options.smearing)
    directory = '/%s' % directory
        
    return (directory, PDFset, 
            options.nmembers, options.bootstrap, options.smearing)
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== saveQCDhists.py ==="

    # number of PDF members = number of subdirectories under each PDF
    # directory
    dirname, PDFset, ndirs, bootstrap, smearing=decodeCommandLine()
    
    gSystem.Load("../CI/lib/libCI.so")
    from ROOT import hutil, QCDSpectrum
    
    # -------------------------------------
    # determine which files to use
    # -------------------------------------
    if dirname == 'JEC':
        first = 0
        ndirs = 1
    else:
        first = 1
        ndirs+= 1
      
    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    rootfile = '../fastNLO/%s/000%s/qcd.root' % (PDFset, dirname)
    if not os.path.exists(rootfile):
        hutil.error("saveQCDhists.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("saveQCDhists.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()

    # pick only nlo histograms
    histnames = filter(lambda x: x[:3] == 'nlo', histnames)
    for ii, histname in enumerate(histnames):
        print "%2d\t%s" % (ii+1, histname)

    # -------------------------------------                     
    # get list of QCD directories, each
    # associated with a different PDF
    # member
    # -------------------------------------            
    QCDdirs = []
    for member in xrange(first, ndirs):
        filename = '../fastNLO/%s/%3.3d%s' % (PDFset, member, dirname)
        QCDdirs.append(filename)
    QCDdirs.sort()

    # -------------------------------------                
    # make a sample of histograms in which
    # the renormalization and factorization
    # scales are randomly selected.
    #
    # Each histogram corresponds to a
    # randomly selected PDF set, pair of
    # renormalization and factorization
    # scales,PDFs, jet energy scale, and
    # jet energy resolution.
    # -------------------------------------

    nspectra = len(QCDdirs)*len(histnames)
    spectra  = [None]*nspectra
    nsample  = 500

    if bootstrap:
        for jj in xrange(nspectra):

            # randomly pick a PDF member from the PDF set
            member = randint(0, len(QCDdirs)-1)        
            qcddir = QCDdirs[member]

            # randomly pick a histogram (each corresponding
            # to different choices of renormalization and
            # factorization scales)
            ii = randint(0, len(histnames)-1)
            histname = histnames[ii]
            spectra[jj] = (qcddir, histname)
    else:
        jj = 0
        for qcddir in QCDdirs:
            for histname in histnames:
                spectra[jj] = (qcddir, histname)
                jj += 1

    # pick first nsample spectra
    shuffle(spectra)
    spectra =spectra[:NSAMPLE]
    
    # insert nominal QCD histogram with nominal
    # scales at position 0 in list of spectra
    qcddir = '../fastNLO/%s/000%s' % (PDFset, dirname)
    fqcdnom = glob(qcddir)
    if len(fqcdnom) == 1:
        fqcdnom = fqcdnom[0]
    else:
        hutil.error("saveQCDhists.py",
                    "can't find directory %s" % qcddir)      
    histname = histnames[3] #(mur, muf=1, 1)
    spectra.insert(0, (qcddir, histname))
    nspectra = len(spectra)
    
    print "="*80
    print "==> nominal spectrum:  %s, %s" % (spectra[0])
    print "==> number of spectra: %d" % nspectra,
    if bootstrap:
        print "\t(use bootstrap sample)"
    else:
        print "\t(use original sample)"

    # --------------------------------------
    # open output histogram file
    # --------------------------------------
    outfilename = 'qcd_%s.root' % PDFset
       
    # -------------------------------------            
    # get data set
    # -------------------------------------
    hdfile = TFile(DATAFILE)
    hdata  = hdfile.Get('hdata')
    hdata.GetXaxis().SetTitle('Jet #font[12]{p}_{T} (GeV)')
    hdata.GetYaxis().SetTitle('count / bin')
    nbins = hdata.GetNbinsX()
    ptlow = hdata.GetBinLowEdge(1)

    hdout = TFile('qcd_data.root', 'recreate')
    hdata.Write()
    hdout.Close()
    
    print "==> saving %d spectra to file %s" % (nspectra, outfilename)
    hfile = TFile(outfilename, 'recreate')
    
    qcdspectrum = []
    qcdrecords  = []

    for index, (QCDdir, histname) in enumerate(spectra):

        if index % 100 == 0:
            print "%5d %s\t%s" % (index, QCDdir, histname)

        name = "hqcd%3.3d" % index
        QCD  = QCDSpectrum(name, name, QCDdir, histname, nbins, ptlow)
        hqcd = QCD()
        hqcd.SetName(name)
        hqcd.GetXaxis().SetTitle('Jet #font[12]{p}_{T} (GeV)')        
        hqcd.GetYaxis().SetTitle('#sigma_{i} (pb)')
        hfile.cd()
        hqcd.Write()

        if index == 0:
            out = open("data_vs_qcd_%s.txt" % PDFset, 'w')
            rec = '%10s %10s %10s %10s' % ('pTlow', 'pTup', 'observed', 'predicted')
            out.write('%s\n' % rec)
            for i in xrange(nbins):
                ptlow = hdata.GetBinLowEdge(i+1)
                ptup  = ptlow + hdata.GetBinWidth(i+1)
                obs   = hdata.GetBinContent(i+1)
                pred  = hqcd.GetBinContent(i+1) * LUMI
                rec = '%10.0f %10.0f %10.0f %10.1f' % \
                  (ptlow, ptup, obs, pred)
                out.write('%s\n' % rec)                
            out.close()
    hfile.Close()

    print "\tdone!\n"    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
