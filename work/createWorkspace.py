#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File: createWorkspace.py
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
from ROOT import gSystem, TFile, kFALSE, kTRUE, \
     RooWorkspace, RooMsgService, RooFit, RooDataSet, RooCmdArg
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
DATAFILE= '../data/data_13TeV_L000.07152ifb.root'
PDFDIR  = 'CT14'
NDIRS   = 200
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '16-May-2016'
    USAGE = '''
    python createWorkspace.py [options] [PDF sets=CT14 NNPDF MMHT]

    options
       -s<smearing>  jes+jer+pdf, jes+jer, pdf [def.=jes+jer+pdf]
       -o<root filename>                       [def.=<PDF>_<s>_workspace.root]
       -n<number of directories>               [def.=%d]
       -h (help)
    ''' % NDIRS
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)
    parser.add_option('-s', '--smearing',
                      action='store',
                      dest='smearing',
                      type='string',
                      default='jes+jer+pdf',
                      help='level of smearing to use')

    parser.add_option('-o', '--output',
                      action='store',
                      dest='filename',
                      type='string',
                      default='',
                      help='output file containing QCD and CI spectra')    

    parser.add_option('-n', '--nspectra',
                      action='store',
                      dest='nspectra',
                      type='int',
                      default=NDIRS,
                      help='number of randomly sampled spectra to use'\
                      ' MC bootstrap integration ')
                                  
    options, PDFsets = parser.parse_args()
    if len(PDFsets) == 0:
        PDFsets = ['CT14', 'MMHT', 'NNPDF']

    # create directory name to contain smeared spectra
    directory = upper(replace(options.smearing, '+', ''))

    # create name of output file
    filename = options.filename
    if filename == '':
        if len(PDFsets) > 1:
            prefix = 'ALL'
        else:
            prefix = PDFsets[0]
        filename = '%s_%s_workspace.root' % (prefix, directory)
        
    return (directory, PDFsets, filename, options.nspectra)
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== createWorkspace.py ==="
    
    dirname, PDFsets, wfilename, ndirs = decodeCommandLine()
    
    print
    print dirname, PDFsets, wfilename, ndirs

    gSystem.Load("libCI")
    from ROOT import hutil, QCDSpectrum, CIXsection, CISpectrum, \
     RooInclusiveJetPdf
    
    # -------------------------------------
    # determine which files to use
    # -------------------------------------
    if dirname == 'JESJER':
        first = 0
        ndirs = 1
    else:
        first = 1
        ndirs+= 1

    # -------------------------------------
    # create a random (bootstrap) sample
    # of spectra that will serve as the
    # sample to approximate the integral of
    # the likelihood over the JES, JER, and
    # PDF nuisance parameters. Use the
    # bootstrap algorithm Roberto implemented
    # since he has shown that it works well.                
    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    rootfile = '../fastNLO/%s/000/%s/qcd.root' % (PDFDIR, dirname)
    if not os.path.exists(rootfile):
        hutil.error("createWorkspace.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("createWorkspace.py",
                    "can't open rootfile %s" % rootfile)
    histnames = hutil.histogramNames(hfile)
    hfile.Close()
    for ii, histname in enumerate(histnames):
        print "%2d\t%s" % (ii+1, histname)
        
    # get list of QCD directories containing spectra
    QCDdirs = []
    for pdfset in PDFsets:
        for member in xrange(first, ndirs):
            filename = '../fastNLO/%s/%3.3d/%s' % (pdfset, member, dirname)
            QCDdirs.append(filename)
    QCDdirs.sort()

    # -------------------------------------                
    # make a sample of histograms in which
    # the renormalization and factorization
    # scales are randomly selected.
    # Each histogram corresponds to a
    # randomly selected PDF set, pair of
    # renormalization and factorization
    # scales,PDFs, jet energy scale, and
    # jet energy resolution.
    # -------------------------------------
    print
    nspectra = len(QCDdirs)*len(histnames)
    print "==> create a bootstrap sample of %d spectra..." % nspectra
    spectra = [None]*nspectra
    
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

    # insert nominal QCD histogram for CT14 with nominal
    # scales at position 0 in list of spectra
    qcddir = '../fastNLO/%s/000/%s' % (PDFDIR, dirname)
    fqcdnom = glob(qcddir)
    if len(fqcdnom) == 1:
        fqcdnom = fqcdnom[0]
    else:
        hutil.error("createWorkspace.py",
                    "can't find directory %s" % qcddir)      
    histname = 'nlo_1.000_1.000_000'
    spectra.insert(0, (qcddir, histname))

    print "="*80
    print "==> number of spectra: %d" % len(spectra)
    
    # --------------------------------------
    # create RooFit workspace
    # --------------------------------------
    print "==> create workspace..."
    RooWorkspace.autoImportClassCode(kTRUE)
    RooWorkspace.addClassDeclImportDir('../CI/include')
    RooWorkspace.addClassImplImportDir('../CI/src')
    
    ws = RooWorkspace("CI")

    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    

    # -------------------------------------
    # create model parameters
    # Note: l = 1/Lambda^2
    # -------------------------------------    
    ws.factory('lambda[0, 0, 0.015]')
    l = ws.var("lambda")
    l.setBins(150)    
    # note upper case "S" in SetTitle; yes this is annoying!
    l.SetTitle('#font[12]{#lambda} = 1/#Lambda^{2} (TeV^{-2})')

    # create kapppa parameters
    record = []
    for ii in xrange(6):
        name = 'kappa%d' % ii
        ws.factory('%s[0,-1,1]' % name)
        record.append(name)
    ws.var('kappa0').setVal(-1) # default: LL model
    record = joinfields(record, ',')
    ws.defineSet('kappaset', record)
        
    # -------------------------------------            
    # create data set
    # -------------------------------------
    hdfile = TFile(DATAFILE)
    hdata  = hdfile.Get('hdata')
    hdata.GetXaxis().SetTitle('Jet #font[12]{p}_{T} (GeV)')
    hdata.GetYaxis().SetTitle('count / bin')
    nbins = hdata.GetNbinsX()
    getattr(ws, 'import')(hdata, 'hdata')

    # cretae a RooFit parameter for each count
    record = [] 
    for ii in xrange(nbins):
        c = hdata.GetBinContent(ii+1)
        j = int(c + 10*sqrt(c))/10000
        cmax = (j+1)*10000           
        name = 'N%2.2d' % ii
        variable = '%s[%10.0f, 0,%10.0f]' % (name, c, cmax)
        ws.factory(variable)
        record.append (name)
        print "%5.0f\t%s" % (hdata.GetBinLowEdge(ii+1), variable)
        record.append (name)
    record = joinfields(record, ',')
    ws.defineSet('Nset', record)    

    # make a RooDataSet from the counts
    data = RooDataSet('data',
                      'counts per bin',
                      ws.set('Nset'))
    data.add(ws.set('Nset'))
    getattr(ws, 'import')(data, RooCmdArg())
        
    # -------------------------------------            
    # finally, create probability model
    # -------------------------------------
    model = RooInclusiveJetPdf('model', '#font[12]{p}(D|#fond[12]{#lambda}, '\
                               '#fonf[12]{#kappa})',
                               ws.set('Nset'),
                               ws.var('lambda'),
                               ws.set('kappaset'))

    print "==> saving %d spectra to workspace %s" % (nspectra, wfilename)

    hfile = TFile(wfilename, "recreate")
    qcdspectrum = []
    cispectrum  = []
    qcdrecords  = []
    cirecords   = []
    for index, (QCDdir, histname) in enumerate(spectra):
        # get CI directory associated with current QCD spectrum
        CIdir = replace(QCDdir, 'fastNLO', 'fastCI')
        if index % 500 == 0:
            print "%5d %s\t%s" % (index, QCDdir, histname)

        name = "QCD%5.5d" % index
        qcdrecords.append(name)
        qcdspectrum.append(QCDSpectrum(name, name, QCDdir, histname))

        name = "CI%5.5d" % index
        cirecords.append(name)
        cispectrum.append( CISpectrum(name, name, CIdir, histname) )

        # also add spectra to probability model
        model.add(qcdspectrum[-1], cispectrum[-1])

    getattr(ws, 'import')(model, RooCmdArg())
    # -------------------------------------                
    print "="*80    
    ws.Print()
    print "="*80    
    print "==> writing workspace to file: %s" % wfilename
    hfile.Close()
    ws.writeToFile(wfilename, kFALSE)
    print "\tdone!\n"    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
