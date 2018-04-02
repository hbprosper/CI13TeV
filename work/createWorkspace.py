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
     RooWorkspace, RooMsgService, RooFit, RooDataSet, RooCmdArg, vector
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
CIPATH  = os.environ['CIPATH']
YMIN    = 1.0e-8  # minimum and maximum y-limits for spectrum
YMAX    = 2.0
LUMI    = 35100.0  # 1/pb
DATAFILE= '../data/data_13TeV_L035.1ifb.root'
PDFDIR  = 'CT14'
NMEMBERS= 200 # number of PDF members/PDF set to use
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '17-Sep-2016'
    USAGE = '''
    python createWorkspace.py [options] [PDF=CT14 NNPDF MMHT]

    options
       -s<smearing>  pdf, jec, jecpdf          [jecpdf]
       -o<root filename>                       [<PDF>_<s>_workspace.root]
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

    parser.add_option('-o', '--output',
                      action='store',
                      dest='filename',
                      type='string',
                      default='',
                      help='output file containing QCD and CI spectra')    

    parser.add_option('-n', '--nmembers',
                      action='store',
                      dest='nmembers',
                      type='int',
                      default=NMEMBERS,
                      help='number of PDF members to use/PDF set')

    parser.add_option('-b', '--bootstrap',
                      action='store_true',
                      dest='bootstrap',
                      help='use a bootstrap sample of spectra')    
                                      
    options, PDFsets = parser.parse_args()
    if len(PDFsets) == 0:
        PDFsets = ['CT14', 'MMHT', 'NNPDF']

    # create directory name to contain smeared spectra
    directory = upper(options.smearing)

    # create name of output file
    filename = options.filename
    if filename == '':
        if len(PDFsets) > 1:
            prefix = 'ALL'
        else:
            prefix = PDFsets[0]
        filename = '%s_%s_workspace.root' % (prefix, directory)

    return (directory, PDFsets, filename,
            options.nmembers, options.bootstrap, options.smearing)
#-----------------------------------------------------------------------------
def main():
    print "\n\t\t\t=== createWorkspace.py ==="

    # number of PDF members = number of subdirectories under each PDF
    # directory
    dirname, PDFsets, wfilename, ndirs, bootstrap, smearing=decodeCommandLine()
    
    gSystem.Load("libCI")
    from ROOT import hutil, QCDSpectrum, CIXsection, CISpectrum, \
     RooInclusiveJetPdf
    
    # -------------------------------------
    # determine which files to use
    # -------------------------------------
    print "dirname", dirname
    if dirname == 'JEC':
        first = 0
        ndirs = 1
    else:
        first = 1
        ndirs+= 1
        
    dirname = '/%s' % dirname
    
    # -------------------------------------                
    # get list of histograms
    # -------------------------------------            
    rootfile = '../fastNLO/%s/000%s/qcd.root' % (PDFDIR, dirname)
    if not os.path.exists(rootfile):
        hutil.error("createWorkspace.py",
                    "can't find rootfile %s" % rootfile)
    hfile = TFile(rootfile)
    if not hfile.IsOpen():
        hutil.error("createWorkspace.py",
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
    print "ndirs", ndirs
    QCDdirs = []
    for pdfset in PDFsets:
        for member in xrange(first, ndirs):
            filename = '../fastNLO/%s/%3.3d%s' % (pdfset, member, dirname)
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
    spectra = [None]*nspectra

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
                
    # insert nominal QCD histogram for CT14 with nominal
    # scales at position 0 in list of spectra
    qcddir = '../fastNLO/%s/000%s' % (PDFDIR, dirname)
    fqcdnom = glob(qcddir)
    if len(fqcdnom) == 1:
        fqcdnom = fqcdnom[0]
    else:
        hutil.error("createWorkspace.py",
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
    # create RooFit workspace. the workspace
    # assumes headers are in include and
    # sources are in src
    # --------------------------------------
    print "==> create workspace..."
    RooWorkspace.autoImportClassCode(kTRUE)
    RooWorkspace.addClassDeclImportDir('../CI')
    RooWorkspace.addClassImplImportDir('../CI')
    
    ws = RooWorkspace("CI")

    # Suppress info messages
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)    

    # -------------------------------------
    # create model parameters
    # Note: l = 1/Lambda^2
    # -------------------------------------    
    ws.factory('lambda[0, 0, 0.02]')
    l = ws.var("lambda")
    l.setBins(200)    
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
    ptlow = hdata.GetBinLowEdge(1)
    
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
        if index % 100 == 0:
            print "%5d %s\t%s" % (index, QCDdir, histname)

        name = "QCD%5.5d" % index
        qcdrecords.append(name)
        qcdspectrum.append(QCDSpectrum(name, name, QCDdir,
                                       histname, nbins, ptlow))
        if index == 0:
            h  = qcdspectrum[-1]()
            nb = h.GetNbinsX()
            print "\tQCD spectrum check: %d, %6.1f, %6.1f" % \
              (nb,
               h.GetBinLowEdge(1),
               h.GetBinLowEdge(nb)+h.GetBinWidth(nb))
            
        name = "CI%5.5d" % index
        cirecords.append(name)
        cispectrum.append(CISpectrum(name, name, CIdir,
                                        histname, nbins, ptlow))
        if index == 0:
            k = vector('double')(6,0)
            k[0]=-1        
            h  = cispectrum[-1](0, k)
            nb = h.GetNbinsX()
            print "\tCI  spectrum check: %d, %6.1f, %6.1f" % \
              (nb,
               h.GetBinLowEdge(1),
               h.GetBinLowEdge(nb)+h.GetBinWidth(nb))

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


    # check calculation of Asimov dataset
    model.setAsimov(True, LUMI)
    print "\n\t== Use Asimov data set (Lumi = %10.1f / pb)" % LUMI
    Asimov = model.Asimov()
    for ii in xrange(Asimov.size()):
        print "\t%4d\t%10.1f" % (ii+1, Asimov[ii])
         
    print "\tdone!\n"    
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
