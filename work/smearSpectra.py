#!/usr/bin/env python
#-----------------------------------------------------------------------------
# smearSpectra.py
# Apply the jet response function given by the CMS Inclusive Jet Team
# (see 6-Sep-2013 presentation by Sanmay Ganguly to the Exotica Group) to the
# specified true spectra 
# created 09-Oct-2014 Harrison B. Prosper - a rewrite
#         23-Jun-2016 HBP update to Run II
#         14-Sep-2016 HBP change EWK correction histogram to histo.
#                         also place correction files in data directory
#-----------------------------------------------------------------------------
import os, sys, re, optparse
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from ROOT import gSystem, gPad, TH1D, TFile, TCanvas, kFALSE, kTRUE, kBlue
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
CIPATH  = os.environ['CIPATH']
YMIN    = 1.0e-10
YMAX    = 1.0e+8
#-----------------------------------------------------------------------------
# default list of histograms to be smeared
HISTNAMES = '''
nlo_0.500_0.500
nlo_0.500_1.000
nlo_1.000_0.500
nlo_1.000_1.000
nlo_1.000_2.000
nlo_2.000_1.000
nlo_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
DATAFILENAME  = '../data/data_13TeV_L035.1ifb.root'
JECUNFILENAME = '../data/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK8PFchs.root'
CORFILENAME   = '../data/data_13TeV_L035.1ifb_plus_corrections.root'
EWKCHISTNAME  = 'ak7/y_0.0-0.5/EwkCor'
NPCHISTNAME   = 'ak7/y_0.0-0.5/NPCor'
JERUNCERTAINTY= 0.1  # 10% relative uncertainty in jet energy resolution
#-----------------------------------------------------------------------------
LUMI = 35100.0 # 1/pb
#-----------------------------------------------------------------------------
def decodeCommandLine():
    VERSION = '16-Sep-2016'
    USAGE = '''
    python smearSpectra.py [options] <pathname to unsmeared spectra>

    options
       -s<smearing>   pdf, jec, jecpdf [jecpdf]
       -m<member>     PDF member       [*]
    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)
    parser.add_option('-s', '--smearing',
                      action='store',
                      dest='smearing',
                      type='string',
                      default='jecpdf',
                      help='level of smearing to apply')

    parser.add_option('-m', '--member',
                      action='store',
                      dest='member',
                      type='string',
                      default='*',
                      help='PDF member')    
    options, args = parser.parse_args()
    if len(args) == 0:
        print USAGE
        sys.exit(0)
    PDFdir = args[0]
    if PDFdir[-1] == '/': PDFdir=PDFdir[:-1]
    return (upper(options.smearing), PDFdir, options.member)
#-----------------------------------------------------------------------------
def makePlot(hname, spectrum, pt, COLOR=kBlue):
    x = array('d')
    nbins = pt.size()-1
    for ii in xrange(nbins+1): x.append(pt[ii])
    h = TH1D(hname, '', nbins, x)
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.GetYaxis().SetTitle('')
    h.SetLineColor(COLOR)
    h.SetLineWidth(1)
    for ii in xrange(nbins):
        # save cross section/bin
        h.SetBinContent(ii+1, spectrum(pt[ii], pt[ii+1]))
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def main():
    print "\n\t<=== smearSpectra.py ===>"
    
    smearingLevel, PDFdir, PDFmember = decodeCommandLine()
    doJEC = find(smearingLevel, 'J') > -1
    doPDF = find(smearingLevel, 'P') > -1

    cmd = '(?<=%s/)[0-9]+' % PDFdir
    getPDFmember = re.compile(cmd)

    gSystem.Load("libCI.so")
    from ROOT import hutil, JECUncertainty, JetSpectrum, JetSpectrumSmeared
    
    setStyle()

    if  doPDF and not doJEC:
        # --------------------------------------------------------------------
        # do PDF smearing only
        # --------------------------------------------------------------------
        member   = '*'
        histnames= HISTNAMES
        nsmears  = 1
        prefix   = 'PDF'
                
    elif doJEC and not doPDF:
        # --------------------------------------------------------------------
        # do JES and JER smearing only.
        # use PDF member zero, with mur=muf=1
        # and sample JES, JER 1000 times.
        # --------------------------------------------------------------------
        member   = '000'
        histnames= ['nlo_1.000_1.000']
        nsmears  = 1000
        prefix   = 'JEC'
        
    else:
        # --------------------------------------------------------------------
        # do PDF, JES and JER smearing.
        # use all PDF members and
        # for each sample JES and JER ONCE.
        # --------------------------------------------------------------------
        member   = PDFmember
        if member != '*': member = '%3.3d' % atoi(member)
        histnames= HISTNAMES
        nsmears  = 1
        prefix   = 'JECPDF'
                
    # ------------------------------------------------------------------------
    # get input rootfiles and construct output file names
    # ------------------------------------------------------------------------
    fdirs = glob('%s/%s' % (PDFdir, member))
    fdirs.sort()
    rootfiles= []
    for fdir in fdirs:
        rootfiles += glob('%s/*.root' % fdir)
        
    outfiles = []
    dirnames = {}
    for rootfile in rootfiles:
        t = split(rootfile, '/')
        filename = t[-1]
        dirname  = '%s/%s' % (joinfields(t[:-1], '/'), prefix)
        outfile  = '%s/%s' % (dirname, filename)
        if not os.path.exists(dirname):
            os.system('mkdir -p %s' % dirname)
        outfiles.append(outfile)
        
    # which is it, fastCI or fastNLO?
    fastNLO = find(rootfiles[0], 'fastNLO') > 0
    
    print "="*80
    print "\tfirst input root file:   %s" % rootfiles[0]
    print "\tlast  input root file:   %s" % rootfiles[-1]
    print "\tnumber of output files:  %d" % len(outfiles)
    print "\tnumber of smearings:     %d" % nsmears
    print "\thistograms to smear:     %s" % histnames
    
    if doJEC:
        print "\tinclude JES+JER uncertainties"
        
    if doPDF:
        print "\tinclude PDF uncertainties"
        
    if fastNLO:
        print "\tfastNLO"
    else:
        print "\tfastCI"
    print "="*80

    # --------------------------------------------------------
    # get correction functions
    # --------------------------------------------------------
    # get jet energy corrections
    if not os.path.exists(JECUNFILENAME):
        hutil.error('smearSpectra.py',
                    "can't open file %s" % JECUNFILENAME)
        
    JESunc = JECUncertainty(JECUNFILENAME)
    JERunc = JERUNCERTAINTY
    
    # get electroweak corrections
    if not os.path.exists(CORFILENAME):
        hutil.error('smearSpectra.py',
                    "can't open file %s" % CORFILENAME)        
    fCor = TFile(CORFILENAME)
    hEWK = fCor.Get(EWKCHISTNAME)
    if hEWK == None:
        hutil.error('smearSpectra.py',
                    "can't get histogram %s" % EWKCHISTNAME)

    # get non-perturbative corrections
    hNPC = fCor.Get(NPCHISTNAME)
    if hNPC == None:
        hutil.error('smearSpectra.py',
                    "can't gethistogram %s" % HPCHISTNAME)
        
    # --------------------------------------------------------
    # read normalvariates.
    # we do this so that we maintain the JES/JER correlation
    # between QCD and the CI cross sections. Any (PDF member,
    # mur, muf) triplet can go with any random choice of the
    # (JES, JER) doublets, but we want to make sure that the
    # same doublet is used for the same QCD and CI triplet.
    # we do this using the variate map (see below).
    # --------------------------------------------------------
    nxy = map(lambda x: map(atof, x),
              map(split, open("normalvariates.txt").readlines()))

    # --------------------------------------------------------      
    # read data
    # --------------------------------------------------------    
    hdfile = TFile(DATAFILENAME)
    if not hdfile.IsOpen():
        hutil.error('smearSpectra.py',
                    "can't open file %s" % DATAFILENAME)
    hdata = hdfile.Get('hdata')
    hdata.GetYaxis().SetTitle('#sigma / bin (pb)')
    hdata.Scale(1.0/LUMI)
    nbins = hdata.GetNbinsX()
    pT    = hutil.binlowedges(hdata)
    pT.push_back(pT.back()+hdata.GetBinWidth(nbins))

    # minimum and maximum pTs of smeared spectrum
    pTmin = pT[0]
    pTmax = pT[-1]
    print "\n\t==> bins = %4d,  "\
      "pT-range = (%-6.1f... %-6.1f) GeV\n" % (nbins, pTmin, pTmax)
    
    # --------------------------------------------------------
    # loop over files and smear selected histograms within
    # each file
    # --------------------------------------------------------              
    cspect = TCanvas('cspec', 'spectra', 10, 10, 500, 500)
    jxy = 0
    variates = {} # map between (PDF member, mur, muf) and (x, y)
    for index, rootfile in enumerate(rootfiles):
        if not os.path.exists(rootfile):
            hutil.error('smearSpectra.py',
                        "can't find rootfile %s" % rootfile)

        member = getPDFmember.findall(rootfile)
        if len(member) == 0:
            hutil.error('smearSpectra.py',
                        "can't get PDF member from %s" % \
                        rootfile)
        member = member[0]

        if index % 100 == 0:
            print "%5d\t%s" % (index, rootfile)
         
        # open an output file for smeared histograms
        hfile = TFile(outfiles[index], 'recreate')
        hist = []
        for histname in histnames:
            
            spectrum = JetSpectrum(rootfile, histname, fastNLO, hNPC, hEWK)
            
            for ii in xrange(nsmears):

                hname = '%s_%3.3d' % (histname, ii)
                
                if doJEC:
                    # apply jet energy correction smearing.
                    #
                    # get the correct mapping between
                    # (PDFmember, mur, muf) and (x, y)
                    key = '%s/%s/%3.3d' % (member, histname, ii)
                    if variates.has_key(key):
                        x, y = variates[key]
                    else:
                        x, y = nxy[jxy]
                        jxy += 1
                        variates[key] = (x, y)

                    sspectrum = JetSpectrumSmeared(spectrum,
                                                   JESunc,
                                                   JERunc,
                                                   x, y,
                                                       pTmin, pTmax)
                    hfile.cd()
                    # create histogram containing cross sections/bin
                    h = makePlot(hname, sspectrum, pT, kBlue)
                else:
                    # no jet energy correction smearing
                    hfile.cd()
                    h = makePlot(hname,  spectrum, pT, kBlue)
                    
                # cache histograms so that they aren't deleted
                hist.append(h)

                cspect.cd()
                if fastNLO:
                    gPad.SetLogy()
                    hdata.Draw('ep')
                    h.Draw('c same')
                else:
                    gPad.SetLogy(kFALSE)
                    h.Draw('c')
                cspect.Update()
            try:
                del spectrum
            except:
                pass
            try:
                del sspectrum
            except:
                pass
        hfile.Write()
        hfile.Close()
    sleep(5)
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
