#!/usr/bin/env python
#-----------------------------------------------------------------------------
# create inclusive jet CI histograms
# created 13-Oct-2014 Harrison B. Prosper
# Updated 17-Sep-2016 HBP - rename from mkCIhists.py to createCIhists.py
#                         - get CIPATH, which points to CI13TeV and which is
#                           set in setup.sh
#-----------------------------------------------------------------------------
import os, sys, re
#from histutil import *
from string import *
from array import array
from ROOT import *
#-----------------------------------------------------------------------------
LHAPATH = os.environ['LHAPDF_DATA_PATH']
CIPATH  = os.environ['CIPATH']

HISTNAMES = '''
_0.500_0.500
_0.500_1.000
_1.000_0.500
_1.000_1.000
_1.000_2.000
_2.000_1.000
_2.000_2.000
'''
HISTNAMES = split(strip(HISTNAMES))
MURMUF = [0, 3, 1, 4, 7, 5, 8]
#-----------------------------------------------------------------------------
def getBins():    
    filename = '%s/data/bins.txt' % CIPATH    
    records  = map(lambda x: map(atof, x),
                   map(split, open(filename).readlines()))
    pt = array('d')
    for ptlow, pthigh in records:
        pt.append(ptlow)
    ptlow, pthigh = records[-1]
    pt.append(pthigh)
    print "\t==> pT_min: %5.1f\tpT_max: %6.1f" % (pt[0], pt[-1])
    return pt

def makeHist(hname, name, pt, data):
    h = TH1D(hname, '', len(pt)-1, pt)
    h.GetXaxis().SetTitle('Jet p_{T} (GeV)')
    h.SetNdivisions(504, "Y")
    h.GetYaxis().SetTitle('d^{2}%s/dp_{T}dy (pb/GeV)' % name)
    h.SetLineWidth(1)
    
    for ii, d in enumerate(data):
        width = pt[ii+1]-pt[ii]
        d = d / (width*0.5)
        h.SetBinContent(ii+1, d)
        h.SetBinError(ii+1, 0)
    return h
#-----------------------------------------------------------------------------
def makeHistograms(dirname, name, number, pt, xsection):
    y = [0]*len(xsection)
    for ii in xrange(number):
        namen = '%s%d' % (name, ii)
        hfilename = '%s/%s.root' % (dirname, namen)
        hfile = TFile(hfilename, 'recreate')
        hist = []
        for which, prefix in [(1, 'nlo'),
                              (0, 'lo')]: # NLO, LO
            for jj, jjj in enumerate(MURMUF):
                hname = '%s%s' % (prefix, HISTNAMES[jj])
                for kk, xsect in enumerate(xsection):
                    cmd = 'xsect.%s(jjj, ii, which)' % name
                    y[kk] = eval(cmd)
                h = makeHist(hname, name, pt, y)
                hist.append(h)
        hfile.Write()
        hfile.Close()
#-----------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        print '''
    ./createCIhists.py PDFset [PDFindex=all]
        '''
        sys.exit(0)
        
    PDFset = argv[0]
    if   PDFset[:2] == 'CT':
        pdfsetdir = '../CT14'
    elif PDFset[:2] == 'MM':
        pdfsetdir = '../MMHT'
    elif PDFset[:2] == 'NN':
        pdfsetdir = '../NNPDF'
    else:
        print "wrong PDFset %s" % PDFset
        sys.exit(0)

    if argc > 1:
        PDFindexMin = atoi(argv[1])
        PDFindexMax = PDFindexMin
    else:
        PDFindexMin =   0
        PDFindexMax = 200

    print
    print "\t==> PDFset:            %s" % PDFset
    print "\t==> PDFset index(min): %d" % PDFindexMin
    print "\t==> PDFset index(max): %d" % PDFindexMax

    gSystem.Load("libCI.so")
    #from ROOT import CIXsection
    
    pt = getBins()
    
    hfile = []
    for index in xrange(PDFindexMin, PDFindexMax+1):
        dirname = '%s/%3.3d' % (pdfsetdir, index)        
        if index % 10 == 0:
            print dirname

        xsection = []
        for pT in pt[10:-1]:
            filename = '%s/data%4.4d.txt' % (dirname, pT)
            xsec = CIXsection(filename)
            if not xsec.good():
                print "*** problem"
                sys.exit()
            xsection.append(xsec) 

        makeHistograms(dirname, 'bi',  6, pt, xsection)
        makeHistograms(dirname, 'aig', 6, pt, xsection)
        makeHistograms(dirname, 'ai',  6, pt, xsection)
        
        makeHistograms(dirname, 'bij', 9, pt, xsection)
        makeHistograms(dirname, 'aijg',9, pt, xsection)
        makeHistograms(dirname, 'aij', 9, pt, xsection)
        
        makeHistograms(dirname, 'bi4', 4, pt, xsection)
        makeHistograms(dirname, 'ai4g',4, pt, xsection)                      
        makeHistograms(dirname, 'ai4', 4, pt, xsection)                      
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
