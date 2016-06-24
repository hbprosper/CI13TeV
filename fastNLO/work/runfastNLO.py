#!/usr/bin/env python
#-----------------------------------------------------------------------------
# Run fast NLO
# created 09-Oct-2012 Harrison B. Prosper
#         15-Apr-2015 HBP - replace call to fnlo-cppread by fnlo-tk-cppread
#         16-May-2016 HBP - adapt to version 2.3.1pre-2163 of fastNLO
#-----------------------------------------------------------------------------
import os, sys, re
from math import *
from histutil import *
from string import *
from glob import glob
from array import array
from time import sleep
from ROOT import *
#-----------------------------------------------------------------------------
TABLE   = 'InclusiveNJets_fnl5332g_v23_fix.tab'
if not os.environ.has_key('LHAPDF_DATA_PATH'):
    sys.exit("\n\t* you need to source setup.sh to define LHAPDF_DATA_PATH *\n")
LHAPATH = os.environ['LHAPDF_DATA_PATH']

getscales = re.compile('[0-9][.][0-9]+', re.M)

def main():
    print "\n\t<=== runfastNLO.py ===>"
    setStyle()

    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        records = map(lambda x: "\t%s\n" % strip(x),
                      os.popen("ls -1 $LHAPDF_DATA_PATH").readlines())
        record  = "%s%s%s" % (RED, joinfields(records), RESETCOLOR) 
        print '''
Usage:
    python runfastnlo.py PDFset [PDFindex=all]

    available PDFsets:
%s
        ''' % record
        sys.exit(0)

    table  = TABLE
    PDFset = argv[0]
    PDFsetfilename = '%s/%s' % (LHAPATH, PDFset)

    if not os.path.exists(PDFsetfilename):
        print "** can't find %s" % PDFset
        sys.exit(0)

    if   PDFset[:2] == 'CT':
        pdfsetdir = '../CT14'
    elif PDFset[:2] == 'MM':
        pdfsetdir = '../MMHT2014'
    elif PDFset[:2] == 'NN':
        pdfsetdir = '../NNPDF30'
    else:
        print "wrong PDFset %s" % PDFset
        sys.exit(0)

    if argc > 1:
        PDFindexMin = atoi(argv[1])
        PDFindexMax = PDFindexMin
    else:
        PDFindexMin =   0
        PDFindexMax = 200

    numberOfScales = 7        # Number of renormalization/factorization scales
    evolutionCode  = 'LHAPDF'
    
    print
    print "\t==> PDFset:            %s" % PDFset
    print "\t==> PDFset index(min): %d" % PDFindexMin
    print "\t==> PDFset index(max): %d" % PDFindexMax

    for index in xrange(PDFindexMin, PDFindexMax+1):
        dirname = '%s/%3.3d' % (pdfsetdir, index)
        cmd = 'mkdir -p %s' % dirname
        os.system(cmd)

        # run fastNLO program
        cmd = 'fnlo-tk-cppread %s %s %d %s %d' % (table,
                                                  PDFset,
                                                  numberOfScales,
                                                  evolutionCode,
                                                  index)
        print "\n%s" % cmd
        records = os.popen(cmd).readlines()

        # decode output
        ii = 0
        while ii < len(records):
            record = strip(records[ii])
            if not find(record, "My Cross Sections") > -1:
                ii += 1
                continue
            # next line should contain the scale factors
            ii += 1
            record = records[ii]
            if not find(record, "scale factor"):
                sys.exit(" ** expected scales in: \n%s\n" % record)

            t = getscales.findall(record)
            if len(t) != 2:
                sys.exit(" ** can't extract scales from: \n%s\n" % record)
            mur, muf = t

            # open output file
            outfilename = '%s/qcd_%s_%s.txt' % (dirname, mur, muf)
            print "\t%s" % outfilename
            out = open(outfilename, 'w')
            rec = '%7s %7s %20s %20s' % \
              ('pT-low', 'pT-high',
               'LO-xsection', 'NLO-xsection')
            out.write('%s\n' % rec)
            print rec
            # skip next three lines
            ii += 3

            # extract cross sections for |y| < 0.5
            all_is_well = True
            nlines = 0
            while all_is_well:
                the_sun_shines = nlines < 100
                ii += 1
                nlines += 1
                record = records[ii]
                t = split(record)
                y = t[4]

                # break if we have run out of data for |y| < 0.5
                all_is_well = (nlines < 100) and (y[:3] == '0.5')
                if not all_is_well: break

                pTlow = t[6]
                pTupp = t[7]
                xLO   = t[9]
                xNLO  = t[10]
                t = (pTlow, pTupp, xLO, xNLO)
                rec = '%7s %7s %20s %20s' % t
                out.write('%s\n' % rec)
                print rec
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print 'ciao!'
