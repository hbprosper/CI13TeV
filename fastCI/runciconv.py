#!/usr/bin/env python
#------------------------------------------------------------------------------
# run ciconv
#------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
#------------------------------------------------------------------------------
PDFSETS = '''
MMHT   MMHT2014nlo68cl_rand1234
NNPDF  NNPDF30_nlo_as_0118_1000
'''
PDFSETS = map(split, split(strip(PDFSETS), '\n'))
#------------------------------------------------------------------------------
def main():
    gridfilenames = glob("fgrid/data*.sum")
    gridfilenames = map(lambda x: split(x,'/')[-1], gridfilenames)
    gridfilenames.sort()
    
    for pdf, pdfset in PDFSETS:
        for ii in xrange(201):
            member = ii
            cmd = 'mkdir -p %s/%3.3d' % (pdf, ii)
            os.system(cmd)
            for gridfilename in gridfilenames:
                a, b = split(gridfilename, '_')
                outfilename = "%s/%3.3d/%s.txt" % \
                  (pdf, ii, split(a, '/data')[-1])
                cmd = 'ciconv %s %d %s %s' % (pdfset,
                                              member,
                                              gridfilename,
                                              outfilename)
                print
                print cmd
                os.system(cmd)
#------------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "choi!"
