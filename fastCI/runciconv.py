#!/usr/bin/env python
#------------------------------------------------------------------------------
# run ciconv
#------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
from time import sleep
#------------------------------------------------------------------------------
PDFSETS = '''
MMHT   MMHT2014nlo68cl_rand1234
NNPDF  NNPDF30_nlo_as_0118_1000
'''
PDFSETS = '''
CT14 CT14nlo_rand1234
'''
PDFSETS = map(split, split(strip(PDFSETS), '\n'))
NPDFS   = 200
BASE    = './'
getnumber = re.compile('[0-9]+')
#------------------------------------------------------------------------------
def main():
    gridfilenames = map(lambda x: split(x, 'fgrid/')[-1],
                        glob("%sfgrid/*.sum" % BASE))
    gridfilenames.sort()
    
    for pdf, pdfset in PDFSETS:
        for ii in xrange(NPDFS+1):
            member = ii
            cmd = 'mkdir -p %s%s/%3.3d' % (BASE, pdf, ii)
            os.system(cmd)
            for gridfilename in gridfilenames:
                number = getnumber.findall(gridfilename)[0]
                outfilename = "%s%s/%3.3d/data%s.txt" % \
                  (BASE, pdf, ii, number)
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
