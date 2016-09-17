#!/usr/bin/env python
#------------------------------------------------------------------------------
# run ciconv
# created: 16-Sep-2016 HBP - a version suitable to be run in a batch job
#------------------------------------------------------------------------------
import os, sys, re
from string import *
from glob import glob
from ROOT import *
from time import sleep
#------------------------------------------------------------------------------
BASE    = './'
getnumber = re.compile('[0-9]+')
#------------------------------------------------------------------------------
def main():
    gridfilenames = map(lambda x: split(x, 'fgrid/')[-1],
                        glob("%sfgrid/*.sum" % BASE))
    gridfilenames.sort()
    pdf = 'CT14'
    pdfset = 'CT14nlo_rand1234'

    if len(sys.argv[1:]) > 0:
	member = atoi(sys.argv[1])
    else:
	member = 0

    cmd = 'mkdir -p %s%s/%3.3d' % (BASE, pdf, member)
    os.system(cmd)
    for gridfilename in gridfilenames:
        number = getnumber.findall(gridfilename)[0]
        outfilename = "%s%s/%3.3d/data%s.txt" % (BASE, pdf, member, number)
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
