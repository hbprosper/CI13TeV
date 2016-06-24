#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        checkWorkspace.py
# Description: check portability of workspace
# Created:     07-Jun-2015 HBP
#-----------------------------------------------------------------
import os,sys
from ROOT import TFile
#-----------------------------------------------------------------
def check(o, message):
    if o == None:
        print "** %s **" % message
        sys.exit(0)
#-----------------------------------------------------------------
WSPACE = 'CI'
def main():
    # --------------------------------------
    # load workspace into memory
    # --------------------------------------
    argv = sys.argv[1:]
    if len(argv) == 0:
        print '''
    Usage:
         ./checkWorkspace.py <workspace-file>
        '''
        sys.exit(0)
        
    filename = argv[0]
    wfile = TFile(filename)
    if not wfile.IsOpen():
        print "** can't open %s" % filename
        sys.exit(0)
        
    ws = wfile.Get(WSPACE)
    check(ws, "can't access workspace %s" % WSPACE)
    
    model = ws.pdf('model')
    check(model, "can't access model")
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
