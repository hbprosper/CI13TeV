#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        computeLimits.py
# Description: read CI workspace and compute limits
# Created:     02-Mar-2014 HBP
#              22-May-2014 HBP slighty modified version of program
#                          used to compute 7 TeV limits
#              09-Jun-2014 HBP add check function
#              15-Nov-2014 HBP use new pdf
#              18-Sep-2016 HBP
#              18-May-2018 HBP permit model-dependent lambda range
#              22-May-2018 HBP attempt to find range automatically
#              26-May-2018 HBP add a log-file
#-----------------------------------------------------------------
import os,sys,re, optparse
from array import array
from time import sleep, ctime
from histutil import *
from string import *
from ROOT import *
#-----------------------------------------------------------------
EXPECTED = True    # if true compute expected limits
#-----------------------------------------------------------------
ENERGY   = 13      # TeV
LUMI     = 35100.0 # 1/fb
BINMIN   =  1      # first bin to use (ROOT bin number convention)
BINMAX   = 28      # last bin to use
LMIN     =  0.0    # lower limit of lambda = 1/Lambda^2
LMAX     = {'LL': {1: 0.0015, -1: 0.0040},
            'RR': {1: 0.0015, -1: 0.0040},
            'VV': {1: 0.0015, -1: 0.0040},
            'AA': {1: 0.0015, -1: 0.0040},
            'V-A':{1: 0.0060, -1: 0.0060}} # upper limit on lambda

YMAX     = 0.060   # maximum Y limit of posterior density plot
WSPACE   = 'CI'    # name of workspace

RELTOL   = 1.e-4   # minimum ratio of min(like)/max(like)
# Models
#                        kappa
#                   0   1   2   3   4   5
KAPPA   = {'LL' : [-1,  0,  0,  0,  0,  0],
           'RR' : [ 0,  0,  0,  0, -1,  0],
           'VV' : [-1,  0, -2,  0, -1,  0],
           'AA' : [-1,  0, +2,  0, -1,  0],
           'V-A': [ 0,  0, -2,  0,  0,  0]}
MODEL   = ['LL', 'RR', 'VV', 'AA', 'V-A']

CL68    = 0.683
P16     = (1-CL68)/2
P50     = 0.5
P84     = (1+CL68)/2
P95     = 0.95
PERCENTILE=[P16, P50, P84, P95]
rand3   = TRandom3()
#-----------------------------------------------------------------
def decodeCommandLine():
    VERSION = '26-May-2018'
    USAGE = '''
    python computeLimits.py [options] <workspace-file> [LL RR etc.]

    options 
       -O              compute observed limit [default compute expected limit]

    options for expected limit
       -L<Lambda>      Lambda value for expected limit [default no CI]
       -l<lumi>        integrated luminosity           [default 35.1/fb]
       -n<ntrials>     number of trials for computing scatter of percentiles

    '''
    parser = optparse.OptionParser(usage=USAGE,
                                   version=VERSION)

    parser.add_option('-L', '--Lambda',
                      action='store',
                      dest='Lambda',
                      type='float',
                      default=-1.0,
                      help='CI mass scale')

    parser.add_option('-l', '--lumi',
                      action='store',
                      dest='lumi',
                      type='float',
                      default=LUMI/1000,
                      help='integrated luminosity in 1/fb')

    parser.add_option('-O', '--observed',
                      action="store_true",
                      dest='observed',
                      default=False,
                      help='compute observed limits')

    parser.add_option('-n', '--ntrials',
                      action="store",
                      type='int',
                      dest='ntrials',
                      default=-1,
                      help='number of trials for computing expected limits')       
        
    options, args = parser.parse_args()
    if len(args) == 0: sys.exit(USAGE)
    filename = args[0]
    models   = args[1:]
    if models == []: models = MODEL
    
    return (filename, not options.observed,
                options.Lambda, 1e3*options.lumi, models,
                options.ntrials)
#-----------------------------------------------------------------
#-----------------------------------------------------------------
def optimizeRange(model, poi, nstep=20):
    print "optimize range of lambda"
    poi_lower = poi.getMin()
    poi_upper = poi.getMax()
    bayes     = Bayes(model, poi_lower, poi_upper)

    swatch = TStopwatch()
    swatch.Start()
    
    step = (poi_upper - poi_lower)/nstep

    big_p =-1;
    small_p=1e300
    for ii in range(nstep+1):
        x = ii * step
        p = model(x)
        if p > big_p:
            big_p = p
            
        if big_p > p:
            if p < small_p:
                small_p = p

        if small_p / big_p < 1.e-3:
            poi_upper = x
            print "\t%8.3f\t%8.2e" % (x, p)
            break
        
        print "\t%8.3f\t%8.2e" % (x, p)
        
    poi.setRange(poi_lower, poi_upper)
    print "\t\t\t==> real time: %8.3f s\n" % swatch.RealTime()            
    print "optimized range: %8.4f ... %8.4f 1/TeV^2\n" %\
      (poi.getMin(), poi.getMax())            
#-----------------------------------------------------------------
def computePercentiles(log, pdf, poi):
    bayes  = Bayes(pdf, poi.getMin(), poi.getMax())
    limits = []
    record = '%10s %10s' % ('percentile', 'limit')
    print record
    log.write('%s\n' % record)
    
    for p in PERCENTILE:
        lb = bayes.percentile(p)            
        if lb > 0:
            limits.append(1.0/sqrt(lb))
        else:
            limits.append(-1)

        record = "%10.3f %10.1f" % (p, limits[-1])
        print record
        log.write('%s\n' % record)            
    return (limits, bayes)

def set_kappa(ws, kappa):
    for ii in xrange(len(kappa)):
        vname = 'kappa%d' % ii
        ws.var(vname).setVal(kappa[ii])
        
def set_lambda_kappa(ws, lam, kappa):
    ws.var('lambda').setVal(lam)
    for ii in xrange(len(kappa)):
        vname = 'kappa%d' % ii
        ws.var(vname).setVal(kappa[ii])            
#-----------------------------------------------------------------
def ExpectedLimits(log, tm, bayes, model, luminosity, Lambda, ntrials=200):

    swatch = TStopwatch()
    swatch.Start()

    tm += 0.05
    i   = int(tm / 0.01)
    tm  = ntrials * (i+1)*0.01 / 60
    
    print "\t=> computing coverage. this will take about %5.2f minutes\n" % tm
    
    lam = 1.0/Lambda**2
    np  = len(PERCENTILE)
    timeleft = TimeLeft(ntrials)

    c   = [0.0]*np
    x1  = [0.0]*np
    x2  = [0.0]*np
    N   = 0
    for trial in range(ntrials):
        # use Asimov dataset and add fluctuations
        model.setAsimov(True, True, luminosity, lam)
        model.initialize(-1) # average over syst. uncert.
        bayes.reset()
        
        limits = [0]*np
        for i, p in enumerate(PERCENTILE):
            limits[i] = bayes.percentile(p)
        if min(limits) <= 0: continue
        N += 1
        
        limits = [ 1.0/sqrt(x) for x in limits ]           
        for ii in range(np):
            if Lambda > limits[ii]: c[ii] += 1.0
            x1[ii] += limits[ii]
            x2[ii] += limits[ii]**2

        if (N-1) % 10 == 0:
            t = [trial] + [ x / N for x in c ]
            t.append(timeleft(trial))
            print "trial: %5d\t%5.2f, %5.2f, %5.2f, %5.2f\t%s" % tuple(t) 
    print

    record = '%10s %16s\t%16s' % ('percentile', 'expected limit', 'coverage')
    log.write('%s\n' % record)    
    
    for ii in range(np):
        x1[ii] /= N
        x2[ii] /= N
        x = x2[ii] - x1[ii]**2
        if x < 0: x = 0.0
        x2[ii]  = sqrt(x2[ii]-x1[ii]**2)
        
        limit = x1[ii]
        dlimit= x2[ii]
        p = PERCENTILE[ii]
        t  = c[ii] / N
        dt = sqrt(t*(1-t)/N)        
        record = "%10.3f %8.1f (%5.1f)\t%8.3f (%5.3f)" % (p, limit, dlimit, t, dt)
        print record
        log.write('%s\n' % record)
        
    print "\t\t\t==> real time: %8.3f s\n" % swatch.RealTime()                    
#-----------------------------------------------------------------
def check(o, message):
    if o == None:
        sys.exit("** %s **" % message)
#-----------------------------------------------------------------        
def makePlot(ws, bayes, poi, hname,
             color=kBlack,
             lstyle=1,
             lwidth=2,
             ymin=0.0,
             ymax=YMAX):

    model  = ws.pdf('model')
 
    xmin   = poi.getMin()
    xmax   = poi.getMax()
    xbins  = 100
    xstep  = (xmax-xmin)/xbins

    h = TH1F(hname, "", xbins, xmin, xmax)

    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.GetXaxis().SetTitle('#lambda (TeV^{-2})')
    h.GetXaxis().SetNdivisions(504)

    h.GetYaxis().SetTitle('p(#lambda | D)')
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetTitleOffset(1.7)
    h.SetLineWidth(lwidth)

    max_f = -1
    max_x = 0.0
    for ii in xrange(xbins):
        x = xmin + (ii+0.5)*xstep
        y = bayes.posterior(x) * xstep
        z = exp(model.logProfileLikelihood(x))*xstep
        f = model.underflowFraction()
        if f > max_f:
            max_f = f
            max_x = x
        h.SetBinContent(ii+1, y)
        h.SetBinError(ii+1, 0)

    total = h.Integral()
    if total == 0:
        sys.exit("total is ZERO!")
    else:
        "Integral: %10.3f" % total

    print '\t=> maximum underflow fraction %10.3f at lambda = %10.5f\n' % \
      (max_f, max_x)
      
    h.Scale(1.0/h.Integral())
    ymax = h.GetMaximum()
    ii   = int(ymax/0.04)
    ymax = (ii+1)*0.04
    h.SetMinimum(0.0)
    h.SetMaximum(ymax)

    return h
#-----------------------------------------------------------------
def setValues(ws, d, scale=1.0):
    if type(d) == type(""): d = ws.set(d)
    iterator = d.iterator()
    ii = 0
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        ws.var(name).setVal(d[name].getVal()/scale)
        print "%5d\t%10.2f" % (ii, ws.var(name).getVal())
        ii += 1
#-----------------------------------------------------------------
def printValues(ws, d):
    if type(d) == type(""):
        print "%5d %10s\t%10.4f" % (1, d, ws.var(d).getVal())
        return

    iterator = d.iterator()
    ii = 0
    while 1:
        ii += 1
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        print "%5d %10s\t%10.2f" % (ii, name, ws.var(name).getVal())
#-----------------------------------------------------------------
def toVector(ws, setname):
    data = ws.set(setname)
    if data == None: return None
    iterator = data.iterator()
    vdata = vector('double')()
    while 1:
        o = iterator.Next()
        if  o == None: break
        name = o.GetName()
        vdata.push_back( ws.var(name).getVal() )
    return vdata
#-----------------------------------------------------------------
def main():

    filename, expected, Lambda, luminosity, models, ntrials = decodeCommandLine()
    
    # --------------------------------------
    # load various codes needed for the
    # calculatioms
    # --------------------------------------
    gSystem.Load('../CI/lib/libCI')

    # set up some standard graphics style
    setStyle()
    
    # --------------------------------------
    # load workspace into memory
    # --------------------------------------
    wfile = TFile(filename)
    if not wfile.IsOpen():
        sys.exit("** can't open %s" % filename)
        
    ws = wfile.Get(WSPACE)
    check(ws, "can't access workspace %s" % WSPACE)
    
    energy ='%d' % ENERGY
    lumin  = luminosity / 1000
    lumi   = '%5.2f' % lumin
        
    prefix = nameonly(filename)
    dname  = split(prefix, '_')[0]
    os.system('mkdir -p figures/%s' % dname)
    prefix = replace(prefix, '_workspace', '')
    
    # --------------------------------------
    # get model etc.
    # --------------------------------------
    model = ws.pdf('model')
    Nset  = ws.set('Nset')
    poi   = ws.var('lambda')        # parameter of interest (lambda)
    
    data  = toVector(ws, 'Nset')
    nbins = model.numberOfBins()
    lam   = None
    kappa = None
    
    # if expected = True then create and use an Asimov data set
    if expected:

        if Lambda <= 0:
            model.setAsimov(True)
            SIGN = [-1, 1]
            
        else:
            if models[0][0] == '~':
                models[0] = models[0][1:]
                sign =-1
            else:
                sign = 1

            lam  = 1.0/Lambda**2                
            key  = models[0]
            SIGN = [sign]
            if KAPPA.has_key(key):
                kappa = [ sign * x for x in KAPPA[key] ]
                set_lambda_kappa(ws, lam, kappa)
                
            printValues(ws, ws.set('kappaset'))
            model.setAsimov(True, False, luminosity, lam)
            
        print "\n\t== use Asimov data set (Lumi = %10.1f / pb)" % luminosity
        
        Asimov = model.Asimov()
        nn = 0
        for ii in xrange(Asimov.size()):
            print "%4d\t%10.1f" % (ii+1, Asimov[ii]),
            nn += 1
            if nn > 3:
                print
                nn = 0
        print        
    else:
        model.setAsimov(False)
        SIGN = [-1, 1]

    # --------------------------------------
    # set range of bins to use.
    # NOTE: use ROOT bin labeling convention
    # first bin is 1 and last is nbins
    # --------------------------------------
    binmin = BINMIN
    binmax = BINMAX
    model.setBinRange(binmin, binmax)
    
    print "="*80
    print "input filename:   %s" % filename
    print "workspace:        %s" % WSPACE
    print "models:           %s" % models
    print "bin range:        [%d ... %d]" % (binmin, binmax)
    print "integrated lumi:  %9.2f / pb" % luminosity
    if expected:
        print 
        print "compute expected limits with"
        print "Lambda:       %9.2f TeV" % Lambda
        
    print "="*80

    postfix = '_l%3.3d' % int(lumin+0.5)
    if Lambda > 0:
        postfix += '_L%3.3d' % int(Lambda)
                
    if expected:
        postfix += '_expected'
    else:
        postfix += '_observed'
        
    # --------------------------------------
    # compute limit for each model
    # --------------------------------------
    # create wrapper for RooFit model
    wrapped_model = PDFWrapper(model, Nset, poi)
    
    for key in models:
        
        for sign in SIGN:

            # --------------------------------------            
            # set kappa values for current model
            # --------------------------------------            
            kappa = [ sign * x for x in KAPPA[key] ]
            set_kappa(ws, kappa)
            printValues(ws, ws.set('kappaset'))

            # --------------------------------------
            # find a reasonable range for lambda
            # --------------------------------------
            poi.setRange(0.0, 0.02)            
            optimizeRange(wrapped_model, poi)

            # --------------------------------------
            # create canvas
            # --------------------------------------            
            if sign > 0:
                print "model(%s, constructive)" % key
                name = '%s_constructive' % key
                fname = 'figures/%s/%s_limit_%s%s' % \
                  (dname, prefix, name, postfix)
                clike1 = TCanvas(fname, fname, 10, 10, 500, 500)
                clike = clike1

            else:
                print "model(%s, destructive)" % key
                name = '%s_destructive' % key
                fname = 'figures/%s/%s_limit_%s%s' % \
                  (dname, prefix, name, postfix)
                clike2 = TCanvas(fname, fname, 515, 10, 500, 500)
                clike = clike2
                
            # --------------------------------------
            # open log file
            # --------------------------------------                            
            logfile = '%s.log' % fname
            print 'logfile: %s' % logfile
            log = open(logfile, 'w')
            log.write('''created: %(time)s\tsqrt(s): %(energy)sTeV\tl: %(lumi)s/fb
model:   %(model)s\tkappa: %(kappa)s
            ''' % {'time':   ctime(),
                   'energy': energy,
                   'lumi':   lumi,
                   'model':  key,
                   'kappa':  kappa})            

            if Lambda > 0:
                log.write('computed with signal mass scale Lambda: %8.1fTeV\n' % Lambda)
                
           # delete histograms to avoid memory leaks
            try:
                del hnom
            except:
                pass
            try:
                del havg
            except:
                pass
            try:
                del hclone
            except:
                pass
            
            # --------------------------------------
            # compute limits with systematic
            # uncertainties
            # --------------------------------------
            record = '\n\tWITH systematic uncertainties'
            print record
            log.write('%s\n' % record)
            
            swatch = TStopwatch()
            swatch.Start()

            model.initialize(-1) # average over syst. uncert.
            limits_avg, bayes_avg = computePercentiles(log, wrapped_model, poi)

            tm = swatch.RealTime()
            print "\t\t=> real time: %10.2f seconds" % tm
            
            havg = makePlot(ws, bayes_avg, poi, "havg",
                                color=kRed, lstyle=1, lwidth=2)            

            # --------------------------------------
            # estimate fluctuations in percentiles
            # --------------------------------------            
            if ntrials > 0:
                set_kappa(ws, kappa)
                ExpectedLimits(log, tm, bayes_avg, model,
                               luminosity, Lambda, ntrials)
            
            # --------------------------------------
            # compute nominal limits
            # --------------------------------------
            record = '\n\tNO systematic uncertainties'
            print record
            log.write('%s\n' % record)

            swatch = TStopwatch()
            swatch.Start()
            
            model.initialize(0) # use nominal cross section            
            limits_nom, bayes_nom = computePercentiles(log, wrapped_model, poi)
            
            tm = swatch.RealTime()
            print "\t\t=> real time: %10.2f seconds" % tm
            
            log.close()
            
            hnom = makePlot(ws, bayes_nom, poi, "hnom",
                                color=kBlue, lstyle=2)
            
            # --------------------------------------
            # update plot
            # --------------------------------------
            clike.cd()
            hnom.Draw('l')
            havg.Draw('l same')
            clike.Update()
            gSystem.ProcessEvents()
            
            # shade 95% region
            hclone = havg.Clone()
            hclone.SetAxisRange(poi.getMin(), limits_avg[-1])
            hclone.SetFillStyle(3001)
            hclone.SetFillColor(30)
            clike.cd()
            hclone.Draw('l same')
            gSystem.ProcessEvents()
            
            # --------------------------------------
            # plot posterior density and limits
            # --------------------------------------
            clike.cd()
            scribe = addTitle('CMS Preliminary  '\
                              '#surds=%sTeV  L=%sfb^{-1}' % \
                              (energy, lumi),
                              0.04)
            scribe.vspace()
            scribe.write("model: %s, #kappa = %s" % (key, kappa), 0.06)
            scribe.write("#color[4]{#Lambda} > %4.1fTeV @ 95%s CL" % \
                         (limits_nom[-1], '%'), 0.06)
            scribe.write("#color[2]{#Lambda} > %4.1fTeV @ 95%s CL" % \
                            (limits_avg[-1], '%'), 0.06)
                            
            for i, limit in enumerate(limits_avg):
                scribe.write("%8.4f\t %4.1f" % (PERCENTILE[i], limit), 0.39)
            clike.Update()
            clike.SaveAs('.png')
            gSystem.ProcessEvents()
            
        sleep(5)
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
