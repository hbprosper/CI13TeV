#!/usr/bin/env python
import os, sys
from string import *

LUMI = 0.07152 # 1/fb

DATA = '''
527.5                          20.5                        10400 +/- 275.801(stat)
570                          22                        6412 +/- 175.186(stat)
615                          23                        4601 +/- 135.547(stat)
662                          24                        3019 +/- 70.5762(stat)
711.5                          25.5                        2083 +/- 62.1048(stat)
763.5                          26.5                        1353 +/- 36.7831(stat)
818                          28                        906 +/- 30.0998(stat)
875.5                          29.5                        600 +/- 24.4949(stat)
936                          31                        399 +/- 19.975(stat)
999.5                          32.5                        294 +/- 17.1464(stat)
1066.5                          34.5                        182 +/- 13.4907(stat)
1136.5                          35.5                        123 +/- 11.0905(stat)
1210                          38                        100 +/- 10(stat)
1287.5                          39.5                        48 +/- 6.9282(stat)
1368.5                          41.5                        32 +/- 5.65685(stat)
1453.5                          43.5                        18 +/- 4.24264(stat)
1542.5                          45.5                        18 +/- 4.24264(stat)
1636                          48                        12 +/- 3.4641(stat)
1734                          50                        10 +/- 3.16228(stat)
1837                          53                        4 +/- 2(stat)
1945                          55                        1 +/- 1(stat)
'''
DATA = map(lambda x: map(atof, (x[0], x[1], x[2])),
           map(split, split(strip(DATA),'\n')))

OLDDATA = '''
Bin = 20, Bin Range (GeV) = 507 - 548, Yield = 722864
Bin = 21, Bin Range (GeV) = 548 - 592, Yield = 468062
Bin = 22, Bin Range (GeV) = 592 - 638, Yield = 296032
Bin = 23, Bin Range (GeV) = 638 - 686, Yield = 186497
Bin = 24, Bin Range (GeV) = 686 - 737, Yield = 120580
Bin = 25, Bin Range (GeV) = 737 - 790, Yield = 76129
Bin = 26, Bin Range (GeV) = 790 - 846, Yield = 48454
Bin = 27, Bin Range (GeV) = 846 - 905, Yield = 31121
Bin = 28, Bin Range (GeV) = 905 - 967, Yield = 19639
Bin = 29, Bin Range (GeV) = 967 - 1032, Yield = 12373
Bin = 30, Bin Range (GeV) = 1032 - 1101, Yield = 7746
Bin = 31, Bin Range (GeV) = 1101 - 1172, Yield = 4670
Bin = 32, Bin Range (GeV) = 1172 - 1248, Yield = 2930
Bin = 33, Bin Range (GeV) = 1248 - 1327, Yield = 1746
Bin = 34, Bin Range (GeV) = 1327 - 1410, Yield = 1111
Bin = 35, Bin Range (GeV) = 1410 - 1497, Yield = 602
Bin = 36, Bin Range (GeV) = 1497 - 1588, Yield = 367
Bin = 37, Bin Range (GeV) = 1588 - 1784, Yield = 311
Bin = 38, Bin Range (GeV) = 1784 - 1890, Yield = 111
Bin = 39, Bin Range (GeV) = 1890 - 2000, Yield = 6
Bin = 40, Bin Range (GeV) = 2000 - 2116, Yield = 0
Bin = 41, Bin Range (GeV) = 2116 - 2238, Yield = 0

'''
OLDDATA = map(lambda x: map(atof, (x[2][:-1], x[7], x[9][:-1])),
           map(split, split(strip(OLDDATA),'\n')))

def main():
    global DATA, OLDDATA
    
    DATA = DATA[2:]
    OLDDATA = OLDDATA[2:]
    
    M = len(DATA)
    out = open('data_13TeV_L000.07152ifb.txt', 'w')
    record = ' %9s %9s %9s' % ('lower', 'upper', 'count')    
    out.write('%s\n' % record)
    
    for ii in xrange(M):
        mid, h, c  = DATA[ii]

        b, l1, u1 = OLDDATA[ii]
        low = mid - h
        upp = mid + h
            
        record = ' %9.0f %9.0f %9.0f' % (low, upp, c)
        out.write('%s\n' % record)
        print record, l1, u1 


main()

