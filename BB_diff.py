import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import indexList, GetAntName, GetPSpecScan, SPW2DATA_DESC_ID
from casatools import quanta as qatool
import analysisUtils as au
import datetime
qa = qatool()
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antList', metavar='antList',
    help='Antenna List (refant at first)', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List e.g. 2, 5, 10', default='')
parser.add_option('-s', dest='spwID', metavar='spwID',
    help='SPW ID   e.g. 13', default='13')
parser.add_option('-M', dest='plotRange', metavar='plotRange',
    help='Max y-axis range in [ps] e.g. 10.0', default='10.0')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
plotAntList = [ant for ant in options.antList.split(',')]
spwID  =  int(options.spwID)
scanList = [int(scan) for scan in options.scanList.split(',')]
plotMax = float(options.plotRange)
'''
prefix = 'uid___A002_X12e95ca_X7383'
plotAntList = ['DA45','DA41','DA48','DA50','DA64','DV09','DV18','DV21','DV22']
scanList = [3,6,8,10]
spwID = 13
plotMax = 0.025
'''
scanNum  = len(scanList)
#-------- Main process
msfile = prefix + '.ms'
msmd.open(msfile)
antList = GetAntName(msfile)
pp = PdfPages('BBD_%s_SPW%d.pdf' % (prefix, spwID))
figBB = plt.figure(0, figsize=(8,11))
figBB.suptitle('%s SQLD Power (w.r.t. %s) SPW=%d' % (prefix, plotAntList[0], spwID))
#figBB.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/BB_scan.py to generate this plot', fontsize=4)
scanNum = len(scanList)
plotAntNum = len(plotAntList) - 1
rowNum = int(min(4, np.floor(np.sqrt(plotAntNum))))
#-------- Reference Power
ant_index = np.where(antList == plotAntList[0])[0][0]
TS, BBX, BBY = [], [], []
for scan_index, scanID in enumerate(scanList):
    timeStamp, RefPower = GetPSpecScan(msfile, ant_index, spwID, scanID)
    TS = TS + timeStamp.tolist()
    #BB = BB + np.mean(RefPower, axis=(0,1)).tolist()
    #TS = TS + [np.mean(timeStamp)]
    BBX = BBX + RefPower[0,0].tolist()
    BBY = BBY + RefPower[1,0].tolist()
DT = [datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f') for mjdSec in TS]
BBXref, BBYref = np.array(BBX), np.array(BBY)
antMap = indexList(plotAntList, antList)[1:]
#-------- Loop for each antenna
for row_index, ant_index in enumerate(antMap):
    BBPL = figBB.add_subplot( int(np.ceil(plotAntNum/rowNum)), rowNum, row_index + 1 )
    BBPL.yaxis.offsetText.set_fontsize(3)
    BBPL.tick_params(labelsize=4); BBPL.tick_params(axis='x')
    BBPL.axis([np.min(DT), np.max(DT), -plotMax, plotMax])
    BBPL.grid(axis='y', linestyle='--', color='gray')
    BBX, BBY = [], []
    for scan_index, scanID in enumerate(scanList):
        timeStamp, BBPower = GetPSpecScan(msfile, ant_index, spwID, scanID)
        #BB = BB + np.mean(BBPower, axis=(0,1)).tolist()
        BBX = BBX + BBPower[0,0].tolist()
        BBY = BBY + BBPower[1,0].tolist()
    BBpowerX = np.array(BBX) / BBXref
    BBpowerX = BBpowerX - BBpowerX[0]
    BBpowerY = np.array(BBY) / BBYref
    BBpowerY = BBpowerY - BBpowerY[0]
    text_sd = antList[ant_index]
    BBPL.plot( DT, BBpowerX, '-', linewidth=0.25, label='pol-X'); BBPL.plot( DT, BBpowerX, '.' )
    BBPL.plot( DT, BBpowerY, '-', linewidth=0.25, label='pol-Y'); BBPL.plot( DT, BBpowerY, '.' )
    BBPL.text(DT[0], 0.8*plotMax, text_sd, fontsize=5)
    if row_index == 0: BBPL.legend()
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()
