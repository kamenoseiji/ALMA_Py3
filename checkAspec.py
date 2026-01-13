import numpy as np
from scipy import stats
from optparse import OptionParser
from interferometry import GetAntName, GetTimerecord, GetChNum, GetPSpecScan, indexList
from Plotters import plotAC
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antName', metavar='antName',
    help='Antennas to plot e.g. PM01,PM03', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List  e.g. 3,6,8', default='3')
#
(options, args) = parser.parse_args()
#-------- checkAspec : power spectrum (autocorrelation)
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
scanList = [int(scan) for scan in options.scanList.split(',')]
spwList = [int(spw) for spw in options.spwList.split(',')]
'''
prefix = 'uid___A002_X133b089_X18b6f'
scanList = [4]
spwList = [19]
'''
#-------- Procedures
msfile = prefix + '.ms'
antList = GetAntName(msfile)
#plotAntList = antList if options.antName == '' else [ant for ant in options.antName.split(',')]
plotAntList = antList
antNum = len(plotAntList)
scanNum = len(scanList)
spwNum = len(spwList)
#-------- Time Records
timeList = []
for scan_index, scan in enumerate(scanList):
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], scan)
    print('scan%d : %d time records' % (scan, len(timeStamp)))
    timeList = timeList + timeStamp.tolist()
#
timeNum = len(timeList)
#-------- SPW records
AC = []
chNumList, freqList = [], []
for spw_index, spw in enumerate(spwList):
    chNum, chWid, Freq = GetChNum(msfile, spw)
    chNumList = chNumList + [chNum]
    freqList = freqList + [Freq* 1.0e-9]
    AC = AC + [np.zeros([antNum, timeNum, 2, chNum])]
#
#-------- Load autocorrelation power spectra
polIndex = [0, -1]
for ant_index, antID in enumerate(indexList(np.array(plotAntList), antList)):
    for spw_index, spw in enumerate(spwList):
        timePointer = 0
        for scan_index, scan in enumerate(scanList):
            timeStamp, Pspec = GetPSpecScan(msfile, antID, spw, scan)    # Pspec[pol, ch, time]
            recNum = len(timeStamp)
            AC[spw_index][ant_index, timePointer:(timePointer + recNum)] = Pspec[polIndex].transpose(2, 0, 1) # AC[spw][ant, time, pol, ch] 
            timePointer += recNum
        #
    #
#
#-------- Plot BP
plotAC(prefix + '-Scan%d' % (scanList[0]), plotAntList, spwList, freqList, AC)
