import numpy as np
from Plotters import plotTsysDic
import analysisUtils as au
from interferometry import GetAntName, GetChNum, GetBandNames
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-t', dest='PLOTTSYS', metavar='PLOTTSYS',
    help='Plot Tsys and Trx spectra', action="store_true")
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
PLOTTSYS= options.PLOTTSYS
'''
prefix = 'uid___A002_X12fb842_X40a0'
PLOTTSYS = True
'''
#-------- Read SYSCAL table
msfile = prefix + '.ms'
tb.open(msfile + '/SYSCAL')
colnameList = tb.colnames()
timeStamp = tb.getcol('TIME')
timeSyscalList = timeStamp[ np.where(np.diff(timeStamp) != 0.0)[0].tolist() + [-1] ]
antID = tb.getcol('ANTENNA_ID')
spwID = tb.getcol('SPECTRAL_WINDOW_ID')
TRX = tb.getcol('TRX_SPECTRUM')
TSY = tb.getcol('TSYS_SPECTRUM')
tb.close()
antList, spwList = GetAntName(prefix + '.ms'), np.unique(spwID)
spwDic = dict(zip(spwList, [[]]*len(spwList)))  # SPW
antNum, spwNum = len(antList), len(spwList)
logFile  = open(prefix + '-TelCal.log', 'w')
msmd.open(msfile)
scanList = msmd.scansforintent("CALIBRATE_ATMOSPHERE*")
timeAMB  = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
timeList = [np.median(msmd.timesforscan(scan)) for scan in scanList]
msmd.close()
BandNames = GetBandNames(msfile, spwList); UniqBands = list(set(BandNames))
NumBands = len(UniqBands)
TsysDic = dict(zip(scanList, [[]]*len(scanList)))  # Scan-based Tsys and Trx
#-------- Read SPW frequency
for spw_index, spw in enumerate(spwList):
    chNum, chWid, freq = GetChNum(msfile, spw)
    spwDic[spw] = {'band': BandNames[spw_index], 'chNum': chNum, 'chWid': chWid, 'freq' : freq}
for scan_index, scan in enumerate(scanList):
    data_index = np.where(timeStamp == timeSyscalList[scan_index])[0].tolist()
    spwsInScan = np.unique(spwID[data_index]).tolist()
    antsInScan = np.unique(antID[data_index]).tolist()
    timeLabel = 'Scan %d : %s' % (scan, au.call_qa_time('%fs' % (timeList[scan_index]), form='fits'))
    logFile.write(timeLabel + '\n'); print(timeLabel)
    TsysDic[scan] = {
        'antList': antList[antsInScan],
        'mjdSec' : timeList[scan_index],
        'spwList': spwsInScan,
        'chNum'  : [spwDic[spw]['chNum'] for spw in spwsInScan],
        'freq'   : [spwDic[spw]['freq'] for spw in spwsInScan],
        'Trx'    : TRX[:,:,data_index],
        'Tsys'   : TSY[:,:,data_index]}
    text_sd = 'Trx  : '
    for spw_index, spw in enumerate(spwsInScan): text_sd = text_sd + ' SPW%03d  %6.1f GHz |' % (spw, np.median(spwDic[spw]['freq'])* 1.0e-9)
    logFile.write(text_sd + '\n'); print(text_sd)
    text_sd = ' Pol : '
    for spw_index, spw in enumerate(spwsInScan): text_sd = text_sd + '     X         Y    |'
    logFile.write(text_sd + '\n'); print(text_sd)
    text_sd =  ' ----:-'
    for spw_index, spw in enumerate(spwsInScan): text_sd = text_sd + '--------------------+'
    logFile.write(text_sd + '\n'); print(text_sd)
    prevAnt, prevSPW = '', -1
    text_sd = ''
    for index in data_index:
        ant = antList[antID[index]]
        spw = spwID[index]
        for pol_index in [0,1]:
            text_sd = text_sd + '%7.1f K ' % (np.median(TRX[pol_index][:,index], axis=0))
        text_sd = text_sd + '|'
        if spw == spwsInScan[-1]:
            text_sd = '%s : ' % (ant) + text_sd
            logFile.write(text_sd + '\n'); print(text_sd)
            text_sd = ''
logFile.close()
if PLOTTSYS: plotTsysDic(prefix, TsysDic)
