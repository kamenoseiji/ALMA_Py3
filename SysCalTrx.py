import numpy as np
from Plotters import plotTsysDic
import analysisUtils as au
from interferometry import GetAntName, GetChNum, GetBandNames, GetAzEl, GetTemp, RADDEG, Tcmb
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
prefix = 'uid___A002_X13a3180_X171db'
PLOTTSYS = True
'''
#-------- Read SYSCAL table
msfile = prefix + '.ms'
tb.open(msfile + '/SYSCAL')
colnameList = tb.colnames()
timeStamp = tb.getcol('TIME')
tempAtm = GetTemp(msfile)
timeSyscalList = timeStamp[ np.where(np.diff(timeStamp) != 0.0)[0].tolist() + [-1] ]
antID = tb.getcol('ANTENNA_ID')
spwID = tb.getcol('SPECTRAL_WINDOW_ID')
TRX = tb.getcol('TRX_SPECTRUM')
TSY = tb.getcol('TSYS_SPECTRUM')
tb.close()
antList, spwList = GetAntName(prefix + '.ms'), np.unique(spwID)
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelRefIndex = np.where(AntID == np.unique(AntID)[0])[0].tolist()
azelTime, AZ, EL = azelTime[azelRefIndex], AZ[azelRefIndex], EL[azelRefIndex]
spwDic = dict(zip(spwList, [[]]*len(spwList)))  # SPW
antNum, spwNum = len(antList), len(spwList)
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
    azelRefIndex = np.argmin(abs(azelTime - timeList[scan_index]))
    TsysDic[scan] = {
        'Band'   : spwDic[spwsInScan[0]]['band'],
        'antList': antList[antsInScan],
        'mjdSec' : timeList[scan_index],
        'Az'     : AZ[azelRefIndex],
        'El'     : EL[azelRefIndex],
        'SecZ'   : 1.0/np.sin(EL[azelRefIndex]),
        'spwList': spwsInScan,
        'chNum'  : [spwDic[spw]['chNum'] for spw in spwsInScan],
        'freq'   : [spwDic[spw]['freq']*1.0e-9 for spw in spwsInScan],
        'Trx'    : TRX[:,:,data_index],
        'Tsys'   : TSY[:,:,data_index]}
for bandName in np.unique(BandNames):
    logFile  = open('%s-%s-TelCal.log' % (prefix, bandName), 'w')
    bandScanList = [scan for scan in TsysDic.keys() if TsysDic[scan]['Band'] == bandName]
    TsysBandDic = {keyScan: dicValue for keyScan, dicValue in TsysDic.items() if keyScan in bandScanList}
    #-------- Save into NPY
    spwNum = len(TsysBandDic[bandScanList[0]]['spwList'])
    np.save('%s-%s.TrxAnt.npy' % (prefix, bandName), antList)     # antList[ant]
    for spw_index, spw in enumerate(TsysBandDic[bandScanList[0]]['spwList']):
        np.save('%s-%s-SPW%d.TrxFreq.npy' % (prefix, bandName, spw), TsysBandDic[bandScanList[0]]['freq'][spw_index])
        TrxList, TskyList = [], []
        for scan_index, scan in enumerate(TsysBandDic.keys()):
            Trx  = TsysBandDic[scan]['Trx'][:,:,range(spw_index, antNum*spwNum, spwNum)]
            Tsys = TsysBandDic[scan]['Tsys'][:,:,range(spw_index, antNum*spwNum, spwNum)]
            TrxList  = TrxList  + [Trx]
            TskyList = TskyList + [np.mean((Tsys - Trx)* tempAtm / (Tsys + tempAtm), axis=0) ]
        np.save('%s-%s-SPW%d.Trx.npy' % (prefix, bandName, spw), np.array(TrxList).transpose(1,2,3,0))
        np.save('%s-%s-SPW%d.Tsky.npy' % (prefix, bandName, spw), np.array(TskyList).transpose(1,2,0))
    #-------- Log Trx 
    atmTimeList, atmElList = [], []
    for scan_index, scan in enumerate(TsysBandDic.keys()):
        atmTimeList = atmTimeList + [TsysBandDic[scan]['mjdSec']]
        atmElList   = atmElList + [TsysBandDic[scan]['El']]
        timeLabel = 'Scan %d : %s EL=%.1f' % (scan, au.call_qa_time('%fs' % (TsysBandDic[scan]['mjdSec']), form='fits'), RADDEG*TsysBandDic[scan]['El'])
        logFile.write(timeLabel + '\n'); print(timeLabel)
        text_sd = 'Trx  : '
        for spw_index, spw in enumerate(TsysBandDic[scan]['spwList']):
            text_sd = text_sd + ' SPW%03d  %6.1f GHz |' % (spw, np.median(TsysBandDic[scan]['freq'][spw_index]))
        logFile.write(text_sd + '\n'); print(text_sd)
        text_sd = ' Pol : '
        for spw_index, spw in enumerate(TsysBandDic[scan]['spwList']): text_sd = text_sd + '     X         Y    |'
        logFile.write(text_sd + '\n'); print(text_sd)
        text_sd =  ' ----:-'
        for spw_index, spw in enumerate(TsysBandDic[scan]['spwList']): text_sd = text_sd + '--------------------+'
        logFile.write(text_sd + '\n'); print(text_sd)
        for ant_index, ant in enumerate(TsysBandDic[scan]['antList']):
            text_sd = '%s : ' % (ant)
            for spw_index, spw in enumerate(TsysBandDic[scan]['spwList']):
                data_index = spwNum* ant_index + spw_index
                for pol_index in [0,1]: text_sd = text_sd + '%7.1f K ' % (np.median(TsysBandDic[scan]['Trx'][pol_index][:,data_index], axis=0))
                text_sd = text_sd + '|'
            logFile.write(text_sd + '\n'); print(text_sd)
    np.save('%s-%s.atmTime.npy' % (prefix, bandName), np.array(atmTimeList))
    np.save('%s-%s.EL.npy' % (prefix, bandName), np.array(atmElList))
    if PLOTTSYS: plotTsysDic(prefix, TsysBandDic)
    logFile.close()
