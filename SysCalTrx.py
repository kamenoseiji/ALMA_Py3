import numpy as np
import analysisUtils as au
from interferometry import GetAntName, GetChNum
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
#-------- Read SYSCAL table
msfile = prefix + '.ms'
tb.open(msfile + '/SYSCAL')
colnameList = tb.colnames()
antID = tb.getcol('ANTENNA_ID')
timeStamp = tb.getcol('TIME')
spwID = tb.getcol('SPECTRAL_WINDOW_ID')
TRX = tb.getcol('TRX_SPECTRUM')
tb.close()
timeList, antList, spwList = np.unique(timeStamp), GetAntName(prefix + '.ms'), np.unique(spwID)
scanNum, antNum, spwNum = len(timeList), len(antList), len(spwList)
logFile  = open(prefix + '-TelCal.log', 'w')
msmd.open(msfile)
scanList = msmd.scansforintent("CALIBRATE_ATMOSPHERE*")
timeAMB  = msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT")
msmd.close()
timeList = timeAMB[(np.where(np.diff(timeAMB) > 10*np.median(np.diff(timeAMB)))[0] - 1).tolist() + [-1]].tolist()
#-------- Read SPW frequency
freqList = []
for spw in spwList:
    chNum, chWid, freq = GetChNum(msfile, spw)
    freqList = freqList + [np.median(freq)* 1.0e-9]
#
for scan_index, scan in enumerate(scanList):
    timeLabel = 'Scan %d : %s' % (scan, au.call_qa_time('%fs' % (timeList[scan_index]), form='fits'))
    logFile.write(timeLabel + '\n'); print(timeLabel)
    text_sd = 'Trx  : '
    for spw_index, spw in enumerate(spwList): text_sd = text_sd + ' SPW%03d  %6.1f GHz |' % (spw, freqList[spw_index])
    logFile.write(text_sd + '\n'); print(text_sd)
    text_sd = ' Pol : '
    for spw_index, spw in enumerate(spwList): text_sd = text_sd + '     X         Y    |'
    logFile.write(text_sd + '\n'); print(text_sd)
    text_sd =  ' ----:-'
    for spw_index, spw in enumerate(spwList): text_sd = text_sd + '--------------------+'
    logFile.write(text_sd + '\n'); print(text_sd)
    for ant_index, ant in enumerate(antList):
        text_sd = '%s : ' % (ant)
        for spw_index, spw in enumerate(spwList):
            for pol_index in [0,1]:
                index = spw_index + spwNum* (ant_index + antNum* scan_index)
                text_sd = text_sd + '%7.1f K ' % (np.median(TRX[pol_index], axis=0)[index])
            text_sd = text_sd + '|'
        logFile.write(text_sd + '\n'); print(text_sd)
logFile.close()
