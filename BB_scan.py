import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import indexList, GetAntName, GetPSpecScan, SPW2DATA_DESC_ID
import datetime
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antName', metavar='antName',
    help='Antenna name', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List e.g. 2, 5, 10', default='')
parser.add_option('-s', dest='spwID', metavar='spwID',
    help='SPW ID   e.g. 13', default='13')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
antName = options.antName
spwID  =  int(options.spwID)
scanList = options.scanList.split(',')
scanList = [int(scan) for scan in scanList]
scanNum  = len(scanList)
#-------- Functions
def medianATM(subscanDic, power):    # clustering sky/amb/hot
    for subScan in subscanDic.keys():
        if len(subscanDic[subScan]) < 3:
            subscanDic[subScan] = [np.nan, np.nan]
        else:
            scanPower = [ np.median(power[0, subscanDic[subScan]]), np.median(power[1, subscanDic[subScan]]) ]
            subscanDic[subScan] = scanPower
        #
    #
    return subscanDic
#
#-------- Main process
msfile = prefix + '.ms'
msmd.open(msfile)
timeOFF, timeON, timeAMB, timeHOT, timeTEST = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#TEST")
subscanDic = dict( zip(['OFF', 'ON', 'AMB', 'HOT', 'TEST'], [[]]* 5))
antList = GetAntName(msfile)
ant_index = np.where(antList == antName)[0][0]
pp = PdfPages('SQLD_%s_%s_SPW%d.pdf' % (prefix, antName, spwID))
figBB = plt.figure(0, figsize=(8,11))
figBB.suptitle('%s %s SQLD Power SPW=%d' % (prefix, antName, spwID))
figBB.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/BB_scan.py to generate this plot', fontsize=4)
scanNum = len(scanList)
panelNum_X = int(math.sqrt(scanNum + 1))
panelNum_Y = math.ceil( scanNum / panelNum_X )
DT = []
for scan_index, scanID in enumerate(scanList):
    timeStamp, BBPower = GetPSpecScan(msfile, ant_index, spwID, scanID)
    subscanDic['OFF'], subscanDic['ON'],subscanDic['AMB'],subscanDic['HOT'],subscanDic['TEST'] = indexList(timeOFF, timeStamp), indexList(timeON, timeStamp), indexList(timeAMB, timeStamp), indexList(timeHOT, timeStamp), indexList(timeTEST, timeStamp)
    subscanDic = medianATM(subscanDic, BBPower[:,0])
    text_X, text_Y = 'Pol-X ', 'Pol-Y '
    for subScan in subscanDic.keys():
        if subscanDic[subScan] == [np.nan, np.nan] : continue
        text_X = text_X + subScan + ':'
        text_Y = text_Y + subScan + ':'
    #
    text_X, text_Y = text_X[:len(text_X)-1] + '=', text_Y[:len(text_Y)-1] + '='
    for subScan in subscanDic.keys():
        if subscanDic[subScan] == [np.nan, np.nan] : continue
        if subScan == 'AMB':
            text_X = text_X + '1:'
            text_Y = text_Y + '1:'
        else:
            text_X = text_X + '%.3f:' % (subscanDic[subScan][0] / subscanDic['AMB'][0])
            text_Y = text_Y + '%.3f:' % (subscanDic[subScan][1] / subscanDic['AMB'][1])
        #
    text_X, text_Y = text_X[:len(text_X)-1], text_Y[:len(text_Y)-1]
    if len(DT) != len(timeStamp):
        DT = []
        for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
    #
    BBPL = figBB.add_subplot( panelNum_Y, panelNum_X, scan_index + 1)
    BBPL.step( DT, BBPower[0,0], where='mid'); BBPL.plot( DT, BBPower[0,0], 'b.', label=text_X)
    BBPL.step( DT, BBPower[1,0], where='mid'); BBPL.plot( DT, BBPower[1,0], 'g.', label=text_Y)
    BBPL.legend(loc='best', prop={'size' :5}, numpoints=1)
    BBPL.tick_params(axis='both', labelsize=4)
    BBPL.set_title('SCAN %d' % scanID)
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()
