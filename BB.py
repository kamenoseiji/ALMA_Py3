import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import indexList, GetAntName, GetPSpecScan, SPW2DATA_DESC_ID
import datetime
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='SPW ID   e.g. 13', default='13')
parser.add_option('-c', dest='scanID', metavar='scanID',
    help='Scan ID  e.g. 2', default='2')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
scanID  = int(options.scanID)
spwID   = int(options.spw)
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
#-------- Parameters
msfile = prefix + '.ms'
msmd.open(msfile)
timeOFF, timeON, timeAMB, timeHOT, timeTEST = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#TEST")
subscanDic = dict( zip(['OFF', 'ON', 'AMB', 'HOT', 'TEST'], [[]]* 5))
#-------- Main process
DT = []
antList = GetAntName(msfile)
antNum = len(antList)
pp = PdfPages('SQLD_%s_SPW%d_Scan%d.pdf' % (prefix, spwID, scanID))
figBB = plt.figure(0, figsize=(8,11))
figBB.suptitle('%s SQLD Power SPW=%d Scan=%d' % (prefix, spwID, scanID))
for ant_index, antName in enumerate(antList):
    timeStamp, BBPower = GetPSpecScan(msfile, ant_index, spwID, scanID)
    if len(DT) == 0:
        for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
    #
    #-------- Subscan Intents
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
    #-------- Plot
    BBPL = figBB.add_subplot( int((antNum + 3)/4), 4, ant_index + 1)
    BBPL.step( DT, BBPower[0,0], where='mid'); BBPL.plot( DT, BBPower[0,0], 'b.', label=text_X)
    BBPL.step( DT, BBPower[1,0], where='mid'); BBPL.plot( DT, BBPower[1,0], 'g.', label=text_Y)
    BBPL.legend(loc='best', prop={'size' :3}, numpoints=1)
    BBPL.set_title(antName, fontsize=5, pad=-14)
    BBPL.tick_params(axis='both', labelsize=4)
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()
