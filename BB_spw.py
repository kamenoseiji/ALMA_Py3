#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix = 'uid___A002_X10d8ccb_X3c7f'
antName = 'DV16'
scanID = 2
spwList = [13, 1, 14, 3]
#-------- Functions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
def GetAntName(msfile):
    tb.open(msfile+'/'+'ANTENNA')
    namelist = tb.getcol("NAME")
    tb.close()
    return namelist
#
def GetPSpecScan(msfile, ant, spwID, scanID):
    data_desc_id = SPW2DATA_DESC_ID(msfile, spwID)
    Out='ANTENNA1 == %d && ANTENNA2 == %d && DATA_DESC_ID == %d && SCAN_NUMBER == %d' % (ant, ant, data_desc_id, scanID)
    tb.open(msfile)
    antXantYspw = tb.query(Out)
    colNameList = antXantYspw.colnames()
    colName = 'DATA'
    if 'FLOAT_DATA' in colNameList: colName = 'FLOAT_DATA'
    timeXY = antXantYspw.getcol('TIME')
    dataXY = antXantYspw.getcol(colName)
    tb.close()
    return timeXY, dataXY.real
#
def SPW2DATA_DESC_ID(msfile, spwID):
    tb.open(msfile + '/' + 'DATA_DESCRIPTION')
    data_desc_id, spw_index = -1,-1
    while( spw_index !=  spwID ):
        data_desc_id += 1
        spw_index = tb.getcell('SPECTRAL_WINDOW_ID', data_desc_id)
    tb.close()
    return data_desc_id
#
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
DT = []
antList = GetAntName(msfile)
ant_index = np.where(antList == antName)[0][0]
pp = PdfPages('SQLD_%s_%s_Scan%d.pdf' % (prefix, antName, scanID))
figBB = plt.figure(0, figsize=(8,11))
figBB.suptitle('%s %s SQLD Power Scan=%d' % (prefix, antName, scanID))
for spw_index, spwID in enumerate(spwList):
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
    if len(DT) == 0:
        for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
    #
    BBPL = figBB.add_subplot( 2, int((len(spwList) + 1) / 2) , spw_index + 1)
    BBPL.step( DT, BBPower[0,0], where='mid'); BBPL.plot( DT, BBPower[0,0], 'bo', label=text_X)
    BBPL.step( DT, BBPower[1,0], where='mid'); BBPL.plot( DT, BBPower[1,0], 'go', label=text_Y)
    BBPL.legend(loc='best', prop={'size' :7}, numpoints=1)
    BBPL.tick_params(axis='both', labelsize=6)
    BBPL.set_title('SPW %d' % spwID)
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()
