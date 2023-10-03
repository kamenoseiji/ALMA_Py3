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
#-------- Parameters
prefix = 'uid___A002_Xff81e8_X253'
msfile = prefix + '.ms'
scanID = 9
spwID = 0
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
    BBPL = figBB.add_subplot( int((antNum + 3)/4), 4, ant_index + 1)
    BBPL.step( DT, BBPower[0,0], where='mid'); BBPL.plot( DT, BBPower[0,0], 'o', label=antName + 'Pol X')
    BBPL.step( DT, BBPower[1,0], where='mid'); BBPL.plot( DT, BBPower[1,0], 'o', label=antName + 'Pol Y')
    BBPL.legend(loc='lower right', prop={'size' :7}, numpoints=1)
    BBPL.tick_params(axis='both', labelsize=6)
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()

