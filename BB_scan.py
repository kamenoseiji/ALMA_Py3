#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix = 'uid___A002_X106ccd1_Xe59d'
antName = 'DV20'
spwID = 13
scanList = [2 , 4 , 6 , 13, 16]
#-------- Functions
import math
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
def kmedianATM(power):    # clustering sky/amb/hot
    num_sample = len(power)
    init_index = [list(range(int(num_sample/3))), list(range(int(num_sample/3), int(2*num_sample/3))), list(range(int(2*num_sample/3), num_sample))]
    prev_index = init_index
    while True:
        new_index = [[],[],[]]
        prev_center = [np.median(power[prev_index[0]]), np.median(power[prev_index[1]]), np.median(power[prev_index[2]]) ]
        for index in list(range(num_sample)):
            nearest_index = np.argmin( abs( power[index] - np.array(prev_center) ) )
            new_index[nearest_index] = new_index[nearest_index] + [index]
        #
        new_center = [np.median(power[new_index[0]]), np.median(power[new_index[1]]), np.median(power[new_index[2]]) ]
        if new_center == prev_center: break
        prev_index = new_index
    #
    return np.array(new_center)
#
#-------- Main process
msfile = prefix + '.ms'
antList = GetAntName(msfile)
ant_index = np.where(antList == antName)[0][0]
pp = PdfPages('SQLD_%s_%s_SPW%d.pdf' % (prefix, antName, spwID))
figBB = plt.figure(0, figsize=(8,11))
figBB.suptitle('%s %s SQLD Power SPW=%d' % (prefix, antName, spwID))
scanNum = len(scanList)
panelNum_X = int(math.sqrt(scanNum + 1))
panelNum_Y = math.ceil( scanNum / panelNum_X )
for scan_index, scanID in enumerate(scanList):
    timeStamp, BBPower = GetPSpecScan(msfile, ant_index, spwID, scanID)
    medianPowerX, medianPowerY = kmedianATM(BBPower[0,0]), kmedianATM(BBPower[1,0])
    text_X = 'Pol-X Po:Pa:Ph=%.3f:1:%.3f' % ( medianPowerX[0]/medianPowerX[1],medianPowerX[2]/medianPowerX[1] )
    text_Y = 'Pol-Y Po:Pa:Ph=%.3f:1:%.3f' % ( medianPowerY[0]/medianPowerY[1],medianPowerY[2]/medianPowerY[1] )
    DT = []
    for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
    #
    BBPL = figBB.add_subplot( panelNum_Y, panelNum_X, scan_index + 1)
    BBPL.step( DT, BBPower[0,0], where='mid'); BBPL.plot( DT, BBPower[0,0], 'b.', label=text_X)
    BBPL.step( DT, BBPower[1,0], where='mid'); BBPL.plot( DT, BBPower[1,0], 'g.', label=text_Y)
    BBPL.legend(loc='lower right', prop={'size' :7}, numpoints=1)
    BBPL.tick_params(axis='both', labelsize=6)
    BBPL.set_title('SCAN %d' % scanID)
#
figBB.savefig(pp, format='pdf')
plt.close('all')
pp.close()



