import os
import numpy as np
from interferometry import ALMA_lat
from casatools import msmetadata as msmdtool
msmd = msmdtool()
#
#-------- Calcurate PA
def AzEl2PA(az, el, lat=ALMA_lat): # Azimuth, Elevation, Latitude (default=ALMA) in [rad]
    cos_lat, sin_lat = np.cos(lat), np.sin(lat)
    return np.arctan2( -cos_lat* np.sin(az), (sin_lat* np.cos(el) - cos_lat* np.sin(el)* np.cos(az)) )
#
#-------- Get Stokes Parameters from AMAPOLA
def GetAMAPOLAStokes(R_DIR, SCR_DIR, sourceList, timeText, FreqGHz):    # 
    # sourceList : e.g. ['J1256-0547', 'J1924-2914']
    # timeText   : e.g. '2017/04/12/11:26:17'
    # FreqGHz    : frequency in GHz
    StokesDic = dict(zip(sourceList, [[]]*len(sourceList)))   # Stokes parameters for each source
    os.system('rm -rf CalQU.data')
    text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, timeText, FreqGHz)
    for source in sourceList: text_sd = text_sd + ' ' + source
    print(text_sd)
    os.system(text_sd)
    fp = open('CalQU.data')
    lines = fp.readlines()
    fp.close()
    for eachLine in lines:
        sourceName = eachLine.split()[0]
        StokesDic[sourceName] = [float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
    #
    return StokesDic
#
#-------- PA and polarization responses
def PolResponse(msfile, StokesDic, BandPA, scanList, AzScanList, ElScanList):
    PAList, CSList, SNList, QCpUSList, UCmQSList = [], [], [], [], []
    msmd.open(msfile)
    #-------- Check AZEL
    print('        Source     :    I     p%     EVPA  QCpUS  UCmQS')
    print('-------+-----------+-------+------+------+------+------')
    for scan_index, scan in enumerate(scanList):
        AzScan, ElScan = AzScanList[scan_index], ElScanList[scan_index]
        PA = AzEl2PA(AzScan, ElScan) + BandPA
        CS, SN, QCpUS, UCmQS = np.cos(2.0* PA), np.sin(2.0* PA), np.zeros(len(PA)), np.zeros(len(PA))
        sourceName = list(StokesDic.keys())[msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])]
        if StokesDic[sourceName] != []:
            QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
            print('Scan%3d %s : %6.2f %6.2f %6.2f %6.2f %6.2f' % (scan, sourceName, StokesDic[sourceName][0], 100.0*np.sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0], 90.0* np.arctan2(StokesDic[sourceName][2], StokesDic[sourceName][1])/np.pi, np.median(QCpUS), np.median(UCmQS)))
        #
        PAList, CSList, SNList, QCpUSList, UCmQSList = PAList + [PA], CSList + [CS], SNList + [SN], QCpUSList + [QCpUS], UCmQSList + [UCmQS]
    msmd.close(); msmd.done()
    return PAList, CSList, SNList, QCpUSList, UCmQSList

