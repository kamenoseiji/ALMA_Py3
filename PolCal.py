import os
import numpy as np
from interferometry import ALMA_lat
import analysisUtils as au
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
        StokesDic[sourceName] = [FreqGHz, float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
    #
    return StokesDic
#
def GetSSOFlux(StokesDic, timeText, FreqGHz):
    # StokesDic  : Stokes parameter dictionlary
    # timeText   : e.g. '2017/04/12/11:26:17'
    # FreqGHz    : frequency list in GHz
    sourceList = list(StokesDic.keys())
    SSOList = [source for source in sourceList if not str.isdigit(source[1])]
    SSODic = dict(zip(SSOList, [[]]*len(SSOList)))
    for SSO in SSOList:
        flux, major, minor, pa = [], [], [], []
        for spw_index, freq in enumerate(FreqGHz):
            text_Freq = '%6.2fGHz' % (freq)
            SSOmodel = au.predictcomp(objname=SSO, standard="Butler-JPL-Horizons 2012", minfreq=text_Freq, maxfreq=text_Freq, nfreqs=1, prefix="", antennalist="aca.cycle3.cfg", epoch=timeText, showplot=True)
            flux  = flux + [SSOmodel['spectrum']['bl0flux']['value']]
        #
        major = SSOmodel['shape']['majoraxis']['value']* np.pi / 21600.0
        minor = SSOmodel['shape']['minoraxis']['value']* np.pi / 21600.0
        pa    = SSOmodel['shape']['positionangle']['value']* np.pi / 180.0
        StokesDic[SSO] = [FreqGHz, flux]
        SSODic[SSO]    = [FreqGHz, flux, [major, minor, pa]]
    #
    return StokesDic, SSODic
#-------- PA and polarization responses
def PolResponse(msfile, StokesDic, BandPA, scanList, AzScanList, ElScanList):
    PAList, CSList, SNList, QCpUSList, UCmQSList = [], [], [], [], []
    scanDic = dict(zip(scanList, [['', 0.0, 0.0, 0.0, 0.0]]*len(scanList))) # scanDict{ scanID: [source, EL, I, BPquality, XYquality]}
    msmd.open(msfile)
    #-------- Check AZEL
    print('------------------- AMAPOLA-based prediction -----------------')
    print('        Source     :    I     p%     EVPA  QCpUS  UCmQS   EL  ')
    print('-------+-----------+-------+------+------+------+------+------')
    for scan_index, scan in enumerate(scanList):
        AzScan, ElScan = AzScanList[scan_index], ElScanList[scan_index]
        PA = AzEl2PA(AzScan, ElScan) + BandPA
        CS, SN, QCpUS, UCmQS = np.cos(2.0* PA), np.sin(2.0* PA), np.zeros(len(PA)), np.zeros(len(PA))
        sourceName = list(StokesDic.keys())[msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])]
        scanDic[scan] = [sourceName, np.median(ElScan), 0.0, 0.0, 0.0]
        if str.isdigit(sourceName[1]) and len(StokesDic[sourceName]) > 0:
            QCpUS = StokesDic[sourceName][2]*CS + StokesDic[sourceName][3]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][3]*CS - StokesDic[sourceName][2]*SN   # Ucos - Qsin
            BPquality = StokesDic[sourceName][1]* np.sin(np.median(ElScan) - 0.5)  # cut at 40 deg
            XYquality = np.median(abs(QCpUS))* np.sin(np.median(ElScan) - 0.5)             # cut at 30 deg
            print('Scan%3d %s : %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (scan, sourceName, StokesDic[sourceName][1], 100.0*np.sqrt(StokesDic[sourceName][2]**2 + StokesDic[sourceName][3]**2)/StokesDic[sourceName][1], 90.0* np.arctan2(StokesDic[sourceName][3], StokesDic[sourceName][2])/np.pi, np.median(QCpUS), np.median(UCmQS), 180.0* np.median(ElScan)/np.pi))
            scanDic[scan] = [sourceName, np.median(ElScan), StokesDic[sourceName][1], BPquality, XYquality]
        #
        PAList, CSList, SNList, QCpUSList, UCmQSList = PAList + [PA], CSList + [CS], SNList + [SN], QCpUSList + [QCpUS], UCmQSList + [UCmQS]
    msmd.close(); msmd.done()
    return PAList, CSList, SNList, QCpUSList, UCmQSList, scanDic
#
