import os
import numpy as np
from interferometry import ALMA_lat, GetAzEl, AzElMatch, AzEl2PA
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
        I, Q, U = float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3])
        if abs(Q) > 0.5*I: Q = 0.0
        if abs(U) > 0.5*I: U = 0.0
        StokesDic[sourceName] = [FreqGHz, I, Q, U, 0.0]
    #
    return StokesDic
#
def GetSSOFlux(StokesDic, timeText, FreqGHz):
    from Grid import SSOCatalog
    # StokesDic  : Stokes parameter dictionlary
    # timeText   : e.g. '2017/04/12/11:26:17'
    # FreqGHz    : frequency list in GHz
    #sourceList = list(StokesDic.keys())
    #SSOList = [source for source in sourceList if not str.isdigit(source[1])]
    SSOList = [source for source in StokesDic.keys() if source in SSOCatalog]
    SSODic = dict(zip(SSOList, [[]]*len(SSOList)))
    if len(SSOList) == 0: return StokesDic, SSODic
    for SSO in SSOList:
        flux, major, minor, pa = [], [], [], []
        for spw_index, freq in enumerate(FreqGHz):
            freqcorr = 1.0
            if freq < 62.0:
                freqcorr = (freq/62.0)**2
                freq = 62.0
                #print('%.1f GHz : corr=%.1f\n' % (freq, freqcorr))
            text_Freq = '%6.2fGHz' % (freq)
            SSOmodel = au.predictcomp(objname=SSO, standard="Butler-JPL-Horizons 2012", minfreq=text_Freq, maxfreq=text_Freq, nfreqs=1, prefix="", antennalist="aca.cycle3.cfg", epoch=timeText, showplot=True)
            flux  = flux + [freqcorr* SSOmodel['spectrum']['bl0flux']['value']]
        #
        major = SSOmodel['shape']['majoraxis']['value']* np.pi / 21600.0
        minor = SSOmodel['shape']['minoraxis']['value']* np.pi / 21600.0
        pa    = SSOmodel['shape']['positionangle']['value']* np.pi / 180.0
        StokesDic[SSO] = [FreqGHz, flux]
        SSODic[SSO]    = [FreqGHz, flux, [major, minor, pa]]
    #
    return StokesDic, SSODic
#-------- PA and polarization responses
def PolResponse(msfile, srcDic, StokesDic, BandPA, scanList, mjdList):
    scanDic = dict(zip(scanList, [[]]* len(scanList)))
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    msmd.open(msfile)
    #-------- Check AZEL
    print('------------------- AMAPOLA-based prediction -----------------')
    print('        Source     :    I     p%     EVPA  QCpUS  UCmQS   EL  ')
    print('-------+-----------+-------+------+------+------+------+------')
    for scan_index, scan in enumerate(scanList):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])
        sourceName = srcDic[sourceID]['Name']
        probeAntID = 0
        while True:
            AzScan, ElScan = AzElMatch(mjdList[scan_index], azelTime, AntID, probeAntID, AZ, EL)
            if np.min(ElScan) < 0.1:
                probeAntID += 1
            else:
                 break
        #
        PA = AzEl2PA(AzScan, ElScan) + BandPA
        CS, SN, QCpUS, UCmQS = np.cos(2.0* PA), np.sin(2.0* PA), np.zeros(len(PA)), np.zeros(len(PA))
        StokesI = 0.0
        if str.isdigit(sourceName[1]) and len(StokesDic[sourceName]) > 0:
            StokesI = StokesDic[sourceName][1]
            QCpUS = StokesDic[sourceName][2]*CS + StokesDic[sourceName][3]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][3]*CS - StokesDic[sourceName][2]*SN   # Ucos - Qsin
            BPquality = StokesDic[sourceName][1]* np.sin(np.median(ElScan) - 0.5)  # cut at 40 deg
            XYquality = np.median(abs(QCpUS))* np.sin(np.median(ElScan) - 0.5)     # cut at 30 deg
            print('Scan%3d %s : %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (scan, sourceName, StokesDic[sourceName][1], 100.0*np.sqrt(StokesDic[sourceName][2]**2 + StokesDic[sourceName][3]**2)/StokesDic[sourceName][1], 90.0* np.arctan2(StokesDic[sourceName][3], StokesDic[sourceName][2])/np.pi, np.median(QCpUS), np.median(UCmQS), 180.0* np.median(ElScan)/np.pi))
        #
        scanDic[scan] = {
            'msfile': msfile,
            'source': sourceName,
            'mjdSec': mjdList[scan_index],
            'EL'    : ElScan,
            'PA'    : PA,
            'SA'    : srcDic[sourceID]['SA'],
            'I'     : StokesI,
            'QCpUS' : QCpUS,
            'UCmQS' : UCmQS}
    msmd.close(); msmd.done()
    return scanDic
#
