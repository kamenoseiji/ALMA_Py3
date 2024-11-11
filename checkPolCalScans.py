import os
SCR_DIR = '/home/skameno/ALMA_Py3/'
import sys
import numpy as np
import analysisUtils as au
from interferometry import BANDPA, BANDFQ, GetBPcalSPWs, GetSourceDic, indexList, GetAzEl, GetTimerecord, AzElMatch, AzEl2PA
from Grid import SSOCatalog, SSOscore, ELshadow
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='uid  e.g. uid___A002_X11adad7_X19ab0', default='')
parser.add_option('-Q', dest='QUMODEL', metavar='QUMODEL',
    help='Use a priori Stokes parameters', action="store_true")
parser.add_option('-B', dest='BPscans', metavar='BPscans',
    help='Bandpass scans for each band e.g. [3,6]', default='')
parser.add_option('-E', dest='EQscans', metavar='EQscans',
    help='Equalizing scans for each band e.g. [3,6]', default='')
#
(options, args) = parser.parse_args()
#
QUMODEL= options.QUMODEL
BPscans = [scan for scan in options.BPscans]
EQscans = [scan for scan in options.EQscans]
prefix = options.prefix
'''
prefix = '2023.1.00675.S_X10ed869_X5fa3'
QUMODEL = True
'''
msfile = prefix + '.ms'
#-------- Check SPWs of atmCal
logfile = open(prefix + '-PolQuery.log', 'w')
print('---Checking spectral windows for ' + prefix)
bpSPWs = GetBPcalSPWs(msfile)
msmd.open(msfile)
bpspwNames = msmd.namesforspws(bpSPWs)
bandNames = []
bandNamePattern = r'RB_..'
for spwName in bpspwNames :
    bandNames = bandNames + re.findall(bandNamePattern, spwName)
UniqBands = np.unique(bandNames).tolist(); NumBands = len(UniqBands)
bpspwLists, bpscanLists, BandPA = [], [], []
for band_index, BandName in enumerate(UniqBands):
    BandPA = BandPA + [(BANDPA[int(BandName[3:5])] + 90.0)*np.pi/180.0]
    bpspwLists  = bpspwLists  + [np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(bandNames))].tolist()]
    bpscanLists = bpscanLists + [msmd.scansforspw(bpspwLists[band_index][0]).tolist()]
    #print(' ', end='')
    text_sd = BandName + ' SPW'
    for spw in bpspwLists[band_index]: text_sd = text_sd + ' %d' % (spw)
    print(text_sd); logfile.write(text_sd + '\n')
#
#-------- Scan List
ONScans = np.sort(np.array(list(set(msmd.scansforintent("*CALIBRATE_AMPLI*")) | set(msmd.scansforintent("*CALIBRATE_BANDPASS*")) | set(msmd.scansforintent("*CALIBRATE_POLARIZATION*")) | set(msmd.scansforintent("*CALIBRATE_FLUX*")) | set(msmd.scansforintent("*CALIBRATE_PHASE*")) | set(msmd.scansforintent("*OBSERVE_CHECK_SOURCE*")) | set(msmd.scansforintent("*OBSERVE_TARGET*")) | set(msmd.scansforintent("*CALIBRATE_APPPHASE*")))))
#-------- Check source list
print('---Checking source list')
srcDic = GetSourceDic(msfile)
srcID  = list(set([msmd.fieldsforscan(scan)[0] for scan in ONScans])); srcID.sort
sourceList = [srcDic[ID]['Name'] for ID in srcID]; numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
msmd.close()
msmd.done()
#-------- Loop for Bands
for band_index, BandName in enumerate(UniqBands):
    msmd.open(msfile)
    bandID = int(BandName[3:5])
    ONScan = np.array(bpscanLists[band_index])[indexList( ONScans, np.array(bpscanLists[band_index]))]
    onsourceScans = ONScan.tolist()
    scanNum = len(onsourceScans)
    SSOScanIndex = []
    StokesDic = dict(zip(sourceList, [[]]*numSource))   # Stokes parameters for each source
    #-------- Check AZEL
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == 0 )[0].tolist() 
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    OnAZ, OnEL, OnPA, sourceScan, BPquality, EQquality, FLscore, refTime = [], [], [], [], [], [], np.zeros(scanNum), []
    #-------- Check QU catalog
    if QUMODEL: # A priori QU model from Rdata
        os.system('rm -rf CalQU.data')
        text_sd = 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (azelTime[0]), form='ymd')[0], BANDFQ[bandID])
        for source in sourceList: text_sd = text_sd + ' ' + source
        print(text_sd); os.system(text_sd)
        fp = open('CalQU.data')
        lines = fp.readlines()
        fp.close()
        for eachLine in lines:
            sourceName = eachLine.split()[0]
            StokesDic[sourceName] = [float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
        #
        for source in [source for source in StokesDic.keys() if len(StokesDic[source]) == 0]: del StokesDic[source]
    #
    text_sd = 'Scan Source        AZ   EL     PA   dPA   pRes   BPqual EQqual'
    print(text_sd); logfile.write(text_sd + '\n')
    for scan_index, scan in enumerate(onsourceScans):
        fieldID = msmd.fieldsforscan(scan)[0]
        interval, timeStamp = GetTimerecord(msfile, 0, 0, bpspwLists[band_index][0], scan)
        AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ = OnAZ + [np.median(AzScan)]
        OnEL = OnEL + [np.median(ElScan)]
        OnPA = OnPA + [np.median(PA)]
        refTime = refTime + [np.median(timeStamp)]
        sourceName = srcDic[fieldID]['Name']
        sourceScan = sourceScan + [sourceName]
        if sourceName not in StokesDic:
            QCpUS, UCmQS = 0.0, 0.0
            BPquality = BPquality + [-9999.9]
            EQquality = EQquality + [-9999.9]
            continue
        if len(StokesDic[sourceName]) == 4:
            CS, SN = np.cos(2.0* OnPA[scan_index]), np.sin(2.0* OnPA[scan_index])
            QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
            if QUMODEL:
                BPquality = BPquality + [1000.0* abs(UCmQS)* np.sin(OnEL[scan_index]) ] # / np.sqrt(StokesDic[sourceName][0])]
            else:
                BPquality = BPquality + [1000.0* abs(UCmQS)* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0])]
            EQquality = EQquality + [StokesDic[sourceName][0]**2 * np.sin(OnEL[scan_index] - ELshadow) / (1.0e-4 + abs(QCpUS))]
        else:
            QCpUS, UCmQS = 0.0, 0.0
            BPquality = BPquality + [-9999.9]
            EQquality = EQquality + [-9999.9]
        #

        text_sd = ' %3d %10s %+6.1f %4.1f %6.1f %5.2f %+5.3f %7.2f %+6.0f' % (onsourceScans[scan_index], sourceName, 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[-1], EQquality[-1])
        print(text_sd); logfile.write(text_sd + '\n')
        if fieldID in SSOList: FLscore[scan_index] = np.exp(np.log(np.sin(OnEL[scan_index])-0.34))* SSOscore[bandID-1][SSOCatalog.index(sourceList)]
    #-------- Select Bandpass Calibrator and Equalization Calibrator
    BPscanIndex = np.argmax(BPquality) if len(BPscans) == 0 else BPscans[band_index]
    EQscanIndex = np.argmax(EQquality) if len(EQscans) == 0 else EQscans[band_index]
    BPScan = onsourceScans[BPscanIndex]; BPcal = sourceScan[BPscanIndex]; BPEL = OnEL[BPscanIndex]
    EQScan = onsourceScans[EQscanIndex]; EQcal = sourceScan[EQscanIndex]; EQEL = OnEL[EQscanIndex]
    text_sd = '%s BPscan %d [%s EL=%4.1f]' % (BandName, BPScan, BPcal, 180.0* BPEL/np.pi); print(text_sd); logfile.write(text_sd + '\n')
    text_sd = '%s EQscan %d [%s EL=%4.1f]' % (BandName, EQScan, BPcal, 180.0* EQEL/np.pi); print(text_sd); logfile.write(text_sd + '\n')
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(srcDic.keys()) )
    if len(BandSSOList) == 0: Apriori = True
    #-------- Polarization setup
    scnspw = bpspwLists[band_index]; scnspwNum = len(scnspw)
    polNum = msmd.ncorrforpol(msmd.polidfordatadesc(scnspw[0]))
    PolList = ['X', 'Y']
    msmd.done()
#
logfile.close()
del msfile, UniqBands
