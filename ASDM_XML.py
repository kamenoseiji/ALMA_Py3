# Module to read XML file in prefix
#
import os
import xml.etree.ElementTree as ET
import numpy as np
from interferometry import indexList
def CheckCorr( prefix ):
    Corr_XML = prefix + '/CorrelatorMode.xml'
    tree = ET.parse(Corr_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #---- Check by Correlat name
        for CorrID in row.findall('correlatorName'): CorrName = CorrID.text
    #
    return CorrName
#
def BandList( prefix ):
    KEY_ALMA ='ALMA_RB'
    RXList = []
    #-------- Check Receiver bands
    RB_XML = prefix + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Find Receiver band name
        for RX in row.findall('frequencyBand'): RXname = RX.text
        if KEY_ALMA in RXname:  RXList = RXList + [RXname.replace(KEY_ALMA, 'RB')]
    #
    return list(set(RXList))
#
def SPW_FULL_RES( prefix ):   # Dictionary for FULL_RES SPWs
    SBsign = {'LSB': -1, 'USB':+1, 'DSB':1}
    #-------- BB-SPWID connection
    SPW_XML = prefix + '/' + 'SpectralWindow.xml'
    tree = ET.parse(SPW_XML)
    root = tree.getroot()
    spwList = []
    spwDic  = dict(zip(spwList, [[]]*len(spwList))) # Dictionary for spw info
    #-------- SPWs with full resolution
    for row in root.findall('row'):
        #---- Check by BB name
        for entry in row.findall('basebandName'): BBname = entry.text
        if BBname == 'NOBB': continue
        BBID = int(BBname.split('_')[1])
        #---- Check by SPW ID
        for entry in row.findall('name'): SPWname = entry.text
        if 'FULL_RES' not in SPWname: continue
        for entry in row.findall('spectralWindowId'): spwID = entry.text; spw = int(entry.text.split('_')[1])
        #---- Check by Band name
        for entry in row.findall('name'): bandName = 'RB_' + entry.text.split('RB_')[1][0:2]
        for entry in row.findall('numChan'): chNum   = int(entry.text)
        for entry in row.findall('refFreq'): refFreq = float(entry.text)
        for entry in row.findall('totBandwidth'): BW = float(entry.text)
        for entry in row.findall('chanFreqStart'): ch0 = float(entry.text)
        for entry in row.findall('chanFreqStep'): chanFreqStep = float(entry.text)
        for entry in row.findall('refChan'): refChan = float(entry.text)
        spwDic[spw] = {
            'ID'     : spwID,
            'Band'   : bandName,
            'BB'     : BBID,
            'chNum'  : chNum,
            'refFreq': refFreq,
            'BW'     : BW,
            'ch0'    : ch0,
            'chStep' : chanFreqStep,
            'refChan': refChan
        }
    #
    RB_XML = prefix + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Avoid WVR and SQLD
        for entry in row.findall('name'): spwName = entry.text
        if 'WVR' in spwName : continue
        for entry in row.findall('receiverSideband'): SBname = entry.text
        if 'NOSB' not in SBname and 'TSB' not in SBname: continue
        #-------- Identify BB from SPWID
        for entry in row.findall('frequencyBand'):
            if '_' not in entry.text: continue
            BandID = int(entry.text.split('_')[-1]) 
        for entry in row.findall('freqLO'): LO1, LO2, LO3 = float(entry.text.split()[-3]), float(entry.text.split()[-2]),float(entry.text.split()[-1])
        for entry in row.findall('sidebandLO'): SB1, SB2, SB3 = SBsign[entry.text.split()[-3]], SBsign[entry.text.split()[-2]], SBsign[entry.text.split()[-1]]
        refFreq = SB1*SB2*SB3*( SB1*LO1 - SB2*LO2 + SB3*LO3 )
        for entry in row.findall('spectralWindowId'): spwID = int(entry.text.split('_')[1])
        spwList = [spw for spw in spwDic.keys() if abs(spwDic[spw]['refFreq'] - refFreq) < spwDic[spw]['BW']/64]
        if spwID in spwList:
            spwDic[spwID]['LO1'], spwDic[spwID]['LO2'], spwDic[spwID]['LO3'] = LO1, LO2, LO3
            spwDic[spwID]['SB1'], spwDic[spwID]['SB2'], spwDic[spwID]['SB3'] = SB1, SB2, SB3
    #
    RB_XML = prefix + '/' + 'Polarization.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        for entry in row.findall('numCorr'): polNum = int(entry.text)
        for spwID in spwDic.keys(): spwDic[spwID]['polNum'] = polNum
    return spwDic
#
def SPW_LO(spwDic, prefix):
    SBsign = {'LSB': -1, 'USB':+1, 'DSB':1}
    RB_XML = prefix + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        for entry in row.findall('spectralWindowId'): spwID = entry.text
        for entry in row.findall('name'): spwName = entry.text
        if 'WVR' in spwName : continue
        for entry in row.findall('freqLO'): LOList = entry.text.split()
        if int(LOList[1]) < 3: continue
        LO1, LO2, LO3 = float(entry.text.split()[-3]), float(entry.text.split()[-2]),float(entry.text.split()[-1])
        for entry in row.findall('sidebandLO'): SB1, SB2, SB3 = SBsign[entry.text.split()[-3]], SBsign[entry.text.split()[-2]], SBsign[entry.text.split()[-1]]
        for spw in spwDic.keys():
            if spwDic[spw]['ID'] == spwID:
                spwDic[spw]['LO1'], spwDic[spw]['LO2'], spwDic[spw]['LO3'] = LO1, LO2, LO3
                spwDic[spw]['SB1'], spwDic[spw]['SB2'], spwDic[spw]['SB3'] = SB1, SB2, SB3
    return spwDic
#
def spwMS(msfile):
    from casatools import msmetadata as msmdtool
    msmd = msmdtool()
    msmd.open(msfile)
    spwList = list(set(msmd.tdmspws()) | set(msmd.fdmspws()))
    spwList = [spw for spw in spwList if len(msmd.scansforspw(spw)) > 0]
    spwList.sort()
    spwDic  = dict(zip(spwList, [[]]*len(spwList))) # Dictionary for spw info
    for spw in spwDic.keys():
        chFreq = msmd.chanfreqs(spw)
        spwDic[spw] = {
            'Band'   : [spwName.split('ALMA_')[1] for spwName in msmd.namesforspws(spw)[0].split('#') if 'ALMA_' in spwName][0],
            'BB'     : msmd.baseband(spw)-1,
            'chNum'  : len(chFreq),
            'refFreq': msmd.reffreq(spw)['m0']['value'],
            'BW'     : msmd.bandwidths(spw),
            'ch0'    : chFreq[0],
            'chStep' : np.median(msmd.chanwidths(spw)),
            'refChan': -0.5,
            'scanList': msmd.scansforspw(spw).tolist() 
        }
    msmd.close()
    return spwDic
def spwIDMS(spwDic, msfile):
    if not os.path.isdir(msfile): return spwDic
    from casatools import msmetadata as msmdtool
    msmd = msmdtool()
    msmd.open(msfile)
    spwsMS = list(set(msmd.tdmspws()) | set(msmd.fdmspws()))
    spwsMS.sort()
    for spw in spwsMS:
        BBid = msmd.baseband(spw)
        refFreqMS = msmd.reffreq(spw)['m0']['value']
        asdmSPWs = [asdmSPW for asdmSPW in spwDic.keys() if spwDic[asdmSPW]['BB'] == BBid and spwDic[asdmSPW]['refFreq'] == refFreqMS]
        if len(asdmSPWs) == 1: spwDic[spw] = spwDic.pop(asdmSPWs[0])
    for spw in spwDic.keys():
        spwDic[spw]['scanList'] = msmd.scansforspw(spw).tolist() 
    msmd.close()
    return spwDic
def SourceList( prefix ):
    SrcList, PosList,  = [], []
    SRC_XML = prefix + '/' + 'Source.xml'
    tree = ET.parse(SRC_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Find Src Name
        for posCoordinate in row.findall('direction'): PosList = PosList + [[float(posCoordinate.text.strip().split()[2]), float(posCoordinate.text.strip().split()[3])]]
        for srcName in row.findall('sourceName'): SrcList = SrcList + [srcName.text.strip()]
    #
    return SrcList, PosList
#
def PolScan( prefix ):
    scanList, SrcList, STList, ETList = [], [], [], []
    KEY_ALMA ='CALIBRATE_POLARIZATION'
    #-------- Check Scan XML
    SC_XML = prefix + '/' + 'Scan.xml'
    tree = ET.parse(SC_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Find Scan Intent
        for ScanIntent in row.findall('scanIntent'):
            SCname = ScanIntent.text
            #print(SCname)
            if KEY_ALMA in SCname:
                for scanNumber in row.findall('scanNumber'): scanList = scanList + [int(scanNumber.text)]
                for StartTime in row.findall('startTime'):   STList   = STList   + [int(StartTime.text)]
                for EndTime   in row.findall('endTime'):     ETList   = ETList   + [int(EndTime.text)]
                for source in row.findall('sourceName'):     SrcList  = SrcList  + [source.text.strip()]
            #
        #
    #
    return scanList, SrcList, STList, ETList
#
