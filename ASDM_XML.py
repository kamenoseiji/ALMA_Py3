# Module to read XML file in prefix
#
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
        BBindex = int(BBname.split('_')[1]) - 1
        #---- Check by SPW ID
        for entry in row.findall('name'): SPWname = entry.text
        if 'FULL_RES' not in SPWname: continue
        for entry in row.findall('spectralWindowId'): spw = int(entry.text.split('_')[1])
        #---- Check by Band name
        for entry in row.findall('name'): bandName = 'RB_' + entry.text.split('RB_')[1][0:2]
        for entry in row.findall('numChan'): chNum   = int(entry.text)
        for entry in row.findall('refFreq'): refFreq = float(entry.text)
        for entry in row.findall('totBandwidth'): BW = float(entry.text)
        for entry in row.findall('chanFreqStart'): ch0 = float(entry.text)
        for entry in row.findall('chanWidth'): chWidth = float(entry.text)
        spwDic[spw] = {
            'Band'   : bandName,
            'BB'     : BBindex,
            'chNum'  : chNum,
            'refFreq': refFreq,
            'BW'     : BW,
            'ch0'    : ch0,
            'chWidth': chWidth
        }
    #
    RB_XML = prefix + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Avoid WVR and SQLD
        for entry in row.findall('receiverSideband'): SBname = entry.text
        if 'NOSB' not in SBname: continue
        #-------- Identify BB from SPWID
        for entry in row.findall('frequencyBand'): BandID = int(entry.text.split('_')[-1])
        for entry in row.findall('freqLO'): LO1, LO2, LO3 = float(entry.text.split()[-3]), float(entry.text.split()[-2]),float(entry.text.split()[-1])
        for entry in row.findall('sidebandLO'): SB1, SB2, SB3 = SBsign[entry.text.split()[-3]], SBsign[entry.text.split()[-2]], SBsign[entry.text.split()[-1]]
        refFreq = SB1*SB2*SB3*( SB1*LO1 - SB2*LO2 + SB3*LO3 )
        spwList = [spw for spw in spwDic.keys() if spwDic[spw]['refFreq'] == refFreq]
        if len(spwList) > 0:
            spw = spwList[0]
            spwDic[spw]['LO1'], spwDic[spw]['LO2'], spwDic[spw]['LO3'] = LO1, LO2, LO3
            spwDic[spw]['SB1'], spwDic[spw]['SB2'], spwDic[spw]['SB3'] = SB1, SB2, SB3
    #
    return spwDic
#
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
