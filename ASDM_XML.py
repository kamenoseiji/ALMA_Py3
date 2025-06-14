# Module to read XML file in ASDM
#
import xml.etree.ElementTree as ET
import numpy as np
from interferometry import indexList
def CheckCorr( ASDM ):
    Corr_XML = ASDM + '/CorrelatorMode.xml'
    tree = ET.parse(Corr_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #---- Check by Correlat name
        for CorrID in row.findall('correlatorName'): CorrName = CorrID.text
    #
    return CorrName
#
def BandList( ASDM ):
    KEY_ALMA ='ALMA_RB'
    RXList = []
    #-------- Check Receiver bands
    RB_XML = ASDM + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Find Receiver band name
        for RX in row.findall('frequencyBand'): RXname = RX.text
        if KEY_ALMA in RXname:  RXList = RXList + [RXname.replace(KEY_ALMA, 'RB')]
    #
    return list(set(RXList))
#
'''
def BBLOfreq( ASDM ):
    #-------- BB-SPWID connection
    SPW_XML = ASDM + '/' + 'SpectralWindow.xml'
    tree = ET.parse(SPW_XML)
    root = tree.getroot()
    spwList, BBList, bandNameList, spwIDList = [], [], [], []
    for row in root.findall('row'):
        #---- Check by BB name
        for BBID in row.findall('basebandName'): BBname = BBID.text
        if BBname == 'NOBB': continue
        BBindex = int(BBname.split('_')[1]) - 1
        #---- Check by SPW name
        for spwName in row.findall('name'): SPWname = spwName.text
        if 'FULL_RES' not in SPWname: continue
        #---- Check SPWID
        for spwID in row.findall('spectralWindowId'): spw = int(spwID.text.split('_')[1])
        BBList = BBList + [BBindex]
        spwList = spwList + [spw]
        for name in row.findall('name'): bandName = 'RB_' + name.text.split('RB_')[1][0:2]
        bandNameList = bandNameList + [bandName]
    #
    UniqBBList = sorted(set(BBList), key=BBList.index) 
    #for BB in UniqBBList: spwIDList = spwIDList + [spwList[max(indexList(np.array([UniqBBList[BB]]), np.array(BBList)))]]
    #
    #-------- Check LO frequencies 
    for bandName in set(bandNameList):


    #spwNum = len(spwIDList)
    RB_XML = ASDM + '/' + 'Receiver.xml'
    tree = ET.parse(RB_XML)
    root = tree.getroot()
    LO2  = np.zeros(spwNum)
    sideBandSign = [1]* len(UniqBBList)
    for row in root.findall('row'):
        #-------- Avoid WVR and SQLD
        for sideBand in row.findall('receiverSideband'): SBname = sideBand.text
        if 'NOSB' not in SBname: continue
        #-------- Identify BB from SPWID
        for spwID in row.findall('spectralWindowId'):
            spwIDnumber = int(spwID.text.split('_')[1])    
            if spwIDnumber in spwIDList:
                BB_index = spwIDList.index(spwIDnumber)
                for freq in row.findall('freqLO'):
                    freqList = freq.text.split()
                    LO1 = float(freqList[2])
                    LO2[BB_index] = float(freqList[3])
                    for sideBand in row.findall('sidebandLO'):
                        if sideBand.text.split()[2] == 'LSB': sideBandSign[BB_index] = -1
                    #
                #
            #
        #
    #
    return LO1, LO2.tolist(), sideBandSign
#
'''
def SourceList( ASDM ):
    SrcList, PosList,  = [], []
    SRC_XML = ASDM + '/' + 'Source.xml'
    tree = ET.parse(SRC_XML)
    root = tree.getroot()
    for row in root.findall('row'):
        #-------- Find Src Name
        for posCoordinate in row.findall('direction'): PosList = PosList + [[float(posCoordinate.text.strip().split()[2]), float(posCoordinate.text.strip().split()[3])]]
        for srcName in row.findall('sourceName'): SrcList = SrcList + [srcName.text.strip()]
    #
    return SrcList, PosList
#
def PolScan( ASDM ):
    scanList, SrcList, STList, ETList = [], [], [], []
    KEY_ALMA ='CALIBRATE_POLARIZATION'
    #-------- Check Scan XML
    SC_XML = ASDM + '/' + 'Scan.xml'
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
