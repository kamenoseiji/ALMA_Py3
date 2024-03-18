import os
INTList = ['POL']
CATList = ['POL_']
#----
def indexList( refArray, keyWord ):     # Compare two arrays and return matched index
    IL = []
    for index in range(len(refArray)):
        if keyWord in refArray[index]: IL = IL + [index]
    return IL
#
def GetAntName(msfile):
    tb.open(msfile+'/'+'ANTENNA')
    namelist = tb.getcol("NAME")
    tb.close()
    return namelist
#
def GetAtmSPWs(msfile):
    msmd.open(msfile)
    atmSPWs = list( (set(msmd.tdmspws()) | set(msmd.fdmspws())) & set(msmd.spwsforintent("CALIBRATE_ATMOSPHERE*")) ); atmSPWs.sort()
    if len(atmSPWs) == 0:
        atmSPWList = msmd.spwsforintent("CALIBRATE_ATMOSPHERE*").tolist()
        tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
        for spwID in atmSPWList:
            if tb.getcell("NUM_CHAN", spwID) > 60: atmSPWs = atmSPWs + [spwID]
        #
        tb.close()
    #
    msmd.close()
    return atmSPWs
#
#-------- Get Bandpass SPWs
def GetBPcalSPWs(msfile):
    msmd.open(msfile)
    bpSPWs  = msmd.spwsforintent("CALIBRATE_POLARIZATION*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_FLUX*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_PHASE*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
    BPspwList, chNumList = [], []
    for spw in bpSPWs:
        chNum, chWid, freq = GetChNum(msfile, spw)
        if chNum > 4:               # Filter out WVR and CHAVG spectral windows
            BPspwList = BPspwList + [spw]   # Filter out WVR and CHAVG spectral windows
            chNumList = chNumList + [chNum]
    msmd.close()
    return BPspwList, chNumList
#
def GetChNum(msfile, spwID):
    tb.open(msfile + '/' + 'SPECTRAL_WINDOW')
    chNum = tb.getcell("NUM_CHAN", spwID)
    chWid = tb.getcell("CHAN_WIDTH", spwID)
    freq  = tb.getcell("CHAN_FREQ", spwID)
    tb.close()
    return chNum, chWid, freq
#
#-------- Export ASDM
prefixList = []
for UID in UIDList:
    prefix = UID.replace("/", "_").replace(":","_").replace(" ","")
    prefixList = prefixList + [prefix]
    if os.path.isdir(prefix): continue
    for index, INT in enumerate(INTList):
        text_sd = 'python /users/skameno/bin/ScanExporterPlus2.py -u %s -i %s' % (UID, INT)
        print(text_sd)
        os.system(text_sd)
    #
#
print(prefixList)
fileNum = len(prefixList)
#-------- asdm2MS and listobs
for prefix in prefixList:
    if os.path.isdir(prefix + '.ms'): continue
    importasdm(prefix)
    listobs(prefix+'.ms', spw='', scan='', verbose=True, listfile=prefix+'.listobs')
#
#-------- Check SPW list
bpsSPWList, chNumList = [], []
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    bpSPWs, chNums = GetBPcalSPWs(prefix + '.ms')
    bpsSPWList = bpsSPWList + [bpSPWs]
    chNumList  = chNumList  + [chNums]
#
#-------- split and concat for each BB
comvis = []
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    spwNum = len(bpsSPWList[file_index])
    #---- Check Flag Antenna
    if 'antFlag' not in locals():    antFlag = []
    removeAnt = ''
    if len(antFlag) > 0:
        removeAnt = '!'
        antList = GetAntName(prefix + '.ms').tolist()
        for antName in antFlag:
            if antName not in antList: continue
            removeAnt = removeAnt + antName + ','
        removeAnt = removeAnt.rstrip(',')
    #
    if len(removeAnt) < 2:  removeAnt = ''
    #---- Channel binning
    chanbin = [1] * spwNum
    if 'chBunch' in locals():
        for spw_index in list(range(spwNum)):
            if chBunch < chNumList[file_index][spw_index]: 
                chanbin[spw_index] = int(chNumList[file_index][spw_index] / int(chNumList[file_index][spw_index] / chBunch))
            #
        #
    #---- scan List
    msmd.open(prefix + '.ms')
    scanList = msmd.scansforintent('*POLARIZATION*').tolist()
    msmd.close()
    #---- split POLcal
    split(prefix+'.ms', outputvis=CATList[0] + prefix + '.ms', scan=",".join(['%d'%(scan) for scan in scanList]), spw=str(bpsSPWList[file_index]).strip('[]'), antenna = removeAnt, width=chanbin, datacolumn='DATA')
    comvis.append(CATList[0] + prefix + '.ms')
    #
#
concat(vis=comvis, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis= Session + '.ms')
listobs(Session +'.ms', spw='', scan='', verbose=True, listfile=Session +'.listobs')
#
#-------- cleanup temporary files
'''
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    os.system('rm -rf *' + prefix+'*.ms*')
#
'''
