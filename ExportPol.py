INTList = ['POLARIZATION']
CATList = ['POL_']
#----
def indexList( refArray, keyWord ):     # Compare two arrays and return matched index
    IL = []
    for index in range(len(refArray)):
        if keyWord in refArray[index]: IL = IL + [index]
    return IL
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
    bpSPWs  = msmd.spwsforintent("CALIBRATE_PHASE*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_POLARIZATION*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_FLUX*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
    BPspwList = []
    for spw in bpSPWs:
        chNum, chWid, freq = GetChNum(msfile, spw)
        if chNum > 4:   BPspwList = BPspwList + [spw]   # Filter out WVR and CHAVG spectral windows
    msmd.close()
    return BPspwList
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
    for index in list(range(len(INTList))):
        text_sd = 'python /users/skameno/bin/ScanExporterPlus2.py -u %s -i %s' % (UID, INTList[index])
        print(text_sd)
        os.system(text_sd)
    #
#
print(prefixList)
fileNum = len(prefixList)
#-------- asdm2MS and listobs
for prefix in prefixList:
    importasdm(prefix)
    listobs(prefix+'.ms', spw='', scan='', verbose=True, listfile=prefix+'.listobs')
#
#-------- Check SPW list
bpsSPWList = []
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    bpsSPWList = bpsSPWList + [GetBPcalSPWs(prefix + '.ms')]
#
#-------- split and concat for each BB
comvis = []
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    #---- split POLcal
    split(prefix+'.ms', outputvis=CATList[0] + prefix + '.ms', spw=str(bpsSPWList[0]).strip('[]'), datacolumn='DATA')
    comvis.append(CATList[0] + prefix + '.ms')
    #
#
concat(vis=comvis, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis= Session + '.ms')
listobs(Session +'.ms', spw='', scan='', verbose=True, listfile=Session +'.listobs')
#
#-------- cleanup temporary files
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    os.system('rm -rf *' + prefix+'*.ms*')
#
