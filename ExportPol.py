import os
import math
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-I', dest='INTList', metavar='INTList',
    help='Intents   e.g. POL, BANDPASS', default='POL')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-u', dest='UIDList', metavar='UIDList',
    help='List of UIDS e.g. uid://A002/X11adad7/X19ab0,uid://A002/X11adad7/X1a981,uid://A002/X11adad7/X1b75e', default='')
parser.add_option('-S', dest='SessionName', metavar='SessionName',
    help='Session Name', default='')
parser.add_option('-f', dest='freqRes', metavar='freqRes',
    help='Frequency resolution in MHz', default='31.25')
#
(options, args) = parser.parse_args()
#
INTList = options.INTList.split(',')
INTList = [intent for intent in INTList]
antFlag = options.antFlag.split(',')
antFlag = [ant for ant in antFlag]
Session = options.SessionName
UIDList = options.UIDList.split(',')
UIDList = [uid for uid in UIDList]
freqRes = float(options.freqRes)
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
#-------- Get Bandpass SPWs
def GetBPcalSPWs(msfile):
    msmd.open(msfile)
    bpSPWs  = msmd.spwsforintent("CALIBRATE_POLARIZATION*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_FLUX*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_BANDPASS*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_PHASE*").tolist(); bpSPWs.sort()
    if len(bpSPWs) == 0: bpSPWs  = msmd.spwsforintent("CALIBRATE_DELAY*").tolist(); bpSPWs.sort()
    BPspwList, chNumList, chWidList = [], [], []
    for spw in bpSPWs:
        chNum, chWid, freq = GetChNum(msfile, spw)
        if chNum > 8:               # Filter out WVR and CHAVG spectral windows
            BPspwList = BPspwList + [spw]   # Filter out WVR and CHAVG spectral windows
            chNumList = chNumList + [chNum]
            chWidList = chWidList + [1.0e-6* abs(chWid[0])]
    msmd.close()
    return BPspwList, chNumList, chWidList
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
    if os.path.isdir(prefix + '.ms'): continue
    text_sd = 'python %s/ALMA_Py3/ScanExporterPlus2.py -u %s -i ' % (os.getenv('HOME'), UID)
    for INT in INTList: text_sd = text_sd + INT + ','
    text_sd = text_sd[:-1]
    print(text_sd)
    os.system(text_sd)
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
bpsSPWList, chNumList, chWidList = [], [], []
for file_index in list(range(fileNum)):
    prefix = prefixList[file_index]
    bpSPWs, chNums, chWids = GetBPcalSPWs(prefix + '.ms')
    bpsSPWList = bpsSPWList + [bpSPWs]
    chNumList  = chNumList  + [chNums]
    chWidList  = chWidList  + [chWids]
#
#-------- split and concat for each BB
comvis = []
for file_index, prefix in enumerate(prefixList):
    chNums = chNumList[file_index]
    #---- Check Flag Antenna
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
    chWids = chWidList[file_index]
    chanbin = [1] * len(chWids)
    for spw_index, spw in enumerate(bpsSPWList[file_index]):
        chBunch = max(1, round(freqRes / chWids[spw_index]))
        chBunch = int(chNums[spw_index] / math.ceil(chNums[spw_index]/chBunch))
        if chBunch < chNums[spw_index]: chanbin[spw_index] = int(chNums[spw_index] / int(chNums[spw_index] / chBunch))
        print('%s SPW=%d chWid=%.2f targetRes=%.2f chNum=%d chBin=%d' % (prefix, spw, chWids[spw_index], freqRes, chNums[spw_index], chanbin[spw_index]))
    #---- scan List
    msmd.open(prefix + '.ms')
    scanList = []
    for intent in INTList:
        scanList = scanList + msmd.scansforintent('*%s*' % (intent)).tolist()
    msmd.close()
    scanList.sort()
    #---- split POLcal
    split(prefix+'.ms', outputvis=INTList[0] + prefix + '.ms', scan=",".join(['%d'%(scan) for scan in scanList]), spw=str(bpsSPWList[file_index]).strip('[]'), antenna = removeAnt, width=chanbin, datacolumn='DATA')
    comvis.append(INTList[0] + prefix + '.ms')
#
concat(vis=comvis, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis= Session + '.ms')
listobs(Session +'.ms', spw='', scan='', verbose=True, listfile=Session +'.listobs')
#
