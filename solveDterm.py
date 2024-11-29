import os
SCR_DIR = os.getenv('HOME') + '/ALMA_Py3/'
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, GetAntName, GetSourceDic, indexList, BANDPA, GetTimerecord, GetPolQuery, BANDFQ, ANT0, ANT1, Ant2BlD, GetAzEl, GetChNum, bunchVec, loadXspecScan, AzElMatch, AzEl2PA, ALMA_lat, CrossPolBL, gainComplex, XXYY2QU, XY2Phase, polariGain, XY2Stokes, XY2PhaseVec, VisMuiti_solveD, InvMullerVector, InvPAVector, get_progressbar_str, RADDEG
import pickle
from Plotters import plotXYP, plotBP, plotSP, lineCmap, plotQUXY, plotXYVis
'''
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List e.g. 3,6,9,42,67', default='')
parser.add_option('-Q', dest='QUmodel', metavar='QUmodel',
    help='Use QU model', action='store_true')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA45', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='SPW to process, e.g. 0', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix
refant  = options.refant
QUmodel = options.QUmodel
antFlag = [ant for ant in options.antFlag.split(',')]
scanList= [int(scan) for scan in options.scanList.split(',')]
spw     = int(options.spw)
'''
prefix = '2018.1.00457.S_Xd54982_X128a'
refant = 'DA44'
QUmodel = True
antFlag = []
spw = 0
BPscan = 0
scanList = [   3, 6, 9, 21, 26, 27, 34, 35, 45]
#----------------------------------------- Procedures
def flagOutLier(value, thresh=5.0):
    return np.where(abs(value - np.median(value)) > thresh* np.std(value))[0].tolist()
#
#if 'antFlag' not in locals():   antFlag = []
polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
scansFile = []
pattern = r'RB_..'
sourceList = []
msfile = prefix + '.ms'
#-------- Check antnna and baseline list
Antenna1, Antenna2 = GetBaselineIndex(msfile, spw, scanList[0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]                    # antList in MS
refAntID = np.where(antList == refant)[0][0]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)  # antNum in MS
BPantList = np.load('%s-REF%s.Ant.npy' % (prefix, refant)).tolist()  # antList in bandpass table
UseAntList= [ant for ant in  BPantList if ant not in antFlag]        # remove flagged antnnas
antMap = indexList(UseAntList, antList) 
UseAntNum = len(UseAntList); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
np.save('%s-SPW%d-%s.Ant.npy' % (prefix, spw, refant), antList[antMap])
#-------- Check source list
srcDic = GetSourceDic(msfile)
srcIDList = list(srcDic.keys())
msmd.open(msfile)
for sourceID in srcIDList:
    if len(set(msmd.scansforfield(sourceID)) & set(scanList)) < 1 : srcDic.pop(sourceID)        # Filter source by scan list
for SSOkey in [key for key in srcDic.keys() if srcDic[key]['RA'] == 0.0]: del srcDic[SSOkey]    # Remove SSO
sourceList = np.unique([srcDic[ID]['Name'] for ID in srcDic.keys()]).tolist(); numSource = len(sourceList)
sourceScan = []
scanDic   = dict(zip(sourceList, [[]]*numSource)) # Scan list index for each source
StokesDic = dict(zip(sourceList, [[]]*numSource)) # Stokes parameters for each source
StokesDicCat = dict(zip(sourceList, [[]]*numSource)) # Stokes parameters in AMAPOLA
#-------- Check scan list
if 'scanList' in locals():
    scanList.sort()
    scanLS = scanList
else:
    scanLS = msmd.scannumbers().tolist()
    scanLS.sort()
#
scanVisDic = dict(zip(scanLS, [[]]* len(scanLS)))
#-------- Time-based flagged data
FGantList = np.load('%s.Ant.npy' % (prefix))
UseAntIndexInFG = indexList(UseAntList,FGantList)
unFlaggedAntNum = len(UseAntIndexInFG)
if os.path.isfile('%s-SPW%d.TS.npy' % (prefix, spw)): TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d.FG.npy' % (prefix, spw)): FG = np.load('%s-SPW%d.FG.npy' % (prefix, spw))
UseTimeList = np.where( np.quantile(FG[UseAntIndexInFG], 3.0/unFlaggedAntNum, axis=0) == 1)[0].tolist()
#-------- Check band and BandPA
spwName = msmd.namesforspws(spw)[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*np.pi/180.0
#-------- Check Stokes Parameters for each source
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
#-------- AZ, EL, PA
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == refAntID )[0].tolist()
if len(azelTime_index) == 0: azelTime_index = np.where(AntID == 0)[0].tolist()
SPW_StokesDic = StokesDicCat
mjdSec = []
#-------- time-independent spectral setups
chNum, chWid, Freq = GetChNum(msfile, spw); chRange = list(range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum))); Freq = 1.0e-9* bunchVecCH(Freq)
for sourceID in srcDic.keys():
    sourceName = srcDic[sourceID]['Name']
    sourceIDscan = list(set(msmd.scansforfield(sourceID)) & set(scanList))
    scanDic[sourceName] = scanDic[sourceName] + sourceIDscan 
    interval, timeStamp = GetTimerecord(msfile, 0, 1, spw, sourceIDscan[0])
    IQU = GetPolQuery(sourceName, timeStamp[0], np.median(Freq), SCR_DIR)
    if len(IQU[0]) == 0: continue
    StokesDicCat[sourceName] = [IQU[0][sourceName], IQU[1][sourceName], IQU[2][sourceName], 0.0]
    print('-- %s I=%.1f p=%.1f%% (est)' % (sourceName, StokesDicCat[sourceName][0], 100.0*np.sqrt(StokesDicCat[sourceName][1]**2 + StokesDicCat[sourceName][2]**2)/StokesDicCat[sourceName][0]))
DxSpec, DySpec = np.zeros([UseAntNum, int(math.ceil(chNum/bunchNum))], dtype=complex), np.zeros([UseAntNum, int(math.ceil(chNum/bunchNum))], dtype=complex)
if 'BPprefix' not in locals():  BPprefix, BPscan = prefix, 0
BPantList, BP_ant = np.load(BPprefix + '-REF' + refant + '.Ant.npy'), np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, refant, BPscan, spw))
BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
if 'XYprefix' not in locals(): XYprefix = prefix
XYspec = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (XYprefix, refant, BPscan, spw))
print('Apply XY phase into Y-pol Bandpass.'); BP_ant[:,1] *= XYspec  # XY phase cal
BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
#-------- For visibilities in each scan
scanVisDic = loadXspecScan(scanVisDic, prefix, spw, bunchNum, TS[UseTimeList])
for scan_index, scan in enumerate(scanVisDic.keys()): mjdSec = mjdSec + scanVisDic[scan]['mjdSec'].tolist()
timeNum, mjdSec = len(mjdSec), np.array(mjdSec)
np.save('%s-SPW%d-%s.TS.npy' % (prefix, spw, refant), mjdSec )
Az, El, PA = np.zeros(timeNum), np.zeros(timeNum), np.zeros(timeNum)
for scan_index, scan in enumerate(scanVisDic.keys()):
    scanVisDic[scan]['flag'] = indexList(TS[UseTimeList], scanVisDic[scan]['mjdSec'])
    Az[scanVisDic[scan]['index']], El[scanVisDic[scan]['index']] = AzElMatch(scanVisDic[scan]['mjdSec'], azelTime, AntID, refAntID, AZ, EL)
    scanVisDic[scan]['EL'] = El[scanVisDic[scan]['index']]
    scanVisDic[scan]['PA'] = AzEl2PA(Az[scanVisDic[scan]['index']], El[scanVisDic[scan]['index']], ALMA_lat) + BandPA
    scanVisDic[scan]['source'] = [source for source in scanDic.keys() if scan in scanDic[source]][0]
    scanVisDic[scan]['visSpec'] = (CrossPolBL(scanVisDic[scan]['visSpec'][:,:,blMap], blInv).transpose(3, 2, 0, 1) / BP_bl).transpose(2,3,1,0) # bandpass cal
    scanVisDic[scan]['visChav'] = np.mean(scanVisDic[scan]['visSpec'][:,chRange], axis=1)
    PA[scanVisDic[scan]['index']] = scanVisDic[scan]['PA']
np.save('%s-SPW%d-%s.Azel.npy' % (prefix, spw, refant), np.array([mjdSec, Az, El, PA]))
caledVis, QCpUS, UCmQS, StokesI = np.ones([4, UseBlNum, timeNum], dtype=complex), np.ones(timeNum), np.ones(timeNum), np.ones(timeNum)
#-------- Parallel-hand phase calibration
print('---- Antenna-based gain solution using tracking antennas')
for scan_index, scan in enumerate(scanVisDic.keys()): caledVis[:,:,scanVisDic[scan]['index']] = scanVisDic[scan]['visChav'] # caledVis[pol, BL, time] : raw chav
Gain = np.array([ np.apply_along_axis(gainComplex, 0, caledVis[0]), np.apply_along_axis(gainComplex, 0, caledVis[-1]) ])
Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))                                   # average in dual parallel polarization
Gain = Gamp* Gain/abs(Gain)                                                     # polarization-averaged gain
caledVis = caledVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())    # caledVis : apply pol-averaged gain correction 
for scan_index, scan in enumerate(scanVisDic.keys()):
    time_index = scanVisDic[scan]['index']
    blWeight = np.mean(caledVis[[0,3]][:,:,time_index].real, axis=(0,2))**2 / (np.var(caledVis[0,:,time_index].real, axis=0) + np.var(caledVis[-1,:,time_index].real, axis=0))
    scanVisDic[scan]['blWeight'] = blWeight / np.sum(blWeight)
    scanVisDic[scan]['scanVis'] = caledVis[:,:,time_index].transpose(0,2,1).dot(scanVisDic[scan]['blWeight'])
    scanVisDic[scan]['Gain'] = Gain[:,:,time_index]
#-------- Coarse estimation of Q and U using XX and YY
print('  -- Solution for Q and U')
if 'QUmodel' not in locals(): QUmodel = False
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    if len(scanLS) < 1 : continue
    if QUmodel:
        QUsol = np.array(SPW_StokesDic[sourceName])[[1,2]]/SPW_StokesDic[sourceName][0]
    else:
        PAList, scanVisXList, scanVisYList = [], [], []
        for scan in scanLS:
            PAList = PAList + scanVisDic[scan]['PA'].tolist()
            scanVisXList = scanVisXList + scanVisDic[scan]['scanVis'][0].tolist()
            scanVisYList = scanVisYList + scanVisDic[scan]['scanVis'][-1].tolist()
        QUsol   = XXYY2QU(np.array(PAList,), np.array([scanVisXList, scanVisYList]))             # XX*, YY* to estimate Q, U
        text_sd = '[XX,YY] %s: Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/np.pi); print(text_sd)
        del PAList, scanVisXList, scanVisYList
    for scan in scanLS:
        time_index = scanVisDic[scan]['index']
        CS, SN = np.cos(2.0* scanVisDic[scan]['PA']), np.sin(2.0* scanVisDic[scan]['PA'])
        scanVisDic[scan]['QCpUS'] = QUsol[0]* CS + QUsol[1]* SN; QCpUS[time_index] = scanVisDic[scan]['QCpUS']
        scanVisDic[scan]['UCmQS'] = QUsol[1]* CS - QUsol[0]* SN; UCmQS[time_index] = scanVisDic[scan]['UCmQS']
#-------- XY phase determination
scanVisXYList, scanVisYXList = [], []
for scan_index, scan in enumerate(scanVisDic.keys()):
    scanVisXYList = scanVisXYList + scanVisDic[scan]['scanVis'][1].tolist()
    scanVisYXList = scanVisYXList + scanVisDic[scan]['scanVis'][2].tolist()
XYphase = XY2Phase(UCmQS, np.array([scanVisXYList, scanVisYXList]))    # XY*, YX* to estimate X-Y phase
XYsign = np.sign(np.cos(XYphase))
text_sd = '  -- Degenerating pi-ambiguity in XY phase : %6.2f [deg] sign = %3.0f' % (XYphase* 180.0 / np.pi, XYsign); print(text_sd)
#-------- Polarized gain adjustment
print('  -- Polarized gain calibration')
GainX, GainY = polariGain(caledVis[0], caledVis[-1], QCpUS)
Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
for scan_index, scan in enumerate(scanVisDic.keys()):
    time_index = scanVisDic[scan]['index']
    caledVis[:,:,time_index] = scanVisDic[scan]['visChav'] / (Gain[polYindex][:,ant0][:,:,time_index]* Gain[polXindex][:,ant1][:,:,time_index].conjugate())
    scanVisDic[scan]['scanVis'] = caledVis[:,:,time_index].transpose(0,2,1).dot(scanVisDic[scan]['blWeight'])
    scanVisDic[scan]['Gain'] = Gain[:,:,time_index]
#-------- Fine estimation of Q and U using XY and YX
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    if len(scanLS) < 1 : continue
    PAList, scanVisXYList, scanVisYXList = [], [], []
    for scan in scanLS:
        PAList = PAList + scanVisDic[scan]['PA'].tolist()
        scanVisXYList = scanVisXYList + scanVisDic[scan]['scanVis'][1].tolist()
        scanVisYXList = scanVisYXList + scanVisDic[scan]['scanVis'][2].tolist()
    Qsol, Usol = XY2Stokes(np.array(PAList), np.array([scanVisXYList, scanVisYXList]))
    text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/np.pi); print(text_sd)
    SPW_StokesDic[sourceName][1], SPW_StokesDic[sourceName][2] = Qsol* SPW_StokesDic[sourceName][0], Usol* SPW_StokesDic[sourceName][0]
    for scan in scanLS:
        CS, SN = np.cos(2.0* scanVisDic[scan]['PA']), np.sin(2.0* scanVisDic[scan]['PA'])
        scanVisDic[scan]['QCpUS'] = Qsol* CS + Usol* SN
        scanVisDic[scan]['UCmQS'] = Usol* CS - Qsol* SN
        QCpUS[scanVisDic[scan]['index']] = scanVisDic[scan]['QCpUS']
        UCmQS[scanVisDic[scan]['index']] = scanVisDic[scan]['UCmQS']
del PAList, scanVisXYList, scanVisYXList
#-------- 2nd polarized gain adjustment
GainX, GainY = polariGain(caledVis[0], caledVis[-1], np.array(QCpUS))
Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
for scan_index, scan in enumerate(scanVisDic.keys()):
    time_index = scanVisDic[scan]['index']
    caledVis[:,:,time_index] = scanVisDic[scan]['visChav'] / (Gain[polYindex][:,ant0][:,:,time_index]* Gain[polXindex][:,ant1][:,:,time_index].conjugate())
    scanVisDic[scan]['Gain'] = Gain[:,:,time_index]
    scanVisDic[scan]['scanVis'] = caledVis[:,:,time_index].transpose(0,2,1).dot(scanVisDic[scan]['blWeight'])
#-------- XY phase correction
XYV = np.mean(caledVis, axis=1)[[1,2]]; XYV = 0.5*(XYV[0] + XYV[1].conjugate())
XYphase, DdotP, DdotM = XY2PhaseVec(mjdSec - np.median(mjdSec), np.mean(caledVis[[1,2]], axis=1), UCmQS, QCpUS, 1000)
for scan_index, scan in enumerate(scanVisDic.keys()):
    time_index = scanVisDic[scan]['index']
    twiddle = np.exp((1.0j)* XYphase[time_index])
    scanVisDic[scan]['XYphase'] = XYphase[time_index]
    scanVisDic[scan]['Gain'][1] *= twiddle
#-------- Display XY cross correlation
ArrayDx, ArrayDy = 0.5* (DdotP - DdotM), 0.5* (DdotP + DdotM)
text_Dx, text_Dy = 'Array Dx = %+.4f %+.4fi' % (ArrayDx.real, ArrayDx.imag), 'Array Dy = %+.4f %+.4fi' % (ArrayDy.real, ArrayDy.imag)
print(text_Dx + ' ' + text_Dy)
pp = PdfPages('XY_%s-REF%s-SPW%d.pdf' % (prefix, refant, spw))
plotXYVis(pp, scanVisDic, DdotP, DdotM)
#-------- SEFD amplitude calibration
FluxList, antGainList = [], []
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    for scan in scanLS:
        antGainList = antGainList + [np.median(abs(scanVisDic[scan]['Gain']), axis=2)]
        FluxList   = FluxList + [StokesDicCat[sourceName][0]]
SEFD = np.sum(np.array(FluxList))/np.sum(np.array(antGainList)**2, axis=0)
for ant_index, ant in enumerate(UseAntList): print('%s SPW%d : SEFD = %.2f Jy / %.2f Jy' % (ant, spw, SEFD[0, ant_index], SEFD[1, ant_index]))
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    fluxList = []
    for scan in scanLS:
        fluxList = fluxList + [np.median(SEFD* np.median(abs(scanVisDic[scan]['Gain']), axis=2)**2)]
    SPW_StokesDic[sourceName] = (np.mean(np.array(fluxList))/SPW_StokesDic[sourceName][0] * np.array(SPW_StokesDic[sourceName])).tolist()
    for scan in scanLS:
        scanVisDic[scan]['Gain'] /= np.sqrt(SPW_StokesDic[sourceName][0])
        CS, SN = np.cos(2.0* scanVisDic[scan]['PA']), np.sin(2.0* scanVisDic[scan]['PA'])
        scanVisDic[scan]['QCpUS'] = SPW_StokesDic[sourceName][1]* CS + SPW_StokesDic[sourceName][2]* SN
        scanVisDic[scan]['UCmQS'] = SPW_StokesDic[sourceName][2]* CS - SPW_StokesDic[sourceName][1]* SN
        StokesI[scanVisDic[scan]['index']] = SPW_StokesDic[sourceName][0]
        QCpUS[scanVisDic[scan]['index']] = scanVisDic[scan]['QCpUS']
        UCmQS[scanVisDic[scan]['index']] = scanVisDic[scan]['UCmQS']
del FluxList, antGainList
for scan_index, scan in enumerate(scanVisDic.keys()):
    caledVis[:,:,scanVisDic[scan]['index']] = scanVisDic[scan]['visChav'] / (scanVisDic[scan]['Gain'][polYindex][:,ant0]* scanVisDic[scan]['Gain'][polXindex][:,ant1].conjugate())
#-------- Antenna-based on-axis D-term (chAvg)
Dx, Dy = VisMuiti_solveD(caledVis, QCpUS, UCmQS, np.repeat(ArrayDx, UseAntNum), np.repeat(ArrayDy, UseAntNum), StokesI)
#-------- D-term-corrected Stokes parameters
Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(UseBlNum, dtype=complex))
print('  -- D-term-corrected visibilities')
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    PAList, srcTimeList = [], []
    for scan in scanLS:
        PAList = PAList + scanVisDic[scan]['PA'].tolist()
        srcTimeList = srcTimeList + scanVisDic[scan]['index']
    srcTimeNum = len(srcTimeList)
    PS = InvPAVector(np.array(PAList), np.ones(srcTimeNum))
    StokesVis = PS.reshape(4, 4*srcTimeNum).dot(Minv.reshape(4, 4*UseBlNum).dot(caledVis[:,:,srcTimeList].reshape(4*UseBlNum, srcTimeNum)).reshape(4*srcTimeNum)) / (srcTimeNum* UseBlNum)
    SPW_StokesDic[sourceName] = StokesVis.real.tolist()
    text_sd = '%s: I= %6.3f+-%6.3f Jy  Q= %6.3f+-%6.3f Jy  U= %6.3f+-%6.3f Jy  EVPA = %6.2f+-%.2f deg' % (sourceName, StokesVis[0].real, abs(StokesVis[0].imag), StokesVis[1].real, abs(StokesVis[1].imag), StokesVis[2].real, abs(StokesVis[2].imag), np.arctan2(StokesVis[2].real, StokesVis[1].real)*90.0/np.pi, 90.0* abs(StokesVis[1].imag* StokesVis[2].real + StokesVis[2].imag* StokesVis[1].real) / (np.pi* abs(StokesVis[1]**2) + abs(StokesVis[2])**2)) 
    print(text_sd)
#-------- get D-term spectra
print('  -- Determining D-term spectra')
for scan_index, scan in enumerate(scanVisDic.keys()):
    scanVisDic[scan]['visSpec'] = (scanVisDic[scan]['visSpec'].transpose(1,0,2,3) / (scanVisDic[scan]['Gain'][polYindex][:,ant0]* scanVisDic[scan]['Gain'][polXindex][:,ant1].conjugate())).transpose(1,0,2,3)
for ch_index in list(range(chNum)):
    GainCaledVis = np.ones([4, UseBlNum, timeNum], dtype=complex)
    for scan_index, scan in enumerate(scanVisDic.keys()):
        time_index = scanVisDic[scan]['index']
        GainCaledVis[:,:,time_index] = scanVisDic[scan]['visSpec'][:,ch_index]
    DxSpec[:,ch_index], DySpec[:,ch_index] = VisMuiti_solveD(GainCaledVis, QCpUS, UCmQS, Dx, Dy, StokesI)
if 'Dsmooth' in locals():
    node_index = list(range(3, chNum/bunchNum, Dsmooth))
    bunchedFreq = bunchVecCH(Freq)
    for ant_index in list(range(UseAntNum)):
        DX_real, DX_imag = scipy.interpolate.splrep(bunchedFreq, DxSpec[ant_index].real, k=3, t=bunchedFreq[node_index]), scipy.interpolate.splrep(bunchedFreq, DxSpec[ant_index].imag, k=3, t=bunchedFreq[node_index])
        DY_real, DY_imag = scipy.interpolate.splrep(bunchedFreq, DySpec[ant_index].real, k=3, t=bunchedFreq[node_index]), scipy.interpolate.splrep(bunchedFreq, DySpec[ant_index].imag, k=3, t=bunchedFreq[node_index])
        DxSpec[ant_index] =scipy.interpolate.splev(bunchedFreq, DX_real) + (0.0 + 1.0j)* scipy.interpolate.splev(bunchedFreq, DX_imag)
        DySpec[ant_index] =scipy.interpolate.splev(bunchedFreq, DY_real) + (0.0 + 1.0j)* scipy.interpolate.splev(bunchedFreq, DY_imag)
del GainCaledVis
for ant_index, ant in enumerate(antList[antMap]):
    DtermFile = np.array([Freq, DxSpec[ant_index].real, DxSpec[ant_index].imag, DySpec[ant_index].real, DySpec[ant_index].imag])
    np.save('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, ant), DtermFile)
#-------- D-term-corrected visibilities (invD dot Vis = PS)
print('  -- Applying D-term spectral correction')
M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([UseBlNum, chNum])).transpose(0,3,1,2) # M[pol, ch, pol, BL]
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    for scan in scanLS:
        PAnum = len(scanVisDic[scan]['PA'])
        scanvisSpec = np.zeros([4, PAnum, chNum])
        scanvisChav = np.zeros([4, PAnum], dtype=complex)
        PS = InvPAVector(scanVisDic[scan]['PA'], np.ones(PAnum))
        for time_index in list(range(PAnum)):
            scanvisSpec[:,time_index] = PS[:,:,time_index].dot(np.sum( M* (scanVisDic[scan]['visSpec'].transpose(3,1,0,2)[time_index]* scanVisDic[scan]['blWeight']), axis=(2,3))).real
            scanvisChav[:,time_index] = np.mean(np.sum( M* (scanVisDic[scan]['visSpec'].transpose(3,1,0,2)[time_index]* scanVisDic[scan]['blWeight']), axis=(2,3)), axis=1)
        scanVisDic[scan]['scanVis'] = np.mean(scanvisSpec, axis=1)
        scanVisDic[scan]['scanVisErr'] = np.std(scanvisSpec, axis=1)
        scanVisDic[scan]['visChav'] = scanvisChav
del scanvisSpec, scanvisChav
pp = PdfPages('%s-SPW%d-%s-QUXY.pdf' % (prefix, spw, refant))
plotQUXY(pp, scanVisDic)
#-------- Save Results
np.save('%s-SPW%d-%s.GA.npy' % (prefix, spw, refant), Gain )
np.save('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refant), XYphase )
XYC = np.zeros([2, timeNum], dtype=complex)
for scan_index, scan in enumerate(scanVisDic.keys()): XYC[:,scanVisDic[scan]['index']] = scanVisDic[scan]['visChav'][[1,2]]
np.save('%s-SPW%d-%s.XYV.npy' % (prefix, spw, refant), XYV )
np.save('%s-SPW%d-%s.XYC.npy' % (prefix, spw, refant), XYC )
#-------- Store and Plot Stokes spectra
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    StokesSpec, StokesSpecErr, WeightSum = np.zeros([4, chNum]), np.zeros([4, chNum]), np.zeros(4)
    for scan in scanLS:
        scanWeight = 1.0 / np.median(scanVisDic[scan]['scanVisErr'], axis=1)**2
        StokesSpecErr += 1.0/scanVisDic[scan]['scanVisErr']**2
        WeightSum += scanWeight
        StokesSpec += (scanVisDic[scan]['scanVis'].T * scanWeight).T
    StokesSpec = (StokesSpec.T / WeightSum).T
    StokesSpecErr = 1.0/np.sqrt(StokesSpecErr)
    np.save('%s-REF%s-%s-SPW%d.StokesSpec.npy' % (prefix, refant, sourceName, spw), StokesSpec)
    np.save('%s-REF%s-%s-SPW%d.StokesErr.npy' % (prefix, refant, sourceName, spw), StokesSpecErr)
    np.save('%s-REF%s-%s-SPW%d.Freq.npy' % (prefix, refant, sourceName, spw), Freq)
    SPW_StokesDic[sourceName] = np.mean(StokesSpec, axis=1).tolist()
#---- Save Stokes parameters of the calibraors
fileDic = open('Stokes.%s-SPW%d.dic' % (prefix, spw), mode='wb')
pickle.dump(SPW_StokesDic, fileDic)
fileDic.close()
#-------- summarize Stokes parameters
StokesTextFile = open('Stokes.%s-SPW%d.txt' % (prefix, spw), mode='w')
text_sd = 'Source       I      Q      U      V      p%    EVPA'
print(text_sd); StokesTextFile.write(text_sd + '\n')
for sourceName in SPW_StokesDic.keys():
    if SPW_StokesDic[sourceName] == []: continue
    text_sd = '%s %6.3f %6.3f %6.3f %6.3f %5.2f %5.2f' % (sourceName, SPW_StokesDic[sourceName][0], SPW_StokesDic[sourceName][1], SPW_StokesDic[sourceName][2], SPW_StokesDic[sourceName][3], 100.0* np.sqrt(SPW_StokesDic[sourceName][1]**2 + SPW_StokesDic[sourceName][2]**2)/SPW_StokesDic[sourceName][0], 90.0* np.arctan2(SPW_StokesDic[sourceName][2], SPW_StokesDic[sourceName][1]) / np.pi)
    print(text_sd); StokesTextFile.write(text_sd + '\n')
#
StokesTextFile.close()
msmd.done()
