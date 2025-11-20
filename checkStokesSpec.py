import os
SCR_DIR = os.getenv('HOME') + '/ALMA_Py3/'
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import scipy
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
parser.add_option('-r', dest='BPscan', metavar='BPscan',
    help='Bandpass scan e.g. 3', default='0')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List e.g. 3,6,9,42,67', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA45', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='SPW to process, e.g. 0', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix
refant  = options.refant
antFlag = [ant for ant in options.antFlag.split(',')]
scanList= [int(scan) for scan in options.scanList.split(',')]
spw     = int(options.spw)
BPscan  = int(options.BPscan)
'''
prefix = '2025.1.00003.CSV_V_star_W_Hya_a_02_TM1'
refant = 'DA41'
antFlag = []
spw = 0
BPscan = 0
scanList = [3,61,104]
#scanList = [3]
#----------------------------------------- Procedures
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
#-------- Load Bandpass and XY-phase table 
BPantList, BP_ant = np.load(prefix + '-REF' + refant + '.Ant.npy'), np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, BPscan, spw))
BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
XYspec = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, BPscan, spw))
BP_ant[:,1] *= XYspec
BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
#-------- Load XY phase table
if os.path.isfile('%s-SPW%d-%s.TS.npy' % (prefix, spw, refant)): TS = np.load('%s-SPW%d-%s.TS.npy' % (prefix, spw, refant))
if os.path.isfile('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refant)): XYphase = np.load('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refant))
SP_PH = scipy.interpolate.splrep(TS, XYphase, k=3)
#-------- Load Gain and Flag table
FGantList = np.load('%s.Ant.npy' % (prefix))
UseAntIndexInFG = indexList(UseAntList,FGantList)
unFlaggedAntNum = len(UseAntIndexInFG)
if os.path.isfile('%s-SPW%d.TS.npy' % (prefix, spw)): TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d.GA.npy' % (prefix, spw)): GA = np.load('%s-SPW%d.GA.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d.FG.npy' % (prefix, spw)): FG = np.load('%s-SPW%d.FG.npy' % (prefix, spw))
twiddle = np.exp((1.0j)* scipy.interpolate.splev(TS, SP_PH))
GA[:,1] *= twiddle
UseTimeList = np.where( np.quantile(FG[UseAntIndexInFG], 2.0/unFlaggedAntNum, axis=0) == 1)[0].tolist()
#-------- Load SEFD table
SEFD = np.load('%s-SPW%d.SEFD.npy' % (prefix, spw))
blWeight = 1.0/SEFD[ant0,0]*SEFD[ant1,1]
blWeight = blWeight / np.sum(blWeight)
GA = GA/abs(GA); GA = GA.transpose(2,0,1)*np.sqrt(SEFD)
#-------- Check band and BandPA
spwName = msmd.namesforspws(spw)[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*np.pi/180.0
#-------- Check Stokes Parameters for each source
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
#-------- AZ, EL, PA
azelTime, AntID, AZ, EL = GetAzEl(msfile)
SPW_StokesDic = StokesDicCat
#-------- time-independent spectral setups
chNum, chWid, Freq = GetChNum(msfile, spw); chRange = list(range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum))); Freq = 1.0e-9* bunchVecCH(Freq)
for sourceID in srcDic.keys():
    sourceName = srcDic[sourceID]['Name']
    sourceIDscan = list(set(msmd.scansforfield(sourceID)) & set(scanList))
    scanDic[sourceName] = scanDic[sourceName] + sourceIDscan 
#-------- Load D-term files 
DxSpec, DySpec = [], []
for ant_index, antName in enumerate(antList[antMap]):
    DtermFile = np.load('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, antName))
    DxSpec = DxSpec + [DtermFile[1] + (0.0 + 1.0j)*DtermFile[2]]
    DySpec = DySpec + [DtermFile[3] + (0.0 + 1.0j)*DtermFile[4]]
FreqList = DtermFile[0]
DxSpec = np.array(DxSpec)   # DxSpec[ant, ch]
DySpec = np.array(DySpec)
M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,chNum])).transpose(0,3,1,2)
#-------- For visibilities in each scan
scanVisDic = loadXspecScan(scanVisDic, prefix, spw, bunchNum, TS[UseTimeList])
for scan_index, scan in enumerate(scanVisDic.keys()):
    timeNum = len(scanVisDic[scan]['mjdSec'])
    StokesSpec = np.zeros([timeNum, 4, chNum], dtype=complex)
    scanVisDic[scan]['AZ'], scanVisDic[scan]['EL'] = AzElMatch(scanVisDic[scan]['mjdSec'], azelTime, AntID, refAntID, AZ, EL)
    scanVisDic[scan]['PA'] = AzEl2PA(scanVisDic[scan]['AZ'], scanVisDic[scan]['EL'], ALMA_lat) + BandPA
    scanVisDic[scan]['source'] = [source for source in scanDic.keys() if scan in scanDic[source]][0]
    scanVisDic[scan]['visSpec'] = (CrossPolBL(scanVisDic[scan]['visSpec'][:,:,blMap], blInv).transpose(3, 2, 0, 1) / BP_bl).transpose(3,0,1,2) # [ch,time,bl,pol]
    GAscan = GA[indexList(scanVisDic[scan]['mjdSec'], TS)]
    scanVisDic[scan]['visSpec'] *= (GAscan[:,ant1][:,:,polXindex]*GAscan[:,ant0][:,:,polYindex].conjugate())  # [ch,time,bl,pol]
    PS = InvPAVector(scanVisDic[scan]['PA'], np.ones(len(scanVisDic[scan]['PA'])))
    for time_index in list(range(timeNum)):
        StokesSpec[time_index] = PS[:,:,time_index].dot(np.sum(M* (scanVisDic[scan]['visSpec'][:,time_index].transpose(0,2,1)*blWeight), axis=(2,3)))
    scanVisDic[scan]['StokesSpec'] = np.mean(StokesSpec.real, axis=0)   # [pol, ch]
    scanVisDic[scan]['StokesErr'] = np.std(StokesSpec.imag, axis=0)/np.sqrt(timeNum)
#-------- Store and Plot Stokes spectra
for sourceName in SPW_StokesDic.keys():
    scanLS = scanDic[sourceName]
    StokesSpec, StokesSpecErr, WeightSum = np.zeros([4, chNum]), np.zeros([4, chNum]), 0
    for scan in scanLS:
        scanWeight = len(scanVisDic[scan]['mjdSec'])
        WeightSum += scanWeight
        StokesSpec += (scanWeight* scanVisDic[scan]['StokesSpec'])
    StokesSpec = StokesSpec / WeightSum
    #StokesSpecErr = 1.0/np.sqrt(StokesSpecErr)
    np.save('%s-REF%s-%s-SPW%d.StokesSpec.npy' % (prefix, refant, sourceName, spw), StokesSpec)
    #np.save('%s-REF%s-%s-SPW%d.StokesErr.npy' % (prefix, refant, sourceName, spw), StokesSpecErr)
    np.save('%s-REF%s-%s-SPW%d.Freq.npy' % (prefix, refant, sourceName, spw), Freq)
    #SPW_StokesDic[sourceName] = np.mean(StokesSpec, axis=1).tolist()
msmd.done()
