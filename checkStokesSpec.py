import os
SCR_DIR = os.getenv('HOME') + '/ALMA_Py3/'
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import scipy
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, GetAntName, GetSourceDic, indexList, BANDPA, GetTimerecord, GetPolQuery, BANDFQ, ANT0, ANT1, Ant2BlD, GetAzEl, GetChNum, bunchVec, loadXspecScan, AzElMatch, AzEl2PA, ALMA_lat, CrossPolBL, gainComplex, splineInterp, XXYY2QU, XY2Phase, polariGain, XY2Stokes, XY2PhaseVec, VisMuiti_solveD, InvMullerVector, InvPAVector, get_progressbar_str, RADDEG
import pickle
from Plotters import plotXYP, plotBP, plotSP, lineCmap, plotQUXY, plotXYVis
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
prefix = '2025.1.00004.CSV_HR5907_a_02_TM1'
refant = 'DV08'
antFlag = []
spw = 0
BPscan = 0
scanList = [15]
'''
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
#-------- Check Stokes Parameters for each source
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
chNum, chWid, Freq = GetChNum(msfile, spw); chRange = list(range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum))); Freq = 1.0e-9* bunchVecCH(Freq)
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
#-------- Check band and BandPA
spwName = msmd.namesforspws(spw)[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*np.pi/180.0
#-------- AZ, EL, PA
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == refAntID )[0].tolist()
if len(azelTime_index) == 0: azelTime_index = np.where(AntID == 0)[0].tolist()
#-------- Calibration tables
FGantList = np.load('%s.Ant.npy' % (prefix))
UseAntIndexInFG = indexList(UseAntList,FGantList)
unFlaggedAntNum = len(UseAntIndexInFG)
if os.path.isfile('%s-SPW%d.TS.npy' % (prefix, spw)): TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d.FG.npy' % (prefix, spw)): FG = np.load('%s-SPW%d.FG.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d.GA.npy' % (prefix, spw)): GA = np.load('%s-SPW%d.GA.npy' % (prefix, spw)).transpose(1,0,2)
if os.path.isfile('%s-SPW%d.SEFD.npy' % (prefix, spw)): SEFD = np.load('%s-SPW%d.SEFD.npy' % (prefix, spw))
if os.path.isfile('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refant)): XYP = np.load('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refant))
blWeight = 1.0/SEFD[ant0,0]*SEFD[ant1,1]
#-------- Load D-term files 
DxSpec, DySpec = [], []
for ant_index, antName in enumerate(antList[antMap]):
    DtermFile = np.load('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, antName))
    DxSpec = DxSpec + [DtermFile[1] + (0.0 + 1.0j)*DtermFile[2]]
    DySpec = DySpec + [DtermFile[3] + (0.0 + 1.0j)*DtermFile[4]]
FreqList = DtermFile[0]
DxSpec = np.array(DxSpec)   # DxSpec[ant, ch]
DySpec = np.array(DySpec)
M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([UseBlNum,chNum])).transpose(0,3,1,2)
#-------- Check scan list
if 'scanList' in locals():
    scanList.sort()
    scanLS = scanList
else:
    scanLS = msmd.scannumbers().tolist()
    scanLS.sort()
#
for sourceID in srcDic.keys():
    sourceName = srcDic[sourceID]['Name']
    sourceIDscan = list(set(msmd.scansforfield(sourceID)) & set(scanLS))
    sourceIDscan.sort()
    scanDic[sourceName] = scanDic[sourceName] + sourceIDscan
scanVisDic = dict(zip(scanLS, [[]]* len(scanLS)))
scanVisDic = loadXspecScan(scanVisDic, prefix, spw, bunchNum)
#-------- Load Bandpass and XY-phase table 
BPantList, BP_ant = np.load(prefix + '-REF' + refant + '.Ant.npy'), np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, BPscan, spw))
BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
XYspec = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, BPscan, spw))
BP_ant[:,1] *= XYspec
BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
for scan_index, scan in enumerate(scanVisDic.keys()):
    scanVisDic[scan]['AZ'], scanVisDic[scan]['EL'] = AzElMatch(scanVisDic[scan]['mjdSec'], azelTime, AntID, refAntID, AZ, EL)
    scanVisDic[scan]['PA'] = AzEl2PA(scanVisDic[scan]['AZ'], scanVisDic[scan]['EL'], ALMA_lat) + BandPA
    if len(indexList(scanVisDic[scan]['mjdSec'], TS)) > 0 :
        scanAntFlag = FG[:, indexList(scanVisDic[scan]['mjdSec'], TS)]
        scanVisDic[scan]['flag'] = scanAntFlag[ant0]* scanAntFlag[ant1]
    else:
        scanVisDic[scan]['flag'] = np.outer(np.min(FG, axis=1)[ant0]* np.min(FG, axis=1)[ant1], np.ones(len(scanVisDic[scan]['mjdSec'])))
    scanVisDic[scan]['blWeight'] = np.mean(scanVisDic[scan]['flag'], axis=1)* blWeight
    scanVisDic[scan]['blWeight'] /= np.sum(scanVisDic[scan]['blWeight'])
    scanVisDic[scan]['visSpec'] = (CrossPolBL(scanVisDic[scan]['visSpec'][:,:,blMap], blInv).transpose(3, 2, 0, 1) / BP_bl).transpose(2,3,1,0)
    scanVisDic[scan]['visChav'] = np.mean(scanVisDic[scan]['visSpec'][:,chRange], axis=1)* scanVisDic[scan]['flag']
#-------- Gain Calibration
denominator = (1.0 + 0.0j)* abs(GA)* FG
GAPH = np.divide(GA, denominator, out=np.zeros_like(denominator), where=(abs(denominator) > 1.0e-20))
for scan_index, scan in enumerate(scanVisDic.keys()):
    SMwindow = int(len(scanVisDic[scan]['mjdSec'])/2)
    Gain = np.ones([2, antNum, len(scanVisDic[scan]['mjdSec'])], dtype=complex)
    for ant_index in list(range(antNum)):
        Gain[0,ant_index] = splineInterp(TS, GAPH[0,ant_index], scanVisDic[scan]['mjdSec'], SMwindow)/ np.sqrt(SEFD[ant_index,0])
        Gain[1,ant_index] = splineInterp(TS, GAPH[1,ant_index], scanVisDic[scan]['mjdSec'], SMwindow)/ np.sqrt(SEFD[ant_index,1])
    Gain[1] *= splineInterp(XYP[0], np.exp((0.0 - 1.0j)*XYP[1]), scanVisDic[scan]['mjdSec'], SMwindow)
    scanVisDic[scan]['Gain'] = Gain
    denominator = Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate()
    scanVisDic[scan]['visSpec'] = np.divide(scanVisDic[scan]['visSpec'].transpose(1,0,2,3), denominator, out=np.zeros_like(scanVisDic[scan]['visSpec'].transpose(1,0,2,3)), where=(abs(denominator) > 1.0e-20)).transpose(1,0,2,3)
    #-------- D-term correction
    PAnum = len(scanVisDic[scan]['PA'])
    scanvisSpec = np.zeros([4, PAnum, chNum])
    scanvisChav = np.zeros([4, PAnum], dtype=complex)
    PS = InvPAVector(scanVisDic[scan]['PA'], np.ones(PAnum))
    for time_index in list(range(PAnum)):
        scanvisSpec[:,time_index] = PS[:,:,time_index].dot(np.sum( M* (scanVisDic[scan]['visSpec'].transpose(3,1,0,2)[time_index]* scanVisDic[scan]['blWeight']), axis=(2,3))).real
        scanvisChav[:,time_index] = np.mean(np.sum( M* (scanVisDic[scan]['visSpec'].transpose(3,1,0,2)[time_index]* scanVisDic[scan]['blWeight']), axis=(2,3)), axis=1)
    scanVisDic[scan]['StokesSpec'] = np.mean(scanvisSpec.real, axis=1)   # [pol, ch]
    scanVisDic[scan]['StokesErr'] = np.std(scanvisSpec.real, axis=1)/np.sqrt(PAnum)
#
StokesTextFile = open('StokesScan.%s-SPW%d.txt' % (prefix, spw), mode='w')
text_sd = '  Source     Scan            UTC       Frequency    Stokes I [Jy]  Stokes Q [Jy]    Stokes U [Jy]    Stokes V [Jy]    p [%]             EVPA [deg]'
print(text_sd); StokesTextFile.write(text_sd + '\n')
np.save('%s-REF%s-%s-SPW%d.Freq.npy' % (prefix, refant, sourceName, spw), Freq)
for sourceName in StokesDicCat.keys():
    scanLS = scanDic[sourceName]
    for scan in scanLS:
        timeLabel = qa.time('%fs' % np.median(scanVisDic[scan]['mjdSec']), form='ymd')[0]
        StokesFlux = np.mean(scanVisDic[scan]['StokesSpec'][:,chRange], axis=1)
        StokesErr =  np.median(scanVisDic[scan]['StokesErr'][:,chRange], axis=1)/np.sqrt(len(chRange))
        np.save('%s-REF%s-%s-SPW%d-Scan%d.StokesSpec.npy' % (prefix, refant, sourceName, spw, scan), scanVisDic[scan]['StokesSpec'])
        np.save('%s-REF%s-%s-SPW%d-Scan%d.StokesErr.npy' % (prefix, refant, sourceName, spw, scan), scanVisDic[scan]['StokesErr'])
        text_sd = '%12s %4d %s %6.2f GHz   %6.3f (%.3f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %+7.2f (%.2f) ' % (sourceName, scan, timeLabel, np.median(Freq), StokesFlux[0], StokesErr[0], StokesFlux[1], StokesErr[1], StokesFlux[2], StokesErr[2], StokesFlux[3], StokesErr[3], 100.0*np.sqrt(StokesFlux[1]**2 + StokesFlux[2]**2)/StokesFlux[0], 100.0*np.sqrt(StokesFlux[1]**2 + StokesFlux[2]**2)*np.sqrt( (StokesErr[0]/StokesFlux[0])**2 + (StokesErr[1]/StokesFlux[1])**2 + (StokesErr[2]/StokesFlux[2])**2 ) / StokesFlux[0]  , 0.5* RADDEG* np.arctan2(StokesFlux[2], StokesFlux[1]), 0.5* RADDEG* abs(StokesErr[1]* StokesFlux[2] + StokesErr[2]* StokesFlux[1]) / (abs(StokesFlux[1])**2 + abs(StokesFlux[2])**2))
        print(text_sd); StokesTextFile.write(text_sd + '\n')
StokesTextFile.close()
msmd.done()
