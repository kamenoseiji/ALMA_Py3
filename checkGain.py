#---- Script for Band-3 Astroholograpy Data
import numpy as np
import matplotlib.pyplot as plt
from interferometry import GetBaselineIndex, CrossCorrAntList, GetAntName, GetUVW, GetVisAllBL, indexList, Ant2Bl, Ant2BlD, ANT0, ANT1, bestRefant, specBunch, bunchVec, ParaPolBL, gainComplexErr
from Plotters import plotGain
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-b', dest='BPprefix', metavar='BPprefix',
    help='Bandpass uid and scan to apply e.g. uid___A002_X10dadb6_X18e6,3', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan ID  e.g. 3,5,7', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='SPW e.g. 26', default='')
parser.add_option('-t', dest='timeBunch', metavar='timeBunch',
    help='Time average', default='1')
parser.add_option('-T', dest='threshold', metavar='threshold', type="float",
    help='SNR threshold to flag', default=0.0)
parser.add_option('-P', dest='PLOTPDF', metavar='PLOTPDF',
    help='Plot PDF', action="store_true")
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
antFlag = [ant for ant in options.antFlag.split(',')]
timeBunch = int(options.timeBunch)
scanList= [int(scan) for scan in options.scanList.split(',')]
spw = int(options.spw)
SNR_THRESH = options.threshold
PLOTPDF = options.PLOTPDF
#-------- Initial Settings
#if 'SNR_THRESH' not in locals(): SNR_THRESH = 0.0
#if 'antFlag' not in locals(): antFlag = []
#if 'msfile' not in locals(): msfile = wd + prefix + '.ms'
msfile = prefix + '.ms'
Antenna1, Antenna2 = GetBaselineIndex(msfile, spw, scanList[0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
msmd.open(msfile)
spwName = msmd.namesforspws(spw)[0]
if spwName != 'none':
    BandName = re.findall(r'RB_..', spwName)[0]; BandID = int(BandName[3:5])
#-------- Array Configuration
print('---Checking array configuration')
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print(text_sd)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
if options.refant not in antList[UseAnt]: options.refant = ''
if options.refant == '':
    timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
else:
    refantID = indexList(np.array([options.refant]), antList[UseAnt])[0]
print('  Use ' + antList[UseAnt[refantID]] + ' as the refant.')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print('  %d baselines are inverted.' % (len(np.where( blInv )[0])))
#-------- Bandpass Table
if options.BPprefix != '':
    BPfileName = '%s-REF%s-SC%d-SPW%d-BPant.npy' % (options.BPprefix.split(',')[0], antList[UseAnt[refantID]], int(options.BPprefix.split(',')[1]), spw)
    print('---Loading bandpass table : ' + BPfileName)
    BP_ant = np.load(BPfileName)
#
#-------- Loop for Scan
GainList, timeList, flagList, fieldList = [], [], [], []
for scan_index, scan in enumerate(scanList):
    field_names = msmd.fieldsforscan(scan, True); fieldList = fieldList + [field_names[0]]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
    print('Loading Visibilities: Scan %d : %s : %d records' % (scan, field_names[0], timeNum))
    if polNum == 4: polIndex = [0, 3]
    if polNum == 2: polIndex = [0, 1]
    if polNum == 1: polIndex = [0]
    parapolNum = len(polIndex)
    if chNum == 1:
        chRange = [0]
    else:
        chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
    #
    tempSpec = ParaPolBL(Xspec[polIndex][:,:,blMap], blInv).transpose(3,2,0,1)  # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    if 'BP_ant' in locals():
        print('Applying bandpass calibration...')
        BPCaledXspec = (tempSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    if timeBunch > 1:
        chAvgVis = np.array([specBunch(chAvgVis[0], 1, timeBunch), specBunch(chAvgVis[1], 1, timeBunch)])
        timeNum = chAvgVis.shape[2]
        timeStamp = bunchVec(timeStamp, timeBunch)
    Gain, antSNR, gainFlag = np.ones([UseAntNum, parapolNum, len(timeStamp)], dtype=complex), np.zeros([UseAntNum, parapolNum, len(timeStamp)]), np.ones([UseAntNum,len(timeStamp)])
    for pol_index in list(range(parapolNum)):
        Gain[:, pol_index], tempErr = np.apply_along_axis(gainComplexErr, 0, chAvgVis[pol_index])
        antSNR[:, pol_index] = abs(Gain[:, pol_index])**2 / abs(tempErr)**2
        gainFlag = gainFlag* ((np.sign( antSNR[:, pol_index] - SNR_THRESH ) + 1.0)/2)
    #
    timeList = timeList + timeStamp.tolist()
    GainList = GainList + [Gain]
    flagList = flagList + [gainFlag]
#
msmd.done(); msmd.close()
timeNum = len(timeList)
GAarray, FGarray = np.ones([UseAntNum, parapolNum, timeNum], dtype=complex), np.ones([UseAntNum, timeNum])
time_index = 0
for scan_index in list(range(len(scanList))):
    timeNum = GainList[scan_index].shape[2]
    GAarray[:,:,time_index:(time_index + timeNum)] = GainList[scan_index]
    FGarray[:,time_index:(time_index   + timeNum)] = flagList[scan_index]
    time_index += timeNum
#
newFlagIndex = np.where( np.median(FGarray, axis=1) == 0)[0].tolist()
newAntFlag = antList[newFlagIndex].tolist()
if len(newAntFlag) > 0: print('Flagged by SNR : %s' % (str(newAntFlag)))
np.save(prefix + '.Ant.npy', antList[antMap]) 
np.save(prefix + '.Field.npy', np.array(fieldList))
np.save('%s-SPW%d.TS.npy' % (prefix, spw), np.array(timeList)) 
np.save('%s-SPW%d.GA.npy' % (prefix, spw), GAarray)
np.save('%s-SPW%d.FG.npy' % (prefix, spw), FGarray) 
if PLOTPDF: plotGain(prefix, spw)
