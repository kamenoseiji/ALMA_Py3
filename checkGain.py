#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#-------- Initial Settings
if 'SNR_THRESH' not in locals(): SNR_THRESH = 0.0
if 'antFlag' not in locals(): antFlag = []
msfile = wd + prefix + '.ms'
Antenna1, Antenna2 = GetBaselineIndex(msfile, spw, scanList[0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
msmd.open(msfile)
spwName = msmd.namesforspws(spw)[0]
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
timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList[UseAnt])[0]
else: refantID = bestRefant(uvDist)
print('  Use ' + antList[UseAnt[refantID]] + ' as the refant.')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print('  %d baselines are inverted.' % (len(np.where( blInv )[0])))
#-------- Bandpass Table
if 'BPprefix' in locals():
    BPfileName = '%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, antList[UseAnt[refantID]], BPscan, spw)
    print('---Loading bandpass table : ' + BPfileName)
    BP_ant = np.load(BPfileName)
#
#-------- Loop for Scan
GainList, timeList, flagList, fieldList = [], [], [], []
for scan_index in list(range(len(scanList))):
    scan = scanList[scan_index]
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
    if 'timeBunch' in locals():
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
plotGain(prefix, spw)
