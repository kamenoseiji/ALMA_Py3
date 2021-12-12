#---- Script for Band-3 Astroholograpy Data
import sys
from scipy import stats
import matplotlib.pyplot as plt
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#-------- Initial Settings
if 'SNR_THRESH' not in locals(): SNR_THRESH = 3.0
msfile = wd + prefix + '.ms'; msmd.open(msfile)
antList = GetAntName(msfile)
antNum = len(antList)
blNum = int(antNum* (antNum - 1) / 2)
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
#timeStamp, UVW = GetUVW(msfile, spw, msmd.scansforspw(spw)[0])
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
GainAP0, GainAP1, timeList, SNRList, flagList, fieldList = [], [], [], [], [], []
scan_index = 0
for scan in scanList:
    field_names = msmd.fieldsforscan(scan, True); fieldList = fieldList + [field_names[0]]
    print('Loading Visibilities: Scan %d : %s' % (scan, field_names[0]))
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
    if polNum == 4: polIndex = [0, 3]
    if polNum == 2: polIndex = [0, 1]
    if polNum == 1: polIndex = [0]
    polNum = len(polIndex)
    if chNum == 1:  chRange = [0]
    if 'chRange' not in locals(): chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
    print('Baseline and polarization mapping... ')
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
    for time_index in list(range(timeNum)):
        #------ Progress bar
        progress = (time_index + 1.0) / timeNum
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
        #------
        gainFlag = np.ones(UseAntNum)
        tempGain, tempErr = gainComplexErr(chAvgVis[0, :, time_index]); GainAP0 = GainAP0 + [tempGain]; tempSNR = abs(tempGain) / tempErr
        SNRList = SNRList + [tempSNR]; gainFlag[np.where(tempSNR  < SNR_THRESH)[0]] = 0.0
        #------
        if polNum == 2:
            tempGain, tempErr = gainComplexErr(chAvgVis[1, :, time_index]); GainAP1 = GainAP1 + [tempGain]; tempSNR = abs(tempGain) / tempErr
            SNRList = SNRList + [tempSNR]; gainFlag[np.where(tempSNR  < SNR_THRESH)[0]] = 0.0
        #
        flagList = flagList + [gainFlag]
    #
    timeList.extend(timeStamp.tolist())
    scan_index += 1
#
msmd.done(); msmd.close()
timeNum = len(timeList)
antFG  = np.array(flagList).T                          # [ant, time]
antSNR = np.array(SNRList).reshape(timeNum, polNum, UseAntNum).transpose(2,1,0)  # [ant, pol, time]
if polNum == 2:
    Gain = np.array([GainAP0, GainAP1]).transpose(2,0,1)    # [ant, pol, time]
else:
    Gain = np.array([GainAP0]).transpose(2,0,1)    # [ant, pol, time]
np.save(prefix + '.Ant.npy', antList[antMap]) 
np.save(prefix + '.Field.npy', np.array(fieldList))
np.save('%s-SPW%d.TS.npy' % (prefix, spw), np.array(timeList)) 
np.save('%s-SPW%d.GA.npy' % (prefix, spw), Gain) 
np.save('%s-SPW%d.FG.npy' % (prefix, spw), antFG) 
plotGain(prefix, spw)
