import os
import analysisUtils as au
import numpy as np
from interferometry import RADDEG, GetAntName,GetSPWnames,GetPHchavSPWs,GetBPcalScans,GetVisAllBL,gainComplexErr,AllanVarPhase,PhaseDiff,specBunch,bunchVec
from Plotters import plotGain
qa = qatool()
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-t', dest='timeBunch', metavar='timeBunch',
    help='Time average', default='2')
parser.add_option('-j', dest='jumpTH', metavar='jumpTH',
    help='Phase jump threshold [deg]', default='30')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
timeBunch = int(options.timeBunch)
jumpTH    = float(options.jumpTH)
'''
prefix = 'uid___A002_X11fc8d6_X5d9b'
timeBunch = 2
jumpTH = 30.0
'''
msfile = prefix + '.ms'
spwList = GetPHchavSPWs(msfile)
antList = GetAntName(msfile)
def AllanVarPhase30(phase): return AllanVarPhase(phase, 30)
BPlogfile = open(prefix + '-PhaseJump.bplog', 'w')
#-------- Antenna-based gain solutions
for spw_index, spw in enumerate(spwList):
    msmd.open(msfile)
    BPscanList = list(set(msmd.scansforintent('*BANDPASS*')) & set(msmd.scansforspw(spw)))
    msmd.done()
    #-------- CHECK Bandpas Scan
    for scan_index, scan in enumerate(BPscanList):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
        timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
        print('Loading Visibilities: SPW%d Scan%d : %d records' % (spw, scan, timeNum))
        chAvgVis = Xspec[:,0]   # chAvgVis[pol, bl, time]
        if timeBunch > 1:
            chAvgVis = np.array([specBunch(chAvgVis[0], 1, timeBunch), specBunch(chAvgVis[1], 1, timeBunch)])
            timeStamp = bunchVec(timeStamp, timeBunch)
        temp = [np.apply_along_axis(gainComplexErr, 0, chAvgVis[0]), np.apply_along_axis(gainComplexErr, 0, chAvgVis[-1])]
        Gain, Gerr = np.array([temp[0][0], temp[1][0]]), np.array([temp[0][1], temp[1][1]]) # Gain[pol, ant, time], Gerr[pol, ant, time]
        SNR = np.sum(np.mean(abs(Gain), axis=2) / np.median(Gerr.real, axis=2), axis=0)
        #-------- Phase continuity
        phaseAvAnt  = np.array([np.apply_along_axis(AllanVarPhase30, 1, np.angle(Gain[0])), np.apply_along_axis(AllanVarPhase30, 1, np.angle(Gain[1]))])    # [pol,ant]
        phaseDiffAnt = np.array([np.apply_along_axis(PhaseDiff, 1, np.angle(Gain[0])), np.apply_along_axis(PhaseDiff, 1, np.angle(Gain[0]))])   # phaseDiffAnt[pol, ant, time]
        phaseMed = np.median(abs(phaseDiffAnt))
        jump_index = np.where( 900* phaseAvAnt > phaseMed )
        jumpPolList = np.unique(jump_index[0]).tolist()
        jumpAntList = np.unique(jump_index[1]).tolist()
        #-------- candidate antennas of phase jump
        RealJumpAnts = []
        for ant_index in jumpAntList:
            if SNR[ant_index] < 3: continue     # Screening by SNR
            jumpTimingIndex = np.where(abs(phaseDiffAnt[0,ant_index]) > 7.0* phaseMed)[0].tolist() + np.where(abs(phaseDiffAnt[1,ant_index]) > 7.0* phaseMed)[0].tolist()
            jumpTimingIndex = np.sort(np.unique(np.array(jumpTimingIndex)))
            if len(jumpTimingIndex) > 0:
                jumpGap= np.array([
                PhaseDiff(np.array([np.angle(Gain[0,ant_index,jumpTimingIndex[0]]), np.angle(Gain[0,ant_index,jumpTimingIndex[-1]+1])]))[0],
                PhaseDiff(np.array([np.angle(Gain[1,ant_index,jumpTimingIndex[0]]), np.angle(Gain[1,ant_index,jumpTimingIndex[-1]+1])]))[0]])
                if np.max(abs(jumpGap))* RADDEG > jumpTH:
                    RealJumpAnts = RealJumpAnts + [ant_index]
                    text_sd = '%s %d %s %+6.1f %+6.1f %s' % (prefix, spw, antList[ant_index], RADDEG* jumpGap[0], RADDEG*jumpGap[1], au.call_qa_time('%fs' % (timeStamp[jumpTimingIndex[0]]), form='fits', prec=6))
                    print(text_sd); BPlogfile.write(text_sd + '\n')
                #
            #
        #
        #-------- Plot phase jumps
        if len(RealJumpAnts) > 0:
            np.save(prefix + '.Ant.npy', antList[RealJumpAnts])
            np.save('%s-SPW%d.TS.npy' % (prefix, spw), timeStamp)
            np.save('%s-SPW%d.GA.npy' % (prefix, spw), Gain.transpose(1,0,2)[RealJumpAnts])
            np.save('%s-SPW%d.FG.npy' % (prefix, spw), np.ones([len(RealJumpAnts), len(timeStamp)]))
            plotGain(prefix, spw)
            os.system('mv GA_%s-SPW%d.pdf BP_%s-SPW%d.pdf' % (prefix, spw, prefix, spw))
        #
    #
#
BPlogfile.close()
os.system('rm %s*.npy' % (prefix))
PHlogfile = open(prefix + '-PhaseJump.phlog', 'w')
#-------- CHECK Phasecal Scan
for spw_index, spw in enumerate(spwList):
    msmd.open(msfile)
    PHscanList = np.sort(np.array(list(set(msmd.scansforintent('*PHASE*')) & set(msmd.scansforspw(spw))))).tolist()
    msmd.done()
    scanVis, scanTime = [], []
    for scan_index, scan in enumerate(PHscanList):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
        timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]
        print('Loading Visibilities: SPW%d Scan%d : %d records' % (spw, scan, timeNum))
        scanVis = scanVis + [np.mean(Xspec[:,0], axis=2)]   # chAvgVis[pol, bl, time]
        scanTime = scanTime + [np.median(timeStamp)]
    #
    scanVis  = np.array(scanVis).transpose(1,2,0)
    scanTime = np.array(scanTime)
    temp = [np.apply_along_axis(gainComplexErr, 0, scanVis[0]), np.apply_along_axis(gainComplexErr, 0, scanVis[-1])]
    Gain, Gerr = np.array([temp[0][0], temp[1][0]]), np.array([temp[0][1], temp[1][1]])
    SNR = np.sum(np.mean(abs(Gain), axis=2) / np.median(Gerr.real, axis=2), axis=0)
    jumpTHant = [max(jumpTH/RADDEG, 5.0/snr) for snr in SNR.tolist()]   # phase jump threshold, cut by SNR
    phaseDiffAnt = np.array([np.apply_along_axis(PhaseDiff, 1, np.angle(Gain[0])), np.apply_along_axis(PhaseDiff, 1, np.angle(Gain[0]))])   # phaseDiffAnt[pol, ant, time]
    RealJumpAnts = []
    for ant_index, ant in enumerate(antList):
        jumpTimingIndex = np.where(abs(phaseDiffAnt[0,ant_index]) > jumpTHant[ant_index])[0].tolist() + np.where(abs(phaseDiffAnt[1,ant_index]) > jumpTHant[ant_index])[0].tolist()    # phase diff > 5 sigma
        jumpTimingIndex = np.sort(np.unique(np.array(jumpTimingIndex)))
        if len(jumpTimingIndex) > 0: RealJumpAnts = RealJumpAnts + [ant_index]
        for jumpTiming in jumpTimingIndex:
            jumpGap = phaseDiffAnt[:,ant_index][:,jumpTiming]
            text_sd = '%s %d %s %+6.1f %+6.1f %s' % (prefix, spw, ant, RADDEG* jumpGap[0], RADDEG*jumpGap[1], au.call_qa_time('%fs' % (scanTime[jumpTiming]), form='fits', prec=6))
            print(text_sd); PHlogfile.write(text_sd + '\n')
        #
    #-------- Plot phase jumps
    if len(RealJumpAnts) > 0:
        np.save(prefix + '.Ant.npy', antList[RealJumpAnts])
        np.save('%s-SPW%d.TS.npy' % (prefix, spw), scanTime)
        np.save('%s-SPW%d.GA.npy' % (prefix, spw), Gain.transpose(1,0,2)[RealJumpAnts])
        np.save('%s-SPW%d.FG.npy' % (prefix, spw), np.ones([len(RealJumpAnts), len(scanTime)]))
        plotGain(prefix, spw)
        os.system('mv GA_%s-SPW%d.pdf PH_%s-SPW%d.pdf' % (prefix, spw, prefix, spw))
    #
PHlogfile.close()
os.system('rm %s*.npy' % (prefix))
