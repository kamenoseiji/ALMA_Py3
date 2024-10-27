import analysisUtils as au
import numpy as np
from interferometry import RADDEG, GetAntName,GetSPWnames,GetBPchavSPWs,GetBPcalScans,GetVisAllBL,gainComplexErr,AllanVarPhase,PhaseDiff,specBunch,bunchVec
from optparse import OptionParser
qa = qatool()
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-t', dest='timeBunch', metavar='timeBunch',
    help='Time average', default='2')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
timeBunch = int(options.timeBunch)
'''
prefix = 'uid___A002_X11f70c4_Xe2f2'
timeBunch = 1
'''
msfile = prefix + '.ms'
spwList = GetBPchavSPWs(msfile)
antList = GetAntName(msfile)
scanList = GetBPcalScans(msfile)
def AllanVarPhase30(phase): return AllanVarPhase(phase, 30)
logfile = open(prefix + '-PhaseJump.log', 'w')
#-------- Antenna-based gain solutions
for spw_index, spw in enumerate(spwList):
    for scan_index, scan in enumerate(scanList):
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
        jump_index = np.where( 100* phaseAvAnt > phaseMed )
        jumpPolList = np.unique(jump_index[0]).tolist()
        jumpAntList = np.unique(jump_index[1]).tolist()
        for ant_index in jumpAntList:
            if SNR[ant_index] < 3: continue
            jumpTimingIndex = np.where(abs(phaseDiffAnt[0,ant_index]) > 9.0* phaseMed)[0].tolist() + np.where(abs(phaseDiffAnt[1,ant_index]) > 9.0* phaseMed)[0].tolist()
            jumpTimingIndex = np.sort(np.unique(np.array(jumpTimingIndex)))
            if len(jumpTimingIndex) > 0:
                jumpGap= np.array([
                PhaseDiff(np.array([np.angle(Gain[0,ant_index,jumpTimingIndex[0]]), np.angle(Gain[0,ant_index,jumpTimingIndex[-1]+1])]))[0],
                PhaseDiff(np.array([np.angle(Gain[1,ant_index,jumpTimingIndex[0]]), np.angle(Gain[1,ant_index,jumpTimingIndex[-1]+1])]))[0]])
                text_sd = '%s %d %s %+6.1f %+6.1f %s' % (prefix, spw, antList[ant_index], RADDEG* jumpGap[0], RADDEG*jumpGap[1], au.call_qa_time('%fs' % (timeStamp[jumpTimingIndex[0]]), form='fits', prec=6))
                print(text_sd); logfile.write(text_sd + '\n')
            #
        #
    #
logfile.close()

