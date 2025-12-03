import os
import math
from casatools import quanta as qatool
import analysisUtils as au
import numpy as np
import scipy
import datetime
from interferometry import indexList, GetChNum, bunchVec, delay_search, Bl2Ant, Ant2Bl, RADDEG
from Grid import tauSMTH, SSOCatalog, lmStokes
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
qa = qatool()
polColor = ['b', 'g', 'm', 'y']
polName = ['X', 'Y', 'XY', 'YX']
#
#-------- Set Color Map
lineCmap = plt.get_cmap('Set1')
#-------- Plot optical depth
def plotTauSpec(prefix, spwList, freqList, Tau0spec):
    figTauSP = plt.figure(0, figsize = (11,8))
    figTauSP.suptitle(prefix + ' Zenith Opacity')
    figTauSP.text(0.45, 0.05, 'Frequency [GHz]')
    figTauSP.text(0.03, 0.45, 'Optical Depth', rotation=90)
    spwNum = len(spwList)
    plotMax = 0.01
    for spw_index in range(spwNum): plotMax = max(plotMax, np.max(Tau0spec[spw_index]))
    for spw_index in range(spwNum):
        chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
        TauPL = figTauSP.add_subplot(1, spwNum, spw_index + 1 )
        TauPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, 1.05* plotMax])
        TauPL.tick_params(axis='both', labelsize=6)
        TauPL.step(freqList[spw_index][chRange], Tau0spec[spw_index][chRange], color='k', where='mid')
        text_sd = 'SPW = %d' % (spwList[spw_index])
        TauPL.text(np.min(freqList[spw_index]), 1.01* plotMax, text_sd, fontsize='8')
    #
    figTauSP.savefig('TAUS_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot Tau fit
def plotTauFit(prefix, antList, spwList, secZ, tempAtm, Tau0, TantN, TskyList, scanFlag):
    antNum = len(antList)
    spwNum = len(spwList)
    airmass = np.arange( 1.0, 1.25*np.max(secZ), 0.01)
    figTauFit = plt.figure(0, figsize = (11,8))
    figTauFit.suptitle(prefix + ' Optical Depth')
    figTauFit.text(0.45, 0.05, 'Airmass')
    figTauFit.text(0.03, 0.45, 'Sky Temperature [K]', rotation=90)
    for spw_index in range(spwNum):
        chAvgTsky = np.median(TskyList[spw_index].transpose(2, 1, 0) - TantN[spw_index], axis=2).T  # chAvgTsky[ant, scan]
        chAvgTau0 = np.median(Tau0[spw_index])
        plotMax = 1.2 * np.max(chAvgTsky)
        TskyPL = figTauFit.add_subplot(1, spwNum, spw_index + 1 )
        TskyPL.axis([1.0, 2.5, 0.0, plotMax])
        TskyPL.plot( airmass, au.Tcmb* np.exp(-chAvgTau0* airmass) + tempAtm* (1.0 - np.exp(-chAvgTau0* airmass)), '-', color='k', alpha=0.5)
        for ant_index in range(antNum):
            rgb = lineCmap(float(ant_index) / antNum )
            plotTsky = chAvgTsky[ant_index]
            TskyPL.scatter( secZ, plotTsky, s=10.0* scanFlag[spw_index, ant_index], color=rgb, label = antList[ant_index])
        #
        text_sd = 'Tau(zenith)=%6.4f' % (chAvgTau0)
        TskyPL.text(1.01, 0.95* plotMax, text_sd, fontsize='9')
        TskyPL.set_title('SPW %d' % (spwList[spw_index]))
    #
    TskyPL.legend(loc = 'lower right', prop={'size' :7}, numpoints = 1)
    figTauFit.savefig('TAUF_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot plotTauEOn
def plotTauEOn(prefix, bandName, spw, onTime, onData, TauEOn):
    DT = []
    for mjdSec in onTime.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle('%s %s SPW %d CHAV data and zenith optical depth' % (prefix, bandName, spw))
    figTau.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    CHAVPL = figTau.add_subplot(2, 1, 1)
    CHAVPL.plot(DT, onData, 'k.')
    CHAVPL.set_ylabel('Channel-Averaded Power')
    TAU0PL = figTau.add_subplot(2, 1, 2)
    TAU0PL.plot(DT, TauEOn, 'b.')
    TAU0PL.set_ylabel('Estimated Zenith Optical Depth')
    figTau.savefig('TauEOn_%s-%s-SPW%d.pdf' % (prefix, bandName, spw))
    plt.close('all')
    return
#-------- Plot Tau0-Excess
def plotTau0E(prefix, atmTime, spwList, Tau0, Tau0Excess, scanFlag):
    spwNum = len(spwList)
    figTauE = plt.figure(0, figsize = (11,8))
    figTauE.suptitle(prefix + ' Zenith Optical Depth')
    DT, DTSpl = [], []
    for mjdSec in atmTime.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    mjdSpl = np.arange(atmTime[0], atmTime[-1], 1)
    for mjdSec in mjdSpl.tolist(): DTSpl.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    figTauE.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    figTauE.text(0.03, 0.45, 'Zenith Optical Depth', rotation=90)
    TauEPL = figTauE.add_subplot(1, 1, 1)
    for spw_index, spw in enumerate(spwList):
        Tau0E = np.median(Tau0[spw_index]) + Tau0Excess[spw_index]
        SP = tauSMTH( atmTime-atmTime[0], Tau0E )
        Tau0ESpl = scipy.interpolate.splev(mjdSpl-atmTime[0], SP)
        TauEPL.plot( DTSpl, Tau0ESpl, '-', color='gray')
        TauEPL.scatter( DT, Tau0E, s=10.0* scanFlag[spw_index], label='SPW %d' % (spw))
    TauEPL.tick_params(axis='x', labelsize=6)
    TauEPL.legend(loc='best', prop={'size' :10}, numpoints = 1)
    figTauE.savefig('TAUE_' + prefix + '.pdf')
    plt.close('all')
    return
#-------- Plot Tsys and Trx Spectrum
def plotTsysDic(prefix, TsysDic):
    pp = PdfPages('TSYS_' + prefix + '.pdf')
    #-------- Plots for Tsys spectra
    scanList = list(TsysDic.keys()); scanNum = len(scanList)
    antList  = TsysDic[scanList[0]]['antList']; antNum = len(antList)
    #-------- Prepare Plots
    figAnt = plt.figure(figsize = (8, 11))
    figAnt.suptitle(prefix)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.42, 'Tsys (X:blue, Y:green) and Trec (X:magenta, Y:yellow) [K]', rotation=90)
    #-------- Page for each antenna
    for ant_index, ant in enumerate(antList):
        if ant_index > 0:
            for PL in TsysPL: figAnt.delaxes(PL)
        TsysPL, TsysMax = [], []
        #-------- To determine scaling max
        for scan_index, scan in enumerate(scanList):
            spwNum = len(TsysDic[scan]['spwList'])
            for spw_index, spw in enumerate(TsysDic[scan]['spwList']):
                index = spwNum* ant_index + spw_index
                TsysMax = TsysMax + [np.percentile(TsysDic[scan]['Tsys'][:,:,index], 75)]
        plotMax = 1.5* np.max(np.array(TsysMax))
        for scan_index, scan in enumerate(scanList):
            spwNum = len(TsysDic[scan]['spwList'])
            for spw_index, spw in enumerate(TsysDic[scan]['spwList']):
                index = spwNum* ant_index + spw_index
                TRX, TSYS, FREQ = TsysDic[scan]['Trx'][:,:,index], TsysDic[scan]['Tsys'][:,:,index], 1.0e-9*TsysDic[scan]['freq'][spw_index]
                polNum = TRX.shape[0]
                currentPL = figAnt.add_subplot(scanNum, spwNum, spwNum* scan_index + spw_index + 1 )
                TsysPL = TsysPL + [currentPL]
                timeLabel = au.call_qa_time('%fs' % (TsysDic[scan]['mjdSec']), form='fits')
                for pol_index in range(polNum):
                    currentPL.step(FREQ, TSYS[pol_index], where='mid', color=polColor[pol_index], label = 'Tsys Pol '+ polName[pol_index])
                    currentPL.step(FREQ, TRX[pol_index],  color=polColor[pol_index+2], where='mid', label = 'Trec Pol ' + polName[pol_index])
                currentPL.axis([np.min(FREQ), np.max(FREQ), 0.0, plotMax])
                currentPL.tick_params(axis='both', labelsize=3, direction='in')
                currentPL.text(0.25*np.min(FREQ) + 0.75*np.max(FREQ), 0.8* plotMax, 'SPW %d' % (spw), fontsize='6')
                if scan_index == 0 and spw_index == 0: currentPL.text(1.2* np.min(FREQ) - 0.2* np.max(FREQ), (0.9 + 0.075*scanNum)*plotMax, ant, fontsize='16')
                if spw_index == 0: currentPL.text(np.min(FREQ), 0.8* plotMax, timeLabel, fontsize='6')
                else: currentPL.set_yticklabels([])
        figAnt.savefig(pp, format='pdf')
    for PL in TsysPL: figAnt.delaxes(PL)
    plt.close('all')
    pp.close()
    del(TsysPL)
    del(figAnt)
    return
#
#-------- Plot Tsys and Trx Spectrum
def plotTsys(prefix, antList, spwList, freqList, atmTime, TrxList, TskyList):
    pp = PdfPages('TSYS_' + prefix + '.pdf')
    #-------- Plots for Tsys spectra
    antNum, spwNum, scanNum, polNum  = len(antList), len(spwList), len(atmTime), TrxList[0].shape[0]
    #-------- Prepare Plots
    figAnt = plt.figure(figsize = (8, 11))
    figAnt.suptitle(prefix)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.42, 'Tsys (X:blue, Y:green) and Trec (X:magenta, Y:yellow) [K]', rotation=90)
    #-------- Plot BP
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in TsysPL: figAnt.delaxes(PL)
        #
        TsysPL, TsysMax = [], []
        for spw_index in range(spwNum):
            TsysMax = TsysMax + [1.7* np.percentile(TrxList[spw_index], 75, axis=(0,1,3))[ant_index] + 1.5* np.median(TskyList[spw_index])]
            #TsysMax = TsysMax + [1.7* np.percentile(TrxList[spw_index], 75, axis=(0,1,3))[ant_index] + np.max(TskyList[spw_index])]
        plotMax = max(TsysMax)
        for spw_index in range(spwNum):
            chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
            for scan_index in range(scanNum):
                currentPL = figAnt.add_subplot(scanNum, spwNum, spwNum* scan_index + spw_index + 1 )
                TsysPL = TsysPL + [currentPL]
                timeLabel = au.call_qa_time('%fs' % (atmTime[scan_index]), form='fits')
                for pol_index in range(polNum):
                    plotTrx  = TrxList[spw_index][pol_index, chRange, ant_index, scan_index]
                    plotTsys = TskyList[spw_index][chRange, ant_index, scan_index] + plotTrx
                    currentPL.step( freqList[spw_index][chRange], plotTsys, where='mid', color=polColor[pol_index], label = 'Tsys Pol '+ polName[pol_index])
                    currentPL.step( freqList[spw_index][chRange], plotTrx,  color=polColor[pol_index+2], where='mid', label = 'Trec Pol ' + polName[pol_index])
                #
                currentPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, plotMax])
                currentPL.tick_params(axis='both', labelsize=6)
                if scan_index == 0: currentPL.set_title('SPW %d' % (spwList[spw_index]))
                if scan_index == 0 and spw_index == 0: currentPL.text(1.2* np.min(freqList[0]) - 0.2* np.max(freqList[0]), (0.9 + 0.075*scanNum)*plotMax, antList[ant_index], fontsize='16')
                if scan_index < scanNum - 1: currentPL.set_xticklabels([])
                if spw_index == 0: currentPL.text(np.min(freqList[spw_index]), 0.8* plotMax, timeLabel, fontsize='8')
                else: currentPL.set_yticklabels([])
                #
            #
        #
        figAnt.savefig(pp, format='pdf')
        #
        #
    #
    for PL in TsysPL: figAnt.delaxes(PL)
    plt.close('all')
    pp.close()
    del(TsysPL)
    del(figAnt)
    return
#
#-------- Plot baseline-AllanVariance map
def plotBLAV(prefix, antList, spw, AV_bl):
    pp = PdfPages('AV_%s-SPW%d.pdf' % (prefix, spw))
    figAV = plt.figure(figsize = (8, 8))
    figAV.suptitle('Allan Variance %s SPW%d' % (prefix, spw))
    axAV = figAV.add_subplot(1, 1, 1)
    antNum, blNum = len(antList), AV_bl.shape[1]
    AVmap = np.zeros([2, antNum, antNum])
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        AVmap[0, ants[1], ants[0]] = AV_bl[0, bl_index]  # pol-X
        AVmap[1, ants[0], ants[1]] = AV_bl[1, bl_index]  # pol-Y
    #
    AVmap = np.sum(AVmap, axis=0)
    for ant_index, ant in enumerate(antList): AVmap[ant_index, ant_index] = 0.0
    imAV = axAV.imshow(AVmap, cmap='coolwarm', vmin=0.0, vmax=0.01)
    figAV.colorbar(imAV, ax=axAV)
    axAV.plot([0,antNum-1], [0,antNum-1], color='w')
    axAV.set_xticks(list(range(antNum))); axAV.set_yticks( list(range(antNum)))
    axAV.set_xticklabels(antList.tolist()); axAV.set_yticklabels(antList.tolist())
    figAV.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot autocorrelation power spectra
def plotAC(prefix, antList, spwList, freqList, AC):
    pp = PdfPages('AC_' + prefix + '.pdf')
    antNum, spwNum, polNum = len(antList), len(spwList), AC[0].shape[2]
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Power Spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.5, 'Median amplitude and variation [dB]', rotation=90)
    #-------- Plot AC
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in ACList: figAnt.delaxes(PL)
            for PL in SDList: figAnt.delaxes(PL)
        #
        ACList, SDList = [], []
        ACMAX = 0.0
        for spw_index in range(spwNum): ACMAX = max(ACMAX, np.max(AC[spw_index][ant_index]))
        for spw_index in range(spwNum):
            Freq = freqList[spw_index]
            for pol_index in range(polNum):
                ACPL = figAnt.add_subplot(4, spwNum, spwNum* pol_index + spw_index + 1)
                SDPL = figAnt.add_subplot(4, spwNum, spwNum* (2+pol_index) + spw_index + 1)
                ACList = ACList + [ACPL]; SDList = SDList + [SDPL]
                avgAC = np.mean(AC[spw_index][ant_index, :, pol_index], axis=0)
                plotAC = 10.0* np.log10(avgAC)
                chNum = len(plotAC)
                chRange = range( int(chNum*0.06), int(chNum* 0.96))
                #chRange = range(3050,3090)
                maxAC, minAC, maxFreq = np.max(plotAC[chRange]), np.min(plotAC[chRange]), Freq[np.argmax(plotAC)]
                text_sd = 'Peak = %.1f dB at %.2f GHz' % (maxAC, maxFreq)
                #plotMax, plotMin = 10.0* np.log10(ACMAX), 10.0* np.log10(ACMAX) - 20.0
                plotMax, plotMin = maxAC + 0.1,  minAC - 0.1
                #plotMax, plotMin = 10.0* np.log10(ACMAX), 118.0
                if spw_index == 0 and pol_index == 0: ACPL.text(1.2* np.min(Freq) - 0.2* np.max(Freq), 1.3*plotMax - 0.3* plotMin, antList[ant_index], fontsize=10)
                ACPL.axis([np.min(Freq[chRange]), np.max(Freq[chRange]), plotMin, 1.1* plotMax - 0.1* plotMin])
                ACPL.tick_params(axis='both', labelsize=6)
                ACPL.set_xticklabels([])
                ACPL.text( np.min(Freq), 1.2*plotMax - 0.2*plotMin, 'AC SPW=%d Pol-%s' % (spwList[spw_index], polName[pol_index]), fontsize=7)
                ACPL.text( np.min(Freq), 1.12*plotMax - 0.12*plotMin, text_sd, fontsize=7)
                ACPL.step(Freq[chRange], plotAC[chRange], where='mid')
                #
                plotSD = 10.0* np.log10(np.std(AC[spw_index][ant_index, :, pol_index]/avgAC, axis=0))
                maxSD, minSD, maxFreq = np.max(plotSD[chRange]), np.min(plotSD[chRange]), Freq[np.argmax(plotSD)]
                text_sd = 'Peak = %.1f dB at %.2f GHz' % (maxSD, maxFreq)
                bgcolor = 'green'
                if maxSD > -30.0: bgcolor = 'orange'
                if maxSD > -20.0: bgcolor = 'red'
                #plotMax, plotMin = max(-20.0, maxSD), min(-40.0, minSD)
                #plotMax, plotMin = max(-20.0, maxSD), min(-20.0, minSD)
                plotMax, plotMin = maxSD, minSD
                SDPL.axis([np.min(Freq[chRange]), np.max(Freq[chRange]), 1.05*plotMin - 0.05*plotMax, 1.05*plotMax - 0.05*plotMin])
                SDPL.axhspan(ymin=-30.0, ymax=plotMax, color=bgcolor, alpha=0.1) 
                SDPL.tick_params(axis='both', labelsize=6)
                SDPL.get_xaxis().get_major_formatter().set_useOffset(False)
                if pol_index == 0: SDPL.set_xticklabels([])
                SDPL.text( np.min(Freq), 1.16* plotMax - 0.16*plotMin, 'SD SPW=%d Pol-%s' % (spwList[spw_index], polName[pol_index]), fontsize=7)
                SDPL.text( np.min(Freq), 1.08* plotMax - 0.08*plotMin, text_sd, fontsize=7)
                SDPL.step(Freq[chRange], plotSD[chRange], where='mid')
            #
        #
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(ACList)
    del(SDList)
    del(ACPL)
    del(SDPL)
    return
#
#-------- Plot Cross power spectrum
def plotSP(pp, prefix, antList, spwList, freqList, BPList, plotMin=0.0, plotMax=1.5, delayMessage=False):
    antNum, spwNum = len(antList), len(spwList)
    delayAnt = np.zeros([antNum, spwNum, 2])
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Bandpass')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #-------- Plot BP
    text_delay = 'Ant  '
    for spw_index, spw in enumerate(spwList):
        text_delay = text_delay + ' |SPW%02d'%spw
        text_delay = text_delay + '  ' + polName[0]
        if BPList[spw_index].shape[1] == 2: text_delay = text_delay + '         ' + polName[1]
    if delayMessage: print(text_delay)
    text_delay = '------+'
    for spw_index, spw in enumerate(spwList): text_delay = text_delay + '--- delay [ns] ----+'
    if delayMessage: print(text_delay)
    for ant_index, antName in enumerate(antList):
        if ant_index > 0:
            for PL in AmpList: figAnt.delaxes(PL)
            for PL in PhsList: figAnt.delaxes(PL)
        #
        AmpList, PhsList = [], []
        text_delay = antName + '  |'
        for spw_index, spw in enumerate(spwList):
            Freq = freqList[spw_index]
            plotFreq = 1.0e-9 * Freq
            chNum = len(Freq); chRange = list(range(int(0.10*chNum), int(0.95*chNum)))
            BW = Freq[chRange[-1]] - Freq[chRange[0]]
            AmpPL = figAnt.add_subplot(2, spwNum, spw_index + 1 )
            PhsPL = figAnt.add_subplot(2, spwNum, spwNum + spw_index + 1 )
            AmpList = AmpList + [AmpPL]
            PhsList = PhsList + [PhsPL]
            ppolNum = BPList[spw_index].shape[1]
            for pol_index in list(range(ppolNum)):
                plotBandpass = BPList[spw_index][ant_index,pol_index]
                delay_ant, delaySNR = delay_search(plotBandpass[chRange])
                delayAnt[ant_index, spw_index, pol_index] = 0.5e9* delay_ant / BW   # delay in [ns]
                AmpPL.fill_between(plotFreq[chRange], plotMin, plotMax, color='yellow', alpha=0.1)
                AmpPL.step(plotFreq, abs(plotBandpass), color=polColor[pol_index], where='mid', label = 'Pol-%s' % (polName[pol_index]))
                text_delay = text_delay + '  %+.4f ' % (0.5e9* delay_ant/ BW)
                PhsPL.fill_between(plotFreq[chRange], -np.pi, np.pi, color='yellow', alpha=0.1)
                PhsPL.plot(plotFreq, np.angle(plotBandpass), '.', color=polColor[pol_index], label = 'Pol-%s %+.3f ns' % (polName[pol_index], -0.5e9* delay_ant / BW))
            #
            if spw_index == 0: AmpPL.set_title(antList[ant_index])
            AmpPL.axis([np.min(plotFreq), np.max(plotFreq), plotMin, plotMax])
            AmpPL.tick_params(axis='both', labelsize=6)
            AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            AmpPL.text( np.min(plotFreq), 0.9* plotMax + 0.1* plotMin, 'SPW=%d Amp' % (spwList[spw_index]))
            PhsPL.axis([np.min(plotFreq), np.max(plotFreq), -np.pi, np.pi])
            PhsPL.tick_params(axis='both', labelsize=6)
            PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            PhsPL.text( np.min(plotFreq), 2.5, 'SPW=%s Phs' % (str(spw)))
        #
        if delayMessage: print(text_delay)
        figAnt.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return delayAnt
#
#-------- Plot Bandpass
def plotBP(pp, prefix, antList, spwList, BPscan, BPList, bunchNum=1, plotMax=1.2, plotMarker=[[]]):
    msfile = prefix + '.ms'
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Scan %d' % (BPscan))
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #-------- Plot BP
    for ant_index in list(range(antNum)):
        if ant_index > 0:
            for PL in AmpList: figAnt.delaxes(PL)
            for PL in PhsList: figAnt.delaxes(PL)
        #
        AmpList, PhsList = [], []
        for spw_index in list(range(spwNum)):
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
            AmpPL = figAnt.add_subplot(2, spwNum, spw_index + 1 )
            PhsPL = figAnt.add_subplot(2, spwNum, spwNum + spw_index + 1 )
            AmpList = AmpList + [AmpPL]
            PhsList = PhsList + [PhsPL]
            for pol_index in list(range(BPList[spw_index].shape[1])):
                plotBandpass = BPList[spw_index][ant_index,pol_index]
                '''
                delayAnt, delaySNR = delay_search(plotBandpass)
                text_sd = '%s %s %+f %.1f' % (antList[ant_index], polName[pol_index], 0.5*delayAnt / (chNum*chWid), delaySNR)
                print(text_sd)
                '''
                AmpPL.step(Freq, abs(plotBandpass), color=polColor[pol_index], where='mid', label = 'Pol-' + polName[pol_index])
                PhsPL.plot( Freq, np.angle(plotBandpass), '.', color=polColor[pol_index], label = 'Pol-' + polName[pol_index])
            #
            if len(plotMarker[0]) > 0: 
                for spurIndex in list(range(len(plotMarker[spw_index]))): AmpPL.vlines(x=1.0e-9 * plotMarker[spw_index][spurIndex], ymin=0.0, ymax=1.25*plotMax, color='gray') 
            if spw_index == 0: AmpPL.set_title(antList[ant_index])
            AmpPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25*plotMax])
            AmpPL.tick_params(axis='both', labelsize=6)
            AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            AmpPL.text( np.min(Freq), 1.1* plotMax, 'SPW=%d Amp' % (spwList[spw_index]))
            PhsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            PhsPL.tick_params(axis='both', labelsize=6)
            PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            PhsPL.text( np.min(Freq), 2.5, 'SPW=%d Phs' % (spwList[spw_index]))
        #
        #plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(AmpList)
    del(PhsList)
    del(AmpPL)
    del(PhsPL)
    del(figAnt)
    return
#
#-------- Plot Gain
def plotGain(prefix, spw, plotAntList=[]):
    #-------- Load tables
    antFile = prefix + '.Ant.npy'; antList = np.load(antFile)
    if plotAntList == []: plotAntList = antList
    plotAntNum = len(plotAntList)
    columnNum  = max(int(np.floor(np.log(plotAntNum))), 1)
    antMap = indexList(plotAntList, antList)
    timeFile = '%s-SPW%d.TS.npy' % (prefix, spw)
    GainFile = '%s-SPW%d.GA.npy' % (prefix, spw)    # Gain[ant, pol, time]
    FlagFile = '%s-SPW%d.FG.npy' % (prefix, spw)    # Flag[ant, time] 
    pp = PdfPages('GA_%s-SPW%d.pdf' %  (prefix, spw))
    DT, timeStamp, Gain, FG = [], np.load(timeFile), np.load(GainFile), np.load(FlagFile)
    polList = polName[:Gain.shape[1]]
    for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    #-------- Prepare Plots
    figAmp, figPhs = plt.figure(figsize = (8, 11)), plt.figure(figsize = (8, 11))
    figAmp.suptitle(GainFile + ' Gain Amplitude'); figPhs.suptitle(GainFile + ' Gain Phase')
    figAmp.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d'))); figPhs.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
    figAmp.text(0.03, 0.45, 'Gain Amplitude = sqrt(correlated flux / SEFD)', rotation=90); figPhs.text(0.03, 0.45, 'Gain Phase [deg]', rotation=90)
    plotMin, plotMax = 0.0, 1.1* np.percentile(np.max(abs(Gain[antMap]), axis=1)* FG[antMap], 90)
    #-------- Plot Gain
    gainXmed, gainYmed = [], []
    for plot_index, ant in enumerate(plotAntList):
        ant_index = np.where(antList == ant)[0][0]
        flag_index = np.where(FG[ant_index] > 0.01)[0].tolist()
        ampMed = [0.0, 0.0]
        AmpPL = figAmp.add_subplot( int(np.ceil(plotAntNum/columnNum)), columnNum, plot_index + 1 )
        PhsPL = figPhs.add_subplot( int(np.ceil(plotAntNum/columnNum)), columnNum, plot_index + 1 )
        for pol_index, pol in enumerate(polList): AmpPL.plot( np.array(DT)[flag_index], abs(Gain[ant_index, pol_index, flag_index]), '.', markersize=3, color=polColor[pol_index], label=pol)
        for pol_index, pol in enumerate(polList): PhsPL.plot( np.array(DT)[flag_index], np.angle(Gain[ant_index, pol_index, flag_index]*Gain[ant_index, pol_index, 0].conjugate())*RADDEG, '.', markersize=3, color=polColor[pol_index], label=pol)
        if len(flag_index) > 0: ampMed = np.median(abs(Gain[ant_index][:, flag_index]), axis=1)
        gainXmed, gainYmed = gainXmed + [ampMed[0]], gainYmed + [ampMed[1]]
        AmpPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
        AmpPL.yaxis.offsetText.set_fontsize(3)
        PhsPL.yaxis.offsetText.set_fontsize(3)
        AmpPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        AmpPL.tick_params(labelsize=4)
        PhsPL.tick_params(labelsize=4)
        AmpPL.axis([np.min(DT), np.max(DT), plotMin, plotMax])
        PhsPL.axis([np.min(DT), np.max(DT), -180.0, 180.0])
        text_sd = '%s : Gain(median) = (' % (antList[ant_index])
        for pol_index, pol in enumerate(polList): text_sd = text_sd + '%.2f%% ' % (100.0* ampMed[pol_index])
        text_sd = text_sd[:-1] + ')'
        print(text_sd)
        AmpPL.text( 0.05, 1.02, text_sd, transform=AmpPL.transAxes, fontsize=5)
        PhsPL.text( 0.05, 1.02, antList[ant_index], transform=PhsPL.transAxes, fontsize=5)
        if plot_index == 0:
            AmpPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
            PhsPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
    #
    text_sd = 'all: Gain(median) = (%.2f%% %.2f%%)' % (100.0* np.median(gainXmed), 100.0* np.median(gainYmed)); print(text_sd)
    figAmp.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/checkGain.py to generate this plot', fontsize=4)
    figPhs.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/checkGain.py to generate this plot', fontsize=4)
    figAmp.savefig(pp, format='pdf'); figPhs.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot visibility amplitude and closure phase
def plotBispec(antList, scanVis, DT, plotFile, labelList, pMax):
    antNum = len(antList)
    polNum, blNum  = scanVis.shape[0], scanVis.shape[1]
    figInch  = max(16,antNum)
    fontSize = min(32, figInch)
    figSPW = plt.figure(figsize=(figInch, figInch))
    figSPW.text(0.475, 0.05, labelList[0], fontsize=fontSize)
    figSPW.text(0.05, 0.5, labelList[1], rotation=90, fontsize=fontSize)
    figSPW.text(0.95, 0.5, labelList[2], rotation=-90, fontsize=fontSize)
    figSPW.suptitle(labelList[3], fontsize=fontSize)
    #
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        #-------- Plot visibility amplitude
        BLamp = figSPW.add_subplot(antNum-1, antNum-1, ants[1]*(antNum -1) + ants[0])
        for pol_index in list(range(polNum)):
            plotVis = scanVis[pol_index, bl_index]
            BLamp.step(DT, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
        #
        BLamp.axis([np.min(DT), np.max(DT), 0.0, 1.25*pMax])
        BLamp.xaxis.set_major_locator(plt.NullLocator())
        if bl_index == 0:
            BLamp.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        if ants[0] - ants[1] == 1:
            BLamp.set_ylabel(antList[ants[1]])
        if ants[1] == 0:    # Antenna label in the top and leftside
            BLamp.set_title( antList[ants[0]] )
        if ants[0] == antNum - 1:    # Antenna at rightside
            BLamp.yaxis.tick_right()
        else:
            BLamp.yaxis.set_major_locator(plt.NullLocator())
        #-------- Plot closure phase
        if bl_index < 2: continue
        if ants[1] > 0:         # plot Closure phase
            BLphs = figSPW.add_subplot(antNum-1, antNum-1, (ants[0] - 1)*(antNum-1) + ants[1])
            BLphs.patch.set_facecolor('lightyellow')
            ant1List, ant0List = [], []
            tri = [Ant2Bl(ants1,ants0) for (ants1,ants0) in zip([ants[0],ants[1],ants[0]],[0, 0, ants[1]])]
            for pol_index in list(range(polNum)):
                plotVis = scanVis[pol_index, tri[0]].conjugate()* scanVis[pol_index, tri[1]]* scanVis[pol_index, tri[2]]
                #plotVis = scanVis[pol_index, tri0].conjugate()* scanVis[pol_index, tri1]* scanVis[pol_index, tri2]
                BLphs.plot(DT, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
            #
            print('%d : %d - %d - %d (ant %s, %s, %s)' % (bl_index, tri[0], tri[1], tri[2], antList[0], antList[ants[1]], antList[ants[0]]))
            #print('%d : %d - %d - %d (ant %s, %s, %s)' % (bl_index, tri0, tri1, tri2, antList[0], antList[ants[1]], antList[ants[0]]))
            BLphs.set_title('%s-%s-%s' % (antList[0], antList[ants[1]], antList[ants[0]] ), fontsize=0.5*fontSize)
            BLphs.axis([np.min(DT), np.max(DT), -math.pi, math.pi])
            BLphs.tick_params(axis='x', labelsize=int(fontSize*0.25), labelrotation=-90)
        if ants[1] == 1: BLphs.set_ylabel(antList[ants[0]] )
        if ants[1] > 1: BLphs.yaxis.set_major_locator(plt.NullLocator())
        if ants[0] < antNum - 1:    # except bottom panel : skip drawing X-axis
            BLphs.xaxis.set_major_locator(plt.NullLocator())
        #
    #
    figSPW.savefig(plotFile + '.pdf', format='pdf', dpi=144)
    plt.close('all')
    return
#
#-------- Plot XY-phase spectra
def plotXYP(pp, prefix, spwList, XYspec, XYdelay, bunchNum=1):
    spwNum = len(spwList)
    figXYP = plt.figure(figsize = (11, 8))
    figXYP.suptitle(prefix + ' XY Phase')
    figXYP.text(0.45, 0.05, 'Frequency [GHz]')
    figXYP.text(0.03, 0.45, 'XY Phase [deg]', rotation=90)
    for spw_index, spw in enumerate(spwList):
        chNum, chWid, Freq = GetChNum(prefix + '.ms', spw); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
        PhsPL = figXYP.add_subplot(1, spwNum, spw_index + 1)
        XYP  = XYspec[spw_index]
        PhsPL.plot( Freq, np.angle(XYP)*RADDEG, '.', label = 'SPW%2d: XYdelay=%+.3f ns' % (spw, -XYdelay[spw_index]))
        PhsPL.axis([np.min(Freq), np.max(Freq), -180.0, 180.0])
        PhsPL.tick_params(axis='both', labelsize=6)
        PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #
    #plt.show()
    figXYP.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot D-term spectra
def plotDSpec(pp, prefix, antList, spwList, FreqList, DxList, DyList):
    plotMax = 0.12
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' D-term spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
    for ant_index in list(range(antNum)):
        if ant_index > 0:
            for PL in DxPList: figAnt.delaxes(PL)
            for PL in DyPList: figAnt.delaxes(PL)
        #
        DxPList, DyPList = [], []
        for spw_index in list(range(spwNum)):
            DxPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            DyPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            DxPList = DxPList + [DxPL]
            DyPList = DyPList + [DyPL]
            #
            plotDx, plotDy = DxList[ant_index + antNum* spw_index], DyList[ant_index + antNum* spw_index]
            DxPL.step( FreqList[spw_index], plotDx.real, where='mid', label = 'reDx')
            DxPL.step( FreqList[spw_index], plotDx.imag, where='mid', label = 'imDx')
            DxPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            DyPL.step( FreqList[spw_index], plotDy.real, where='mid', label = 'reDy')
            DyPL.step( FreqList[spw_index], plotDy.imag, where='mid', label = 'imDy')
            DyPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            #
            if spw_index == 0: DxPL.set_title(antList[ant_index])
            DxPL.tick_params(axis='both', labelsize=6)
            DyPL.tick_params(axis='both', labelsize=6)
        #
        DxPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        DyPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        #plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all'); pp.close()
    del(DxPList); del(DyPList); del(DxPL); del(DyPL); del(figAnt)
    return
#
#-------- Plot Flagging map
def plotFlag(pp, prefix, antList, DT, spwList, FGList):
    antNum = len(antList)
    figFG = plt.figure(figsize = (8, 11))
    figFG.suptitle('%s Flag Map' % (prefix));
    figFG.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    for spw_index, spw in enumerate(spwList):
        FGPL = figFG.add_subplot(len(spwList), 1, spw_index + 1 )
        FG = FGList[spw_index]
        for ant_index, ant in enumerate(antList):
            FGPL.scatter(DT, np.linspace(ant_index, ant_index, len(DT)), vmin=0, vmax=1, marker='s', s=1, c= 1.0 - FG[ant_index], cmap=cm.cool)
        FGPL.set_title('SPW%d' % (spw))
        FGPL.tick_params(labelsize=3)
        FGPL.set_yticks(list(range(antNum)))
        FGPL.set_yticklabels(antList.tolist())
    #
    figFG.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot Stokes parameters vs uv distance
def plotFL(pp, scanDic, SPWDic):
    polLabel = ['I', 'Q', 'U', 'V']
    Pcolor   = ['black', 'blue', 'red', 'green']
    for scan_index, scan in enumerate(scanDic.keys()):
        scanFlag = scanDic[scan]['scanFlag']
        if len(scanFlag) == 0: continue
        figFL, axes = plt.subplots(3, 4, figsize = (11, 8))
        figFL.suptitle(scanDic[scan]['msfile'][:-3])
        figFL.text(0.45, 0.05, 'Projected baseline [m]')
        UVW = scanDic[scan]['UVW'][:,:,scanFlag]; uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        sourceName = scanDic[scan]['source']
        text_src  = 'Scan%02d %010s EL=%4.1f deg' % (scan, sourceName, 180.0* np.median(scanDic[scan]['EL'][scanFlag])/np.pi)
        timeLabel = qa.time('%fs' % np.median(scanDic[scan]['mjdSec'][scanFlag]), form='ymd')[0] + ' SA=%.1f' % (scanDic[scan]['SA']) + ' deg.'
        #---- Display results
        uvMin, uvMax, IMax = min(uvDist), max(uvDist), np.max(scanDic[scan]['I'])
        axes[0,0].text(0.0, 1.15*180, text_src)
        axes[0,3].text(0.0, 1.15*180, timeLabel)
        spwNum = len(scanDic[scan]['I'])
        for spw_index, spw in enumerate(SPWDic['spw']):
            StokesVis = scanDic[scan]['StokesVis'][spw_index]
            ScanFlux, ScanSlope, ErrFlux = lmStokes(StokesVis, uvDist)
            #-------- Plot Stokes visibilities
            axes[0,spw_index].plot( np.array([0.0, uvMax]), np.array([0.0, 0.0]), '-', color='grey')        # phase-0 line
            axes[1,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[0], ScanFlux[0]+ uvMax* ScanSlope[0]]), '-', color=Pcolor[0])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[1], ScanFlux[1]+ uvMax* ScanSlope[1]]), '-', color=Pcolor[1])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[2], ScanFlux[2]+ uvMax* ScanSlope[2]]), '-', color=Pcolor[2])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[3], ScanFlux[3]+ uvMax* ScanSlope[3]]), '-', color=Pcolor[3])
            #
            axes[0,spw_index].plot( uvDist, 180.0* np.angle(StokesVis[0])/ np.pi, '.', label=polLabel[0], color=Pcolor[0])  # Plot visibility phase
            axes[1,spw_index].plot( uvDist, StokesVis[0].real, '.', label=polLabel[0], color=Pcolor[0])     # Plot Stokes I
            axes[2,spw_index].plot( uvDist, StokesVis[1].real, '.', label=polLabel[1], color=Pcolor[1])     # Plot Stokes Q
            axes[2,spw_index].plot( uvDist, StokesVis[2].real, '.', label=polLabel[2], color=Pcolor[2])     # Plot Stokes U
            axes[2,spw_index].plot( uvDist, StokesVis[3].real, '.', label=polLabel[3], color=Pcolor[3])     # Plot Stoees V
            axes[0,spw_index].axis([0.0, uvMax, -180, 180])
            axes[1,spw_index].axis([0.0, uvMax, 0.0, 1.25*IMax])
            axes[2,spw_index].axis([0.0, uvMax, -0.25*IMax, 0.25*IMax])
            axes[0,spw_index].text(0.0, 1.02*180, 'SPW%2d %5.1f GHz' % (spw, 1.0e-9* np.median(SPWDic['freq'][spw_index])))
        #
        axes[0,0].set_ylabel('Phase [deg]')
        axes[1,0].set_ylabel('Stokes I [Jy]')
        axes[2,0].set_ylabel('Stokes Q,U,V [Jy]')
        axes[1,spw_index].legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        axes[2,spw_index].legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        #plt.show()
        figFL.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#-------- Plot Stokes Spectra
def plotStokesSpectra(prefix, refantName, sourceName, spw):
    polLabel, Pcolor = ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
    StokesSpecFile = '%s-REF%s-%s-SPW%d.StokesSpec.npy' % (prefix, refantName, sourceName, spw)
    StokesFreqFile = '%s-REF%s-%s-SPW%d.Freq.npy' % (prefix, refantName, sourceName, spw)
    if not os.path.isfile(StokesSpecFile): return
    StokesSpec = np.load(StokesSpecFile)
    Freq = np.load(StokesFreqFile)
    chNum = len(Freq)
    sdF = 1.0 / np.sqrt(chNum)
    StokesFlux = np.mean(StokesSpec, axis=1)
    StokesErr  = sdF* np.std(StokesSpec, axis=1)
    #text_sd = '%12s %6.2f GHz   %6.3f (%.3f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %+7.2f (%.2f) ' % (sourceName, np.median(Freq), np.mean(StokesSpec[0]), sdF*np.std(StokesSpec[0]), np.mean(StokesSpec[1]), sdF* np.std(StokesSpec[1]),np.mean(StokesSpec[2]), sdF* np.std(StokesSpec[2]),np.mean(StokesSpec[3]), sdF* np.std(StokesSpec[3]), )
    text_sd = '%12s %6.2f GHz   %6.3f (%.3f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %7.4f (%.4f) %+7.2f (%.2f) ' % (sourceName, np.median(Freq), StokesFlux[0], StokesErr[0], StokesFlux[1], StokesErr[1], StokesFlux[2], StokesErr[2], StokesFlux[3], StokesErr[3], 100.0*np.sqrt(StokesFlux[1]**2 + StokesFlux[2]**2)/StokesFlux[0], 100.0*np.sqrt(StokesFlux[1]**2 + StokesFlux[2]**2)*np.sqrt( (StokesErr[0]/StokesFlux[0])**2 + (StokesErr[1]/StokesFlux[1])**2 + (StokesErr[2]/StokesFlux[2])**2 ) / StokesFlux[0]  , 0.5* RADDEG* np.arctan2(StokesFlux[2], StokesFlux[1]), 0.5* RADDEG* abs(StokesErr[1]* StokesFlux[2] + StokesErr[2]* StokesFlux[1]) / (abs(StokesFlux[1])**2 + abs(StokesFlux[2])**2))
    print(text_sd)
    figSP = plt.figure(figsize = (11, 8))
    figSP.suptitle(prefix + ' ' + sourceName)
    figSP.text(0.45, 0.05, 'Frequency [GHz]')
    figSP.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
    StokesI_SP = figSP.add_subplot( 2, 1, 1 )
    StokesP_SP = figSP.add_subplot( 2, 1, 2 )
    IMax = np.max(StokesSpec[0])
    StokesI_SP.step(Freq, StokesSpec[0], where='mid', label='Stokes %s = %.3f (%.3f) Jy' % (polLabel[0], np.mean(StokesSpec[0]), sdF* np.std(StokesSpec[0])), color=Pcolor[0])
    StokesP_SP.step(Freq, StokesSpec[1], where='mid', label='Stokes %s = %.4f (%.4f) Jy' % (polLabel[1], np.mean(StokesSpec[1]), sdF* np.std(StokesSpec[1])), color=Pcolor[1])
    StokesP_SP.step(Freq, StokesSpec[2], where='mid', label='Stokes %s = %.4f (%.4f) Jy' % (polLabel[2], np.mean(StokesSpec[2]), sdF* np.std(StokesSpec[2])), color=Pcolor[2])
    StokesP_SP.step(Freq, StokesSpec[3], where='mid', label='Stokes %s = %.4f (%.4f) Jy' % (polLabel[3], np.mean(StokesSpec[3]), sdF* np.std(StokesSpec[3])), color=Pcolor[3])
    StokesI_SP.tick_params(axis='both', labelsize=6)
    StokesP_SP.tick_params(axis='both', labelsize=6)
    StokesI_SP.axis([np.min(Freq), max(Freq), 0.0, 1.25*IMax])
    StokesP_SP.axis([np.min(Freq), max(Freq), -0.10*IMax, 0.10*IMax])
    StokesI_SP.text(min(Freq), IMax*1.35, sourceName)
    StokesI_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    StokesP_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
    figSP.savefig('SP_%s-REF%s-%s-SPW%d.pdf' % (prefix, refantName, sourceName, spw))
    plt.close('all')
    return
#-------- Plot XY cross correlation diagram
def plotXYVis(pp, scanDic, DdotP, DdotM):
    cmap = plt.get_cmap("tab20")
    figXY = plt.figure(figsize = (11,8))
    XYP = figXY.add_subplot( 1, 1, 1 )
    plotMax = 0.0
    for scan_index, scan in enumerate(scanDic.keys()):
        plotY = DdotP + DdotM* scanDic[scan]['QCpUS'] + scanDic[scan]['UCmQS']* np.exp((1.0j)*scanDic[scan]['XYphase'])
        XYP.plot( scanDic[scan]['UCmQS'], plotY.real, '-', color=cmap(2* scan_index))
        XYP.plot( scanDic[scan]['UCmQS'], plotY.imag, '-', color=cmap(2* scan_index + 1))
        XYP.plot( scanDic[scan]['UCmQS'], np.mean(scanDic[scan]['scanVis'][[1,2]].real, axis=0), '.', color=cmap(2* scan_index), label='Scan%d %s' % (scan, scanDic[scan]['source']))
        XYP.plot( scanDic[scan]['UCmQS'], 0.5*np.diff(scanDic[scan]['scanVis'][[2,1]].imag, axis=0)[0], '.', color=cmap(2* scan_index + 1))
        plotMax = max(plotMax, np.max(abs(scanDic[scan]['UCmQS'])))
    XYP.set_xlabel('U $\cos 2 \psi $ - Q $\sin 2 \psi$'); XYP.set_ylabel('<XY*>'); XYP.set_title(scanDic[scan]['msfile'])
    XYP.axis([-1.25*plotMax, 1.25*plotMax, -1.25*plotMax, 1.25*plotMax]); XYP.grid()
    box = XYP.get_position()
    XYP.set_position([box.x0, box.y0, box.width* 0.75, box.height])
    XYP.legend(loc = 'upper left', prop={'size' :8}, numpoints = 1, bbox_to_anchor = (1,1))
    figXY.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#-------- Plot XY cross correlation versus QU
def plotQUXY(pp, scanDic):
    figQU = plt.figure(figsize = (8, 11))
    figQU.text(0.25, 0.05, 'Linear polarization angle w.r.t. X-Feed [deg]')
    sourceList = list(set([scanDic[scan]['source'] for scan in scanDic.keys()]))
    ParaPolPL  = figQU.add_subplot( 2, 1, 1 )
    CrosPolPL  = figQU.add_subplot( 2, 1, 2 )
    cmap = plt.get_cmap("tab10")
    plotMax = 0.0
    for source_index, sourceName in enumerate(sourceList):
        if sourceName in SSOCatalog: continue
        scanList = [scan for scan in scanDic.keys() if 'scanVis' in scanDic[scan].keys() and scanDic[scan]['source'] == sourceName]
        if len(scanList) < 1: continue
        ThetaPlot, VisXX, VisYY, VisXY = [], [], [], []
        #-------- PA range to draw model
        for scan in scanList:
            Stokes = np.mean(scanDic[scan]['scanVis'], axis=1)
            EVPA = 0.5* np.arctan2(Stokes[2], Stokes[1])
            ThetaPlot = ThetaPlot + np.arctan(np.tan(scanDic[scan]['PA'] - EVPA)).tolist()
        ThetaPlot = np.array(ThetaPlot)
        ThetaMin, ThetaMax = np.min(ThetaPlot), np.max(ThetaPlot)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        PArange = ThetaRange + EVPA
        #-------- Draw model
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCmQS, QCpUS = Stokes[2]*CSrange - Stokes[1]* SNrange, Stokes[1]*CSrange + Stokes[2]* SNrange
        ParaPolPL.plot( RADDEG* ThetaRange, QCpUS, color=cmap(source_index), linestyle='dashed')
        ParaPolPL.plot( RADDEG* ThetaRange,-QCpUS, color=cmap(source_index), linestyle='dashdot')
        CrosPolPL.plot( RADDEG* ThetaRange, UCmQS, color=cmap(source_index), linestyle='solid')
        CrosPolPL.plot( RADDEG* ThetaRange, np.zeros(len(ThetaRange)), color=cmap(source_index), linestyle='dotted')
        #-------- Plot data points
        for scan in scanList:
            visChav = scanDic[scan]['visChav']   # visChav[Stokes, time]
            VisXX = VisXX + (abs(scanDic[scan]['visChav'][0]) - abs(np.mean(scanDic[scan]['visChav'][[0,3]], axis=0))).tolist()
            VisYY = VisYY + (abs(scanDic[scan]['visChav'][3]) - abs(np.mean(scanDic[scan]['visChav'][[0,3]], axis=0))).tolist()
            VisXY = VisXY + np.mean(scanDic[scan]['visChav'][[1,2]], axis=0).tolist()
            text_sd = 'Scan %d : %s' % (scan, qa.time('%fs' % (scanDic[scan]['mjdSec'][0]), form='fits', prec=6)[0][11:21])
            text_PA = RADDEG* np.arctan(np.tan((scanDic[scan]['PA'][0] - EVPA)))
            rotText, posText = 90, 'bottom'
            if text_PA < 0: rotText, posText = -90, 'top'
            CrosPolPL.text(text_PA, 0, text_sd, verticalalignment=posText, fontsize=6, rotation=rotText)
        VisXX, VisYY, VisXY, QCpUS, UCmQS = np.array(VisXX), np.array(VisYY), np.array(VisXY), np.array(QCpUS), np.array(UCmQS)
        ParaPolPL.plot( RADDEG* ThetaPlot, VisXX, 'x', fillstyle='full', color=cmap(source_index), label=sourceName + ' XX')
        ParaPolPL.plot( RADDEG* ThetaPlot, VisYY, '1', fillstyle='none', color=cmap(source_index), label=sourceName + ' YY')
        CrosPolPL.plot( RADDEG* ThetaPlot, VisXY.real, 'o', fillstyle='none', mew=0.2, color=cmap(source_index), label=sourceName + ' re')
        CrosPolPL.plot( RADDEG* ThetaPlot, VisXY.imag, ',', color=cmap(source_index), label=sourceName + ' im')
        plotMax = np.max([plotMax, np.max(abs(VisXX))])
        plotMax = np.max([plotMax, np.max(abs(VisYY))])
        plotMax = np.max([plotMax, np.max(abs(VisXY))])
    #
    ParaPolPL.axis([-90.0, 90.0, -1.2*plotMax, 1.2*plotMax])
    ParaPolPL.grid(linewidth=0.5, linestyle='dotted')
    ParaPolPL.set_ylabel('Parallel-hand correlation - Stokes I [Jy]')
    CrosPolPL.axis([-90.0, 90.0, -1.2*plotMax, 1.2*plotMax])
    CrosPolPL.grid(linewidth=0.5, linestyle='dotted')
    CrosPolPL.set_ylabel('Cross-hand correlation [Jy]')
    box = ParaPolPL.get_position()
    ParaPolPL.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ParaPolPL.legend(loc = 'upper right', prop={'size' :4.5}, ncol=2, numpoints=1, labelspacing=0.05, bbox_to_anchor = (1.3,1))
    box = CrosPolPL.get_position()
    CrosPolPL.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    CrosPolPL.legend(loc = 'upper right', prop={'size' :4.5}, ncol=2, numpoints=1, labelspacing=0.05, bbox_to_anchor = (1.3,1))
    figQU.suptitle('Cross-Pol. Correlations %s' % (scanDic[scan]['msfile']))
    figQU.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#

