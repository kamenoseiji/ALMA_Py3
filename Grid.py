import math
import numpy as np
import scipy
#SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOscore   = [
[ 5.0,     4.0,       1.0,        1.0,        0.1,     0.2,  0.3,      0.2,     0.1,     0.1,    10.0,   10.0,     10.0],   # Band 1
[ 6.0,     5.0,       1.0,        1.0,        0.1,     0.6,  0.5,      0.3,     0.2,     0.2,    10.0,   10.0,     10.0],   # Band 2
[ 7.0,     6.0,       1.0,        1.0,        0.2,     0.7,  0.7,      0.5,     0.3,     0.3,    10.0,   10.0,     10.0],   # Band 3
[ 8.0,     7.0,       2.0,        2.0,        0.3,     0.9,  0.9,      0.6,     0.4,     0.4,     8.0,   10.0,      8.0],   # Band 4
[ 9.0,     8.0,       3.0,        3.0,        0.5,     1.0,  1.0,      0.8,     0.5,     0.5,     4.0,    4.0,      4.0],   # Band 5
[10.0,     9.0,       4.0,        4.0,        5.0,     3.0,  3.0,      1.0,     0.6,     0.6,     2.0,    2.0,      2.0],   # Band 6
[10.0,     9.0,       5.0,        5.0,        7.0,     4.0,  4.0,      3.0,     0.7,     0.7,     1.0,    1.0,      1.0],   # Band 7
[ 5.0,     6.0,       7.0,        7.0,        7.0,     8.0,  4.0,      3.0,     0.7,     0.7,     1.0,    1.0,      1.0],   # Band 8
[ 2.0,     3.0,       7.0,        7.0,        7.0,     9.0,  4.0,     10.0,     7.0,     4.0,     1.0,    1.0,      1.0]]   # Band 9
ELshadow = math.pi* 40.0 / 180.0
SunAngleThresh = 5.0
#-- Acceptable opacity       B3   B4    B5   B6   B7   B8,  B9, B10
TauLimit = [0.0,  0.1, 0.1, 0.2, 0.2, 0.25, 0.5, 0.5, 0.5, 0.5, 0.6]
#-- Source name alias
sourceDic = {
'J0006-063':'J0006-0623',
'0006-063':'J0006-0623',
'J0237+288':'J0237+2848',
'J0238+166':'J0238+1636',
'3c84':'J0319+4130',
'0334-401':'J0334-4008',
'J0334-401':'J0334-4008',
'J0423-013':'J0423-0120',
'J0510+180':'J0510+1800',
'J0519-454':'J0519-4546',
'J0522-364':'J0522-3627',
'J0538-440':'J0538-4405',
'J0750+125':'J0750+1231',
'J0854+201':'J0854+2006',
'J1037-295':'J1037-2934',
'J1058+015':'J1058+0133',
'J1107-448':'J1107-4449',
'J1146+399':'J1146+3958',
'3c273':'J1229+0203',
'3c279':'J1256-0547',
'3C279':'J1256-0547',
'J1337-129':'J1337-1257',
'1337-129':'J1337-1257',
'J1427-421':'J1427-4206',
'1427-421':'J1427-4206',
'J1517-243':'J1517-2422',
'1517-243':'J1517-2422',
'J1550+054':'J1550+0527',
'1550+054':'J1550+0527',
'J1613-586':'J1617-5848',
'1613-586':'J1617-5848',
'J1625-254':'J1625-2527',
'3c345':'J1642+3948',
'J1733-130':'J1733-1304',
'1733-130':'J1733-1304',
'nrao530':'J1733-1304',
'J1751+096':'J1751+0939',
'1751+096':'J1751+0939',
'J1924-292':'J1924-2914',
'1924-292':'J1924-2914',
'J2025+337':'J2025+3343',
'J2056-472':'J2056-4714',
'2056-472':'J2056-4714',
'2157-694':'J2157-6941',
'J2148+069':'J2148+0657',
'J2232+117':'J2232+1143',
'3c454.3':'J2253+1608',
'J2258-279':'J2258-2758',
'BL_Lac':'J2202+4216',
'VLA_J1650+0824':'J1650+0824',
'VLA_J1652+0618':'J1652+0618',
'VLA_J1656+1826':'J1656+1826',
'VLA_J1658+0741':'J1658+0741',
'VLA_J1658+0515':'J1658+0515',
'VLA_J1700+0522':'J1700+0522',
'VLA_J1707+0148':'J1707+0148',
'VLA_J1707+1331':'J1707+1331',
'VLA_J1707+1846':'J1707+1846',
'VLA_J1708+0035':'J1708+0035',
'VLA_J1716+2152':'J1716+2152',
'VLA_J1719+0658':'J1719+0658',
'VLA_J1719+1745':'J1719+1745',
'VLA_J1719+0817':'J1719+0817',
'VLA_J1722+1013':'J1722+1013',
'VLA_J1726+0639':'J1726+0639',
'VLA_J1728+1215':'J1728+1215',
'VLA_J1728+0427':'J1728+0427',
'VLA_J1730+0024':'J1730+0024',
'VLA_J1734+0926':'J1734+0926',
'VLA_J1737+0621':'J1737+0621',
'VLA_J1740+2211':'J1740+2211',
'VLA_J1730+0024':'J1730+0024',
'VLA_J1734+0926':'J1734+0926',
'VLA_J1737+0621':'J1737+0621',
'VLA_J1740+2211':'J1740+2211'}

def sourceRename(sourceList):
    renameList = []
    for srcname in sourceList:
        srcname = srcname.replace('*','')
        renameList = renameList + [sourceDic.get(srcname, srcname)]
    return renameList
#
def AeNominal(msfile, antList):
    msmd.open(msfile)
    antDia = [msmd.antennadiameter(antName)['value'] for antName in antList]
    msmd.close()
    msmd.done()
    return 0.7* 0.25* np.pi* antDia**2      # Nominal Collecting Area
#
def TsysLoad(BandatmSPW, BandName):
    Tau0spec, TrxList, Tau0Coef = [], [], []
    for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
        TrxList = TrxList + [np.median(np.load('%s-%s-SPW%d.Trx.npy' % (prefix, UniqBands[band_index], spw)), axis=3) + np.median(np.load('%s-%s-SPW%d.TantN.npy' % (prefix, UniqBands[band_index], spw)), axis=1)]   # TrxList[spw][pol, ch, ant]
        Tau0spec = Tau0spec + [np.load('%s-%s-SPW%d.Tau0.npy' % (prefix, BandName, spw))]  # Tau0spec[spw][ch]
        Tau0Coef = Tau0Coef + [np.load('%s-%s-SPW%d.Tau0C.npy' % (prefix, BandName, spw))] # Tau0Coef[spw][2] --- intercept + slope
    Tau0E    = np.load(prefix +  '-' + BandName + '.TauE.npy') # Tau0E[spw, atmScan]
    return TrxList, Tau0spec, Tau0Coef, Tau0E
#
def aprioriSEFD(Ae, EL, TrxSpec, Tau0Spec):
    secZ = 1.0 / np.sin(EL)
    zenithTau = Tau0spec + scipy.interpolate.splev(np.median(timeStamp), exTauSP) + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ

    exp_Tau = np.exp(-zenithTau * secZ )
    TsysEQScan = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]

    return 2.0* kb* TsysEQScan.T / Ae
#
#-------- Disk Visibility
def diskVis(diskRadius, u):
    # diskRadius : radius of planet disk [rad]
    # u          : spatial frequency (= baseline / wavelength)
    argument = 2.0* np.pi* u* diskRadius
    return 2.0* scipy.special.jn(1, argument) / argument
#
#-------- Disk Visibility with primary beam correction, u must be smaller than 0.3/diskRadius
def diskVisBeam(diskShape, u, v, primaryBeam):
    # diskShape  : planet disk diameter [MajorAxis, MinorAxis, PA] (rad)
    # u,v        : spatial frequency (= baseline / wavelength)
    # primaryBeam: FWHM of primary beam [rad]
    from interferometry import beamF
    cs, sn = np.cos(diskShape[2]), np.sin(diskShape[2])
    diskRadius = 0.5* np.sqrt(diskShape[0]* diskShape[1])
    DSmaj = 1.0 / np.sqrt( (0.30585 / diskShape[0])**2 + 2.0* np.log(2.0)/(np.pi* primaryBeam)**2 )    # Primary-beam correction
    DSmin = 1.0 / np.sqrt( (0.30585 / diskShape[1])**2 + 2.0* np.log(2.0)/(np.pi* primaryBeam)**2 )    # Primary-beam correction
    uvDisp = (DSmin*(u* cs - v* sn))**2 + (DSmaj*(u* sn + v* cs))**2
    return beamF(diskRadius/primaryBeam)* np.exp(-0.5* uvDisp)
#
#-------- Apertue effciency measurements using Solar System Objects
def SSOAe(antList, antMap, spwDic, uvw, scanDic, SSODic, XSList):
    # antList   : List of antenna name
    # antMap    : Antenna order starting with refant
    # spwDic    : SPW dictionary ['spw', 'freq', 'chNum', 'chRange', 'BW']
    # uvw       : baseline vector [m]
    # scanDic   : scan dictionary ['msfile', 'source', 'mjdSec', 'EL', 'PA', 'I', 'QCpUS', 'Tau', 'Tsys', 'Gain']
    # SSODic    : SSO dictionary
    # XSList    : Cross Correlation XspecList[spw][pol, ch, bl, time]
    from interferometry import GetAntD, Bl2Ant, ANT0, ANT1, kb, gainComplexVec
    SSOname = scanDic['source']
    text_sd = ' Flux Calibrator : %10s EL=%.1f' % (SSOname, 180.0*np.median(scanDic['EL'])/np.pi)
    if np.median(scanDic['EL']) < ELshadow : return              # filter QSO out
    uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    UVlimit = 0.05 / SSODic[SSOname][2][0]          # Maximum usable uv distance [lambda]
    text_sd = ' uv limit = %5.0f klambda' % (UVlimit*1.0e-3); print(text_sd)
    antNum, blNum, antDia = len(antList), len(uvDist), GetAntD(antList)
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    #-------- Check usable antenna/baselines
    uvFlag = np.ones(blNum)
    for spw_index, spw in enumerate(spwDic['spw']):
        chRange = spwDic['chRange'][spw_index]
        centerFreq = np.mean(spwDic['freq'][spw_index][chRange])
        uvFlag[np.where( uvDist > UVlimit* 299792458 / centerFreq)[0].tolist()] *= 0.0
    #
    SAant = [0] + [Bl2Ant(bl_index)[0] for bl_index in np.where(uvFlag > 0.1)[0].tolist() if Bl2Ant(bl_index)[1] == 0]
    SAbl  = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] in SAant and Bl2Ant(bl_index)[1] in SAant]
    if len(SAant) < 4: return #  Too few antennas
    #-------- Baseline map
    print('Subarray : ',); print(antList[np.array(antMap)[SAant]])
    for ant_index in list(range(1, antNum)):
        text_sd = antList[antMap[ant_index]] + ' : '; print(text_sd, end='')
        blList = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] == ant_index]
        for bl_index in blList:
            if uvFlag[bl_index] < 1.0: text_sd = '\033[91m%4.0f\033[0m' % (uvDist[bl_index])
            else: text_sd = '%4.0f' % (uvDist[bl_index])
            print(text_sd, end=' ')
        print('')
    print('       ', end='')
    for ant_index, ant in enumerate(antList[antMap[0:antNum-1]]): print(ant, end=' ')
    print('')
    AeSPW, WgSPW = [], []
    #-------- Aperture Efficiency
    for spw_index, spw in enumerate(spwDic['spw']):
        uvWave = uvw[0:2,:] * centerFreq / 299792458    # UV distance in wavelength
        primaryBeam = 1.13* 299792458 / (np.pi * antDia* centerFreq)
        SSOmodelVis = SSODic[SSOname][1][spw_index]*  diskVisBeam(SSODic[SSOname][2], uvWave[0], uvWave[1], primaryBeam[ant0]* primaryBeam[ant0]* np.sqrt(2.0 / (primaryBeam[ant0]**2 + primaryBeam[ant1]**2)))
        VisChav = np.mean(XSList[spw_index][:,chRange][:,:,SAbl], axis=(1,3)) / SSOmodelVis[SAbl]
        Gain = gainComplexVec(VisChav.T).T  # Gain[pol, ant]
        Aeff = 2.0* kb* abs(Gain)**2 / (0.25* np.pi* antDia[np.array(antMap)[SAant]]**2)
        Ae, Wg = np.zeros([antNum, 2]), np.zeros([antNum, 2])
        for ant_index, SA in enumerate(np.array(antMap)[SAant]):
            Ae[SA] = Aeff[:, ant_index]
            Wg[SA] = np.sign(Aeff[:, ant_index])* np.median(abs(SSOmodelVis))
        #
        AeSPW = AeSPW + [Ae]
        WgSPW = WgSPW + [Wg]
    FscaleDic = {
        'Ae'    : AeSPW,
        'Wg'    : WgSPW}
    return FscaleDic
#-------- Average Ae among multiple SSOs
def averageAe(FscaleDic, antList, spwList):
    AeList, WgList = [], []
    for SSO_index, SSOname in enumerate(FscaleDic.keys()):
        AeSPW = FscaleDic[SSOname]['Ae']; AeList = AeList + [np.array(AeSPW)]
        WgSPW = FscaleDic[SSOname]['Wg']; WgList = WgList + [np.array(WgSPW)]
        text_sd = ' Aeff: '
        for spw_index, spw in enumerate(spwList): text_sd = text_sd + 'SPW%02d-X SPW%02d-Y ' % (spw, spw)
        text_sd = text_sd + '--- %s' % (SSOname)
        print(text_sd)
        for ant_index, ant in enumerate(antList):
            text_sd = '%s : ' % (ant)
            for spw_index, spw in enumerate(spwList):
                for pol_index in [0,1]:
                    if WgSPW[spw_index][ant_index][pol_index] == 0.0:
                        text_sd = text_sd + '  ----- '
                    else:
                        text_sd = text_sd + '  %4.1f%% ' % (100.0* AeSPW[spw_index][ant_index][pol_index])
            print(text_sd)
    #
    return  (np.sum(np.array(WgList) * np.array(AeList), axis=0)/(np.sum(np.array(WgList), axis=0)+1.0e-9)).transpose(1,2,0)   # Aeff[ant, pol, spw]
#-------- Transfer and equalize aperture efficiencies
def AeTransfer(VisChav, Aeff, antDia):
    from interferometry import Bl2Ant, gainComplexVec
    blNum = VisChav.shape[1]
    SAant = np.where(np.min(Aeff, axis=1) > 0.25)[0].tolist()
    SAbl  = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] in SAant and Bl2Ant(bl_index)[1] in SAant]
    GainSA  = gainComplexVec(VisChav[:,SAbl].T)
    GainAll = gainComplexVec(VisChav.T)
    scaleFlux = np.median(abs(GainSA**2).T / (antDia**2* Aeff.T), axis=1)
    return ((abs(GainAll)**2 / scaleFlux).T / (antDia**2)).T
#
#-------- Gain scaling
def GainScale(Aeff, antDia, polVis):
    # Aeff :    Aperure Efficiency Table
    # antDia :  Antenna diameter
    # polVis[pol,bl,time] :    cross power spectra
    from interferometry import kb, ANT0, ANT1
    polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
    blNum = polVis.shape[1]
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    fluxScale = np.sqrt( 2.0* kb / (0.25* np.pi* antDia**2* Aeff.T))
    return (polVis.transpose(2,0,1)* fluxScale[polYindex][:,ant0]* fluxScale[polXindex][:,ant1]).transpose(1,2,0)
#-------- Stokes parameters
def Vis2Stokes(VisChav, Dcat, PA):
    # VisChav   : channel-averaged visibiliities [pol, bl, time]
    # Dcat      : D-term [ant, pol]
    # PA        : Parallactic Angle with respect to X-feed [time]
    from interferometry import ANT0, ANT1, InvMullerMatrix, InvPAVector
    PAnum = len(PA)
    PS = InvPAVector(PA, np.ones(PAnum))
    blNum = VisChav.shape[1]
    Stokes = np.zeros([4,blNum], dtype=complex)
    for bl_index in list(range(blNum)):
        Minv = InvMullerMatrix(Dcat[ANT1[bl_index], 0], Dcat[ANT1[bl_index], 1], Dcat[ANT0[bl_index], 0], Dcat[ANT0[bl_index], 1])
        Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.dot( VisChav[:,bl_index]).reshape(4*PAnum)) / PAnum
    #
    return Stokes
#-------- Smooth time-variable Tau
def tauSMTH( timeSample, TauE ):
    if len(timeSample) > 5:
        SplineWeight = np.ones(len(timeSample) + 4)
        flagIndex = (np.where(abs(TauE - np.median(TauE))/np.std(TauE) > 3.0)[0] + 2).tolist()
        SplineWeight[flagIndex] = 0.01
        tempTime = np.append([timeSample[0]-500.0, timeSample[0]-300.0], np.append(timeSample, [timeSample[-1]+300.0, timeSample[-1]+500.0]))
        tempTauE = np.append([TauE[0], TauE[0]], np.append(TauE, [TauE[-1], TauE[-1]]))
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3, w=SplineWeight, t=tempTime[list(range(1, len(tempTime), 2))] - 60.0 )
    else:
        tempTime = np.arange(np.min(timeSample) - 3600.0,  np.max(timeSample) + 3600.0, 300.0)
        tempTauE = np.repeat(np.median(TauE), len(tempTime))
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3)
    #
    return smthTau
#
'''
#-------- Load Tsys table
Tau0spec, TrxList, Tau0Coef = [], [], []
for spw in spwList:
    TrxList = TrxList + [np.median(np.load('%s-%s-SPW%d.Trx.npy' % (prefix, UniqBands[band_index], spw)), axis=3) + np.median(np.load('%s-%s-SPW%d.TantN.npy' % (prefix, UniqBands[band_index], spw)), axis=1)]   # TrxList[spw][pol, ch, ant]
    Tau0spec = Tau0spec + [np.load('%s-%s-SPW%d.Tau0.npy' % (prefix, UniqBands[band_index], spw))]  # Tau0spec[spw][ch]
    Tau0Coef = Tau0Coef + [np.load('%s-%s-SPW%d.Tau0C.npy' % (prefix, UniqBands[band_index], spw))] # Tau0Coef[spw][2] --- intercept + slope
#
TrxAnts  = np.load(prefix +  '-' + UniqBands[band_index] + '.TrxAnt.npy') # TrxAnts[ant]
Tau0E    = np.load(prefix +  '-' + UniqBands[band_index] + '.TauE.npy') # Tau0E[spw, atmScan]
atmTimeRef = np.load(prefix +  '-' + UniqBands[band_index] + '.atmTime.npy') # atmTimeRef[atmScan]
TrxMap = indexList(TrxAnts, antList); TrxFlag = np.zeros([antNum]); TrxFlag[TrxMap] = 1.0
Tau0E = np.nanmedian(Tau0E, axis=0); Tau0E[np.isnan(Tau0E)] = np.nanmedian(Tau0E); Tau0E[np.isnan(Tau0E)] = 0.0
for spw_index in list(range(spwNum)):
    TrxMed = np.median(TrxList[spw_index], axis=1)  # TrxMed[pol, ant]
    for pol_index in [0,1]:
        Trx2MedianRatio = TrxMed[pol_index] / np.median(TrxMed[pol_index])
        TrxFlag[ np.where(Trx2MedianRatio < 0.1)[0].tolist() ] *= 0.0   # Flagged by negative Trx
        TrxFlag[ np.where(Trx2MedianRatio > 3.0)[0].tolist() ] *= 0.0   # Flagged by too-high Trx 
    if np.median(Tau0spec[spw_index][chRange]) < 0.0: TrxFlag *= 0.0    # Negative Tau(zenith) 
#
#TrxFlag[useAnt] *= np.median(np.min(scanFlag, axis=1), axis=(0,2))
print('Ant: ', end='')
for ant_index in list(range(antNum)): print(antList[ant_index], end=' ')
print(); print('givn', end='')
for ant_index in list(range(antNum)): print('    %.0f' % (flagAnt[ant_index]), end='')
print(); print('Trx ', end='')
for ant_index in list(range(antNum)): print('    %.0f' % (TrxFlag[ant_index]), end='')
print(); print('gain', end='')
for ant_index in list(range(antNum)): print('    %.0f' % (gainFlag[ant_index]), end='')
print()
flagAnt = flagAnt* TrxFlag* gainFlag
UseAnt = np.where(flagAnt > 0.0)[0]; UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
print('%d / %d usable antennas' % (UseAntNum, antNum))
if len(UseAnt) < 4: sys.exit('Too few usable antennas. Reduction failed.')
#-------- Check Scans for atmCal
ingestFile = open(prefix + '-' + UniqBands[band_index] + '-Ingest.log', 'w') 
text_sd = '#source,    RA,eRA,dec,edec, frequency,  flux,  eflux,     %P,    d%P,  EVPA, eEVPA, uvmin, uvmax,         date, fluxCal, ExecBlock\n'; ingestFile.write(text_sd)
logfile = open(prefix + '-' + UniqBands[band_index] + '-Flux.log', 'w') 
logfile.write(BPcalText + '\n'); logfile.write(EQcalText + '\n')
#-------- Tsys measurements
scanList = onsourceScans
msmd.close()
msmd.done()
#-------- Array Configuration
'''
'''
print('---Determining refant')
if 'refant' in locals(): refantID = np.where(antList == refant)[0][0]
if 'refantID' not in locals():
    msmd.open(msfile)
    timeStamp, UVW = GetUVW(msfile, spwList[0], msmd.scansforspw(spwList[0])[0])
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = UseAnt[bestRefant(uvDist, UseAnt)]
print('  Use ' + antList[refantID] + ' as the refant.')
#
antMap = [refantID] + UseAnt[np.where(UseAnt != refantID)].tolist()
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
'''
'''
print('  %d baselines are inverted' % (len(np.where( blInv )[0])))
AeNominal = 0.7* 0.25* np.pi* antDia**2      # Nominal Collecting Area
#-------- Flag table
if 'FGprefix' in locals():
    print('---Checking Flag File')
    FGList = []
    for spw_index in list(range(spwNum)): FG = np.load('%s-SPW%d.FG.npy' % (FGprefix, spwList[spw_index])); FGList = FGList + [np.min(FG, axis=0)]
    FG = np.min( np.array(FGList), axis=0)
    TS = np.load('%s-SPW%d.TS.npy' % (FGprefix, spwList[spw_index]))
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = np.where(FG[indexList(timeStamp, TS)] == 1.0)[0]
else :
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
    flagIndex = list(range(timeNum))
#
#-------- Bandpass Table
Trx2antMap = indexList( antList[antMap], antList[TrxMap] )
BPDone = False
print('---Generating antenna-based bandpass table')
BPList = []
secZ = 1.0 / np.mean(np.sin(BPEL))
for spw_index in list(range(spwNum)):
    #BP_ant, XY_BP, XYdelay, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPScan, blMap, blInv)     # BP_ant[antMap, pol, ch]
    #print('BPscan = %d' % BPscan)
    BP_ant = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, BPscan, spwList[spw_index]))
    XY_BP = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, BPscan, spwList[spw_index]))
    BP_ant[:,1] *= XY_BP
    zenithTau = Tau0spec[spw_index] + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ
    exp_Tau = np.exp(-zenithTau * secZ )
    atmCorrect = 1.0 / exp_Tau
    TsysBPScan = atmCorrect* (TrxList[spw_index].transpose(2,0,1)[Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)) # [antMap, pol, ch]
    TsysBPShape = abs((TsysBPScan.transpose(2,0,1) / np.median(TsysBPScan, axis=2))).transpose(1,2,0)
    BPList = BPList + [BP_ant* np.sqrt(TsysBPShape)]
#
if PLOTBP:
    pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, antList[refantID], BPscan))
    plotBP(pp, prefix, antList[antMap], spwList, BPscan, BPList)
#
BPDone = True
##-------- Equalization using EQ scan
scanList = onsourceScans
GainP, AeSeqX, AeSeqY = [], [], []  # effective area x flux density of the equalizer
polXindex, polYindex, scan_index = (arange(4)//2).tolist(), (arange(4)%2).tolist(), onsourceScans.index(EQScan)
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], EQScan); timeNum = len(timeStamp)
AzScan, ElScan = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refantID, AZ, EL)
PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; PA = np.arctan2( np.sin(PA), np.cos(PA))
if len(StokesDic[EQcal]) > 0:
    StokesEQ = np.array(StokesDic[EQcal])
else: 
    StokesEQ = np.array([1.0, 0.0, 0.0, 0.0])
#
QCpUS = (StokesEQ[1]* np.cos(2.0* PA) + StokesEQ[2]* np.sin(2.0* PA)) / StokesEQ[0]
#
##-------- Smooth excess Tau 
exTauSP = tauSMTH(atmTimeRef, Tau0E)
secZ = 1.0 / np.mean(np.sin(ElScan))
##-------- Scale aperture efficiency
for spw_index in list(range(spwNum)):
    zenithTau = Tau0spec[spw_index] + scipy.interpolate.splev(np.median(timeStamp), exTauSP) + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ
    exp_Tau = np.exp(-zenithTau * secZ )
    #TsysEQScan = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
    atmCorrect = 1.0 / exp_Tau
    TsysEQScan = np.median(atmCorrect* (TrxList[spw_index].transpose(2,0,1)[Trx2antMap] + Tcmb*exp_Tau + tempAtm* (1.0 - exp_Tau)), axis=2) # [antMap, pol]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], EQScan)
    chNum = Xspec.shape[1]; chRange = list(range(int(0.07*chNum), int(0.95*chNum)))
    tempSpec = CrossPolBL(Xspec[:,:,blMap], blInv).transpose(3,2,0,1)[flagIndex]       # Cross Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    BPCaledXspec = (tempSpec / (BPList[spw_index][ant0][:,polYindex]* BPList[spw_index][ant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)[pPol]  # chAvgVis[pPol,bl,time]
    GainP = GainP + [np.array([np.apply_along_axis(clphase_solve, 0, chAvgVis[0]), np.apply_along_axis(clphase_solve, 0, chAvgVis[1])])]
    pCalVisX = np.mean(chAvgVis[0] / (1.0 + QCpUS) / (GainP[spw_index][0,ant0]* GainP[spw_index][0,ant1].conjugate()), axis=1)
    pCalVisY = np.mean(chAvgVis[1] / (1.0 - QCpUS) / (GainP[spw_index][1,ant0]* GainP[spw_index][1,ant1].conjugate()), axis=1)
    #-------- Antenna-based Gain
    AeSeqX = AeSeqX + [2.0* kb* TsysEQScan[:, 0]* abs(gainComplex(pCalVisX))**2] # Ae x Seq (X-pol)
    AeSeqY = AeSeqY + [2.0* kb* TsysEQScan[:, 1]* abs(gainComplex(pCalVisY))**2] # Ae x Seq (Y-pol)
#
GainP = np.array(GainP) # GainP[spw, pol, ant, time]
AeSeqX, AeSeqY = np.array(AeSeqX), np.array(AeSeqY) # AeSeqX[spw, antMap], AeSeqX[spw, antMap] : (relative) aperture efficiency, assuming 1 Jy
##-------- inter-SPW phasing using EQ scan
spwPhase = [0.0]* 2* spwNum
for ant_index in list(range(1,UseAntNum)):
    for pol_index in [0, 1]:
        spwPhase = spwPhase + [0.0]
        for spw_index in list(range(1,spwNum)): spwPhase = spwPhase + [np.angle(GainP[spw_index, pol_index, ant_index].dot(GainP[0, pol_index, ant_index].conjugate()))]
    #
#
spwPhase = np.array(spwPhase).reshape([UseAntNum, 2, spwNum]); spwTwiddle = exp(1.0j *spwPhase)
for spw_index in list(range(1,spwNum)): BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
#-------- Flux models for solar system objects
msmd.done()
exec(open(SCR_DIR + 'SSOflux.py').read()); logfile.write(FLScaleText + '\n')
######## Outputs from SSOflux.py :
#  SSOflux0[SSO, spw] : model flux density of solar system objects at zero spacing
#  SSOmodelVis[SSO, spw, bl] : predicted relative visibility (max=1.0)
#  uvFlag[SSO, spw, bl] : 0=resolved, 1=unresolved
########
flaggedBlList = list(set(range(blNum)) - set(blMap)); uvFlag[:,:,flaggedBlList] = 0.0
atmCorrect = []
for spw_index in list(range(spwNum)): atmCorrect = atmCorrect + [np.median(Tau0spec[spw_index])]
atmCorrect = np.exp(-outer( np.array(atmCorrect),  1.0/np.sin( np.array(OnEL)[indexList(np.array(SSOscanID), np.array(onsourceScans))]))).T
SSOflux = SSOflux0* atmCorrect  # SSOflux[SSO, spw] : attenuated SSO flux
uvFlag = np.min(uvFlag, axis=1) # all-SPW uv flag
##-------- Scaling with the flux calibrator
AeX, AeY = np.zeros([UseAntNum, spwNum, SSONum]), np.zeros([UseAntNum, spwNum, SSONum])
#-------- Sub-array with unflagged antennas (short baselines)
SSO_flag = np.ones(SSONum)
for sso_index in list(range(SSONum)):
    scan_index = indexList(np.array(SSOscanID), np.array(onsourceScans))[sso_index]
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], SSOscanID[sso_index]); timeNum = len(timeStamp)
    AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, refantID, AZ, EL)
    if (not BLCORR) and np.max(ElScan) < 30.0/180.0*pi : continue    # Too low Elevation
    SAantMap, SAblMap, SAblInv = subArrayIndex(uvFlag[sso_index], refantID) # in Canonical ordering
    if len(SAantMap) < 4: continue #  Too few antennas
    print('Subarray : ',); print(antList[SAantMap])
    SAblIndex = indexList(np.array(SAblMap), np.array(blMap))
    SAant0, SAant1 = np.array(ant0)[SAblIndex].tolist(), np.array(ant1)[SAblIndex].tolist()
    bpAntMap = indexList(antList[SAantMap],antList[antMap])
    Trx2antMap = indexList( antList[SAantMap], antList[TrxMap] )
    secZ = 1.0 / np.mean(np.sin(ElScan))
    for spw_index in list(range(spwNum)):
        #-------- Baseline-based cross power spectra
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], SSOscanID[sso_index])
        chNum = Xspec.shape[1]; chRange = list(range(int(0.07*chNum), int(0.95*chNum)))
        tempSpec = CrossPolBL(Xspec[:,:,SAblMap], SAblInv).transpose(3,2,0,1)      # Cross Polarization Baseline Mapping
        BPCaledXspec = (tempSpec / (BPList[spw_index][SAant0][:,polYindex]* BPList[spw_index][SAant1][:,polXindex].conjugate())).transpose(2,3,1,0) # Bandpass Cal
        chAvgVis = (np.mean(BPCaledXspec[:, chRange], axis=1)[pPol].transpose(0,2,1)/SSOmodelVis[sso_index,spw_index,SAblMap]).transpose(0,2,1)
        GainX, GainY = np.apply_along_axis( gainComplex, 0, chAvgVis[0]), np.apply_along_axis( gainComplex, 0, chAvgVis[1])
        #-------- Tsys
        Ta = SSOflux[sso_index, spw_index]* AeNominal[SAantMap] / (2.0* kb)
        zenithTau = Tau0spec[spw_index] + scipy.interpolate.splev(np.median(timeStamp), exTauSP) + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ
        exp_Tau = np.exp(-zenithTau * secZ )
        TsysSPW = np.median(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
        #-------- Aperture efficiency
        AeX[bpAntMap, spw_index, sso_index] = 2.0* kb* np.median(abs(GainX), axis=1)**2 * (Ta + TsysSPW[:, 0]) / SSOflux[sso_index, spw_index]
        AeY[bpAntMap, spw_index, sso_index] = 2.0* kb* np.median(abs(GainY), axis=1)**2 * (Ta + TsysSPW[:, 1]) / SSOflux[sso_index, spw_index]
    #
#
#-------- SSO Flagging
for sso_index in list(range(SSONum)):
    for spw_index in list(range(spwNum)):
        index = np.where(AeX[:, spw_index, sso_index] > 1.0)[0].tolist()
        if len(index) < 4: SSO_flag[sso_index] = 0.0; continue
        FLX_stat, FLY_stat = AeX[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()], AeY[index, spw_index, sso_index]/AeNominal[np.array(antMap)[index].tolist()]
        if np.median(FLX_stat) < 0.4: SSO_flag[sso_index] = 0.0 ; print('FLX < 0.4')
        if np.median(FLY_stat) < 0.4: SSO_flag[sso_index] = 0.0 ; print('FLY < 0.4')
        if np.median(FLX_stat) > 2.5: SSO_flag[sso_index] = 0.0 ; print('FLX > 2.5')
        if np.median(FLY_stat) > 2.5: SSO_flag[sso_index] = 0.0 ; print('FLY > 2.5')
        if np.percentile(FLX_stat, 75) / np.median(FLX_stat) > 2.0: SSO_flag[sso_index] = 0.0   ; print('FLX 75%-percentile / median > 2.0')
        if np.percentile(FLY_stat, 75) / np.median(FLY_stat) > 2.0: SSO_flag[sso_index] = 0.0   ; print('FLY 75%-percentile / median > 2.0')
    #
    try:
        if SSO_flag[sso_index] == 0.0: raise END
    except END:
        text_sd = '%s is not available as a flux calibrator.' % (sourceList[BandSSOList[sso_index]])
        logfile.write(text_sd + '\n'); print(text_sd)
    #
#
SSOUseList = np.where(SSO_flag == 1.0)[0].tolist()
fluxCalText = ''
if len(SSOUseList) > 0:
    fluxCalText = sourceList[BandSSOList[SSOUseList[0]]]
    for sso_ID in SSOUseList:
        fluxCalText = fluxCalText + sourceList[BandSSOList[sso_ID]] + '-'
    #
    fluxCalText = fluxCalText[:-1]
    exec(open(SCR_DIR + 'AmpCalStokes.py').read())
else:
    fluxCalText = 'SEFD'
    print('No available Solar System Objects!! Try a-priori calibration.')
    exec(open(SCR_DIR + 'aprioriStokes.py').read())
#
msmd.close()
msmd.done()
'''
