import numpy as np
import matplotlib.pyplot as plt
from interferometry import GetAntName, GetAntPos, GetAzEl, GetAzOffset, GetChNum, GetPolQuery, GetSourceDic, AzElMatch, indexList, interValue, smoothGain, QUscale, AzEl2PA, simpleGaussFit
SCR_DIR = os.environ['HOME'] + '/ALMA_Py3/'
RADSEC = np.pi/180/3600 # radian per arcsec
CVEL   = 299792458      # spacetime constant
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-p', dest='scanAnt', metavar='scanAnt',
    help='Antenna to plot', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
scanAnt = options.scanAnt
spwList = [] if options.spwList == '' else [int(spw) for spw in options.spwList.split(',')]
'''
prefix = 'uid___A002_X1366160_X4514.WVR'
spwList = [0,1,2,3]
scanAnt = 'DV01'
'''
#-------- Check array configuration
msfile = prefix + '.ms'
mjdSec, AntID, AZ, EL = GetAzEl(msfile)
mjdSec, AntID, Offset = GetAzOffset(msfile)
antList = GetAntName(msfile)
antPos = GetAntPos(msfile)
antNum = len(antList)
stareAntList = [antList[ant_index] for ant_index in list(range(antNum)) if np.std( Offset[0,np.where(AntID == ant_index)[0].tolist()] ) < 1.0]
scanAntList = list( set(antList) - set(stareAntList) )
scanAntID = np.where(antList == scanAnt)[0][0]
if scanAnt not in scanAntList: 
    print('%s is not in scan antenna list : ' % (scanAnt) ); print(scanAntList)
    a=1/0
stareAntIDs = indexList(stareAntList, antList)
refantID = stareAntIDs[np.argmin([(antPos[:,stare_index] - antPos[:,scanAntID]).dot(antPos[:,stare_index] - antPos[:,scanAntID]) for stare_index in stareAntIDs])]
refant = antList[refantID]
#---- AZEL of the refant
time_index = np.where(AntID == refantID)[0].tolist() 
az,el = AZ[time_index], EL[time_index]
#---- Offset of the scan ant
time_index = np.where(AntID == scanAntID)[0].tolist()
mjdSec, Offset = mjdSec[time_index], Offset[:,time_index]
flagIndex = np.where( Offset[0]**2 + Offset[1]**2 < 1.5* np.percentile( Offset[0]**2 + Offset[1]**2, 99 ) )[0]  # remove off position
mjdSec, Offset, az, el = mjdSec[flagIndex], Offset[:,flagIndex], az[flagIndex], el[flagIndex]
srcDic = GetSourceDic(msfile)
GAI = []
for spw in spwList:
    chNum, chWid, freq = GetChNum(msfile, spw)
    #-------- Gain Cal in scan 1
    text_sd = 'casa -c %scheckGain.py -u %s -s %d -c %d -R %s' % (SCR_DIR, prefix, spw, 1, refant) 
    print(text_sd)
    os.system(text_sd)
    antList = np.load('%s.Ant.npy' % (prefix))
    ant_index = np.where(antList == scanAnt)[0][0]
    GA = np.load('%s-SPW%d.GA.npy' % (prefix, spw))[ant_index]
    TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
    azTS, elTS = interValue(mjdSec, az, TS), interValue(mjdSec, el, TS)
    paTS = AzEl2PA(azTS, elTS)
    for sourceID in srcDic.keys():
        IQU = GetPolQuery(srcDic[sourceID]['Name'], TS[0], np.median(freq)*1.0e-9, SCR_DIR)
        srcDic[sourceID]['I'], srcDic[sourceID]['Q'], srcDic[sourceID]['U'] = IQU[0][srcDic[sourceID]['Name']], IQU[1][srcDic[sourceID]['Name']], IQU[2][srcDic[sourceID]['Name']]
    Xscale, Yscale = QUscale(paTS, srcDic[0]['Q']/srcDic[0]['I'], srcDic[0]['U']/srcDic[0]['I'])
    GAXY = np.mean( GA / np.array([np.sqrt(Xscale), np.sqrt(Yscale)]), axis=1)
    #-------- Gain solutions in scan 2
    text_sd = 'casa -c %scheckGain.py -u %s -s %d -c %d -R %s -a ' % (SCR_DIR, prefix, spw, 2, refant) 
    for ant in scanAntList:
        if ant != scanAnt: text_sd = text_sd + ant + ','
    text_sd = text_sd[:-1]
    print(text_sd)
    os.system(text_sd)
    #-------- Antenna gain for the scanning antenna
    antList = np.load('%s.Ant.npy' % (prefix))
    ant_index = np.where(antList == scanAnt)[0][0]
    GA = np.load('%s-SPW%d.GA.npy' % (prefix, spw))[ant_index]
    TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
    GAI = GAI + [np.mean((GA.T/GAXY), axis=1)]
#-------- Dish Mapping : (u, v) stand for the position in the dish surface [m]
azoTS, eloTS = interValue(mjdSec, Offset[0], TS), interValue(mjdSec, Offset[1], TS) 
dishPlace = np.arange(-6,6,0.1)
PI= np.zeros([len(dishPlace), len(dishPlace)], dtype=complex)
for spw_index, spw in enumerate(spwList):
    chNum, chWid, freq = GetChNum(msfile, spw)
    phaseScale = 2.0* np.pi* freq[0]*RADSEC/CVEL
    surfaceScale = CVEL/(4.0* np.pi*freq[0])*1.0e6
    ampScale = 1.0/np.sum(np.abs(GAI))
    for u_index,u in enumerate(dishPlace):
        for v_index,v in enumerate(dishPlace):
            twiddle = np.exp(phaseScale* (0.0+1.0j)* (azoTS*u + eloTS*v))
            PI[u_index,v_index] += GAI[spw_index].dot(twiddle)
PI *= ampScale
figFileName = '%s-%s.amp.png' % (prefix, scanAnt)
plt.imshow(abs(PI), extent=(-6,6,-6,6))
plt.colorbar(); plt.title(prefix + ' ' + scanAnt)
plt.savefig(figFileName, format='png', dpi=72)
plt.close('all')
figFileName = '%s-%s.phase.png' % (prefix, scanAnt)
plt.imshow(np.angle(PI)*surfaceScale, extent=(-6,6,-6,6))
plt.colorbar(); plt.title(prefix + ' ' + scanAnt)
plt.savefig(figFileName, format='png', dpi=72)
plt.close('all')
np.save('AH-%s-%s.npy' % (prefix, ant), PI)
