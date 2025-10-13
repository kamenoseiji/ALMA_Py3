import numpy as np
import sys
import subprocess
from scipy import stats
from ASDM_XML import spwMS, SPW_LO
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List  e.g. 3,6', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
scanList  = [] if options.scanList == '' else [int(scan) for scan in options.scanList.split(',')]
'''
prefix = 'uid___A002_X12dca4b_X1b2'
scanList = [3]
'''
spurLog = open(prefix + '-LO2Spur.log', 'w')
#-------- Get LO1 and LO2 frequencies
msfile = prefix + '.ms'
if not os.path.isdir(msfile): importasdm(prefix)
SPWdic = spwMS(msfile)
SPWLO = SPW_LO(SPWdic, prefix)
BBList = np.sort(np.unique(np.array([SPWdic[spw]['BB'] for spw in SPWdic.keys()]))).tolist()
#-------- Reconfigure SPWs in MS
#
text_sd = 'SPW Band  BB   LO1 [GHz]  LO2 [GHz]  RF_range [GHz]     IF1_range [GHz]    Leakages'; print(text_sd) ; spurLog.write(text_sd + '\n') 
text_sd = '-------------------------------------------------------------------------------------'; print(text_sd) ; spurLog.write(text_sd + '\n') 
for spw_index, spw in enumerate(SPWdic.keys()):
    if set(scanList) & set(SPWdic[spw]['scanList']) == set(): continue 
    LO1, RFrange = SPWdic[spw]['LO1']*1.0e-9, np.sort(np.array([SPWdic[spw]['ch0'], SPWdic[spw]['ch0'] + (SPWdic[spw]['chNum']-1)* SPWdic[spw]['chStep']]))* 1.0e-9
    IFrange = np.sort(abs(RFrange - LO1))
    text_sd = '%02d  %s BB%d  %10.6f %9.6f [%7.3f - %7.3f] [%6.3f - %6.3f]' % (spw, SPWdic[spw]['Band'], SPWdic[spw]['BB']+1, LO1, SPWdic[spw]['LO2']*1.0e-9, RFrange[0], RFrange[1], IFrange[0], IFrange[1])
    for BB in BBList:
        LO2 = [SPWdic[spw2]['LO2']*1.0e-9 for spw2 in SPWdic.keys() if SPWdic[spw2]['BB'] == BB][0]
        if (IFrange[0] < LO2) and (LO2 < IFrange[1]): text_sd = text_sd + ' RF=%.3f ' % (LO1 + SPWdic[spw]['SB1']*LO2)
    print(text_sd) ; spurLog.write(text_sd + '\n')
spurLog.close()
