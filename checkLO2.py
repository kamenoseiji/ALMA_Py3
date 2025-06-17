import sys
import subprocess
from scipy import stats
from ASDM_XML import SPW_FULL_RES
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spurLog = open(prefix + '-LO2Spur.log', 'w')
#-------- Get LO1 and LO2 frequencies
SPWdic = SPW_FULL_RES(prefix)
text_sd = 'SPW Band  BB   LO1 [GHz]  LO2 [GHz]  RF_range [GHz] '; print(text_sd) ; spurLog.write(text_sd + '\n') 
text_sd = '----------------------------------------------------'; print(text_sd) ; spurLog.write(text_sd + '\n') 
for spw_index, spw in enumerate(SPWdic.keys()):
    LO1, RF_S, RF_E = SPWdic[spw]['LO1'], SPWdic[spw]['ch0']*1.0e-9, (SPWdic[spw]['ch0'] + SPWdic[spw]['SB1']*SPWdic[spw]['SB2']*SPWdic[spw]['SB3']*SPWdic[spw]['chNum']* SPWdic[spw]['chWidth'])*1.0e-9
    text_sd = '%02d  %s BB%d  %10.6f %9.6f [%7.3f - %7.3f] ' % (spw, SPWdic[spw]['Band'], SPWdic[spw]['BB']+1, SPWdic[spw]['LO1']* 1.0e-9, SPWdic[spw]['LO2']*1.0e-9, RF_S, RF_E)
    #-------- Possible LO2 leakage
    for anotherSPW in SPWdic.keys():
        if LO1 != SPWdic[anotherSPW]['LO1']: continue
        LO2inRF = SPWdic[anotherSPW]['SB2']* SPWdic[anotherSPW]['SB3']* SPWdic[anotherSPW]['LO1'] - SPWdic[anotherSPW]['SB1']* SPWdic[anotherSPW]['SB3']* SPWdic[anotherSPW]['LO2']
        if (LO2inRF - RF_S)* (LO2inRF - RF_E) < 0: text_sd = text_sd + ' Spurious at RF=%10.6f GHz' % (LO2inRF)
#
    print(text_sd) ; spurLog.write(text_sd + '\n')
spurLog.close()
