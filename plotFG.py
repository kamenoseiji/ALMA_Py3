import os
import numpy as np
import datetime
import analysisUtils as au
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Plotters import plotFlag
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = [int(spw) for spw in options.spwList.split(',')]
'''
prefix = '2023.1.00533.S_X11fe144_X23dc1'
spwList = [0]
'''
#-------- Load Flag data
ANfileName = '%s.Ant.npy' % (prefix)
if not os.path.isfile(ANfileName):
    print('%s does not exist.' % (ANfileName))
    exit()
antList = np.load(ANfileName)
TS = np.load('%s-SPW%d.TS.npy' % (prefix, spwList[0]))
DT = []
for mjdSec in TS.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))

pp = PdfPages('FG_%s.pdf' % (prefix))
logfile = open(prefix + '-PolFlag.log', 'w')
FGList, flagAntList = [], []
for spw_index, spw in enumerate(spwList):
    FGfileName = '%s-SPW%d.FG.npy' % (prefix, spw)  # Flag file [ant, time]
    if not os.path.isfile(FGfileName):
        print('%s does not exist.' % (FGfileName))
        continue
    FG = np.load(FGfileName)
    FGList = FGList + [FG]
    antNum, timeNum = FG.shape[0], FG.shape[1]
    #-------- antennas to flag
    flagAntIndex = np.where( np.quantile(FG, 0.25, axis=1) < 1)[0].tolist()
    flagAntList = antList[flagAntIndex]
    unFlaggedIndex = list(set(range(antNum)) - set(flagAntIndex))
    unFlaggedAntNum = len(unFlaggedIndex)
    for ant_index, ant in enumerate(flagAntList):
        text_sd = 'SPW%d %s' % (spw, ant)
        print(text_sd)
        logfile.write(text_sd + '\n')
    #-------- Time to flag
    flagTimeList = np.where( np.quantile(FG[unFlaggedIndex], 3.0/unFlaggedAntNum, axis=0) < 1)[0].tolist()
    for time_index in flagTimeList:
        text_sd = 'SPW%d %.3f' % (spw, TS[time_index])
        print(text_sd)
        logfile.write(text_sd + '\n')
    #
#
logfile.close()
plotFlag(pp, prefix, antList, DT, spwList, FGList)
