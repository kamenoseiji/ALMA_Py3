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
    #-------- antennas to flag
    flagAntList = antList[np.unique(np.where(FG < 1.0)[0])].tolist()
    for ant in flagAntList:
        ant_index = np.where(antList == ant)[0][0]
        text_sd = 'SPW%d %s' % (spw, ant)
        flaggedTimeIndex = np.where(FG[ant_index] < 1.0)[0]
        if len(np.where(np.diff(flaggedTimeIndex) == 1)[0]) > 0.01* len(TS) : # more than 1% of continuously flagged time
            print(text_sd)
            logfile.write(text_sd + '\n')
#
logfile.close()
plotFlag(pp, prefix, antList, DT, spwList, FGList)
