#-------- Plot Multi-BB Group Delay
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from casatools import quanta as qatool
import analysisUtils as au
import datetime
from interferometry import indexList, GetAntName, GetAntPos
from Plotters import polColor, polName
qa = qatool()
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-M', dest='plotRange', metavar='plotRange',
    help='Max y-axis range in [ps] e.g. 10.0', default='10.0')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
plotMax = float(options.plotRange)
'''
prefix = 'uid___A002_X12e95ca_X7383'
plotMax = 10.0
'''
DL = np.load(prefix + '-GDL.npy')
msfile = prefix + '.ms'
antList = GetAntName(msfile) 
antPos  = GetAntPos(msfile)
DLantList = np.load(prefix + '.Ant.npy')
antNum = len(antList)
antMap = indexList(DLantList, antList)
plotAntList = antList[antMap].tolist()
refVector = antPos[:,antMap[0]]
R = np.sqrt(refVector.dot(refVector))
refVector_n = refVector/R
R_xy = np.sqrt(refVector_n[:2].dot(refVector_n[:2]))
CS, SN = refVector_n[0]/R_xy, refVector_n[1]/R_xy
Rm = np.array([ [refVector_n[2]*CS, refVector_n[2]*SN, -R_xy], [-SN,  CS, 0.0], [R_xy*CS, R_xy*SN, refVector_n[2]]])
ENpos = np.array([[0.0, 1.0], [-1.0, 0.0]]).dot(Rm.dot(antPos)[:2])
#-------- Load tables
DT, TS = [], DL[0]
for mjdSec in TS.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
pp = PdfPages('DLMAP_%s.pdf' %  (prefix))
figDL = plt.figure(figsize = (8, 11))
figDL.suptitle(prefix + ' Multi-BB Group Residual Delay Map')
for scan_index, UTC in enumerate(DT):
    DLPL = figDL.add_subplot(1,1,1)
    DLPL.yaxis.offsetText.set_fontsize(3)
    plotDL = 0.5e3* (DL[1+2*antNum:1+3*antNum, scan_index] + DL[1+3*antNum:1+4*antNum, scan_index])
    mappable = DLPL.scatter(ENpos[0,antMap], ENpos[1,antMap], c=plotDL, cmap='gist_rainbow', vmin=-plotMax, vmax=plotMax)
    DLPL.set_aspect('equal', adjustable='box')
    DLPL.set_title('%s %s' % (prefix, UTC.strftime('%Y-%m-%dT%H:%M:%S')))
    figDL.colorbar(mappable, ax=DLPL, orientation='horizontal')
    figDL.savefig(pp, format='pdf')
    figDL.delaxes(DLPL)
plt.close('all')
pp.close()

