# TsysScanCopy.py : copy Tsys values from the source scan to the destinate scan
import numpy as np
from interferometry import GetAntName
polDic = {'X':0, 'Y':1, 'x':0, 'y':1}
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-f', dest='fileName', metavar='fileName',
    help='Tsys file name  e.g. uid___A002_X12c99be_X55f7T0.tsys', default='')
parser.add_option('-a', dest='antenna', metavar='antenna',
    help='Antenna origin e.g. DA41', default='')
parser.add_option('-s', dest='spws', metavar='spws',
    help='SPWs to copy e.g. 0,1,2', default='')
parser.add_option('-c', dest='scans', metavar='scans',
    help='Scans (from, to) e.g. 2,5 indicates copy scan 2 to 5', default='')
parser.add_option('-p', dest='polList', metavar='polList',
    help='Polarization List e.g. X,Y', default='')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
fileName  = options.fileName
antName   = options.antenna
scans  = [] if options.scans == '' else [int(scan) for scan in options.scans.split(',')]
spws   = [] if options.spws == '' else [int(spw) for spw in options.spws.split(',')]
polList= [] if options.polList == '' else [pol for pol in options.polList.split(',')]
'''
prefix = 'uid___A002_X12cdde9_X59a6'
fileName = 'uid___A002_X12cdde9_X59a6T0.tsys'
antName = 'DA44'
scans = [4,2]
spws = [17,19,21,23]
polList = ['X','Y']
'''
#-------- Open Tsys table
antList = GetAntName(prefix + '.ms'); antNum = len(np.unique(antList))
ant_index = np.where(antList == antName)[0][0]
tb.open(fileName, nomodify=False)
AntID = tb.getcol('ANTENNA1'); antNum = len(np.unique(AntID))
ScanID = tb.getcol('SCAN_NUMBER'); scanNum = len(np.unique(ScanID))
SpwID  = tb.getcol('SPECTRAL_WINDOW_ID'); spwNum = len(np.unique(SpwID))
Tsys = tb.getcol('FPARAM'); polNum, chNum = Tsys.shape[0], Tsys.shape[1]
FG = tb.getcol('FLAG')
SpwID = SpwID.reshape(spwNum, scanNum, antNum)
AntID = AntID.reshape(spwNum, scanNum, antNum)
ScanID = ScanID.reshape(spwNum, scanNum, antNum)
#-------- Identify copy source and target addresses
polIDList = [polDic[polID] for polID in polList]
SpwIDList = np.where(SpwID[:,0,0] == spws)[0].tolist()
antIdIndex = np.where(AntID[0,0] == ant_index)[0][0]
ScanSourceID = np.where(ScanID[0][:,0] == scans[0])[0][0]
ScanDestID   = np.where(ScanID[0][:,0] == scans[1])[0][0]
SrcIndex = [index for index in list(range(spwNum*scanNum*antNum)) if (index%antNum == antIdIndex) & ((index//antNum)%scanNum == ScanSourceID) and ((index//(scanNum*antNum)) in SpwIDList)]
DstIndex = [index for index in list(range(spwNum*scanNum*antNum)) if (index%antNum == antIdIndex) & ((index//antNum)%scanNum == ScanDestID) and ((index//(scanNum*antNum)) in SpwIDList)]
#-------- Copy
for pol in polIDList:
    for index, Dst in enumerate(DstIndex):
        Tsys[pol][:,Dst] = Tsys[pol][:,SrcIndex[index]]
        FG[pol][:,Dst] = FG[pol][:,SrcIndex[index]]
#-------- Save
tb.putcol('FPARAM', Tsys)
tb.putcol('FLAG', FG)
tb.close()
print('Antenna=%s' % (antName))
print('POL=',end=''); print(polList)
print('SPW=',end=''); print(spws)
print('Tsys source scan=%d' % (scans[0]))
print('Tsys destination scan=%d' % (scans[1]))
