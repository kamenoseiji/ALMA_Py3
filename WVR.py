import os
from interferometry import GetBPchavSPWs
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
msfile  = prefix + '.ms'
spwList = [] if options.spwList == '' else [int(spw) for spw in options.spwList.split(',')]
if len(spwList) < 1: spwList = GetBPchavSPWs(msfile)
spwChar = ''
for spw in spwList: spwChar = spwChar + '%d,' % (spw)
spwChar = spwChar[:-1]
print(spwChar)
msmd.open(msfile)
WVRspw = msmd.wvrspws().tolist()
msmd.done()
if len(WVRspw) > 0:
    os.system('rm -rf %s.WVR' % (prefix))
    wvrgcal(vis=prefix + '.ms', caltable=prefix + '.WVR', wvrspw=WVRspw, toffset=0, statsource='', wvrflag=[], minnumants=2 )
    applycal(vis=prefix + '.ms', spw=spwChar, interp='nearest', gaintable=prefix + '.WVR', spwmap=[], calwt=True, flagbackup=False)
    os.system('rm -rf ' + prefix + '.WVR.ms')
    split(vis=prefix + '.ms', spw=spwChar, datacolumn='corrected', outputvis=prefix+'.WVR.ms')
#
