from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = options.spwList
os.system('rm -rf WVR')
wvrgcal(vis=prefix + '.ms', caltable='WVR', wvrspw=[4], toffset=0, statsource='', wvrflag=[])
applycal(vis=prefix + '.ms', spw=spwList, interp='nearest', gaintable='WVR', spwmap=[], calwt=True, flagbackup=False)
os.system('rm -rf ' + prefix + '.WVR.ms')
split(vis=prefix + '.ms', spw=spwList, datacolumn='corrected', outputvis=prefix+'.WVR.ms')
