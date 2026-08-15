import sys
import os
import datetime
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-s', dest='SBcode', metavar='SBcode',
    help='SB code  e.g. PSO_Grid', default='PSO_Grid')
parser.add_option('-b', dest='Subject', metavar='Subject',
    help='Subject  e.g. baseline', default='')
parser.add_option('-d', dest='days', metavar='days',
    help='Backword days', default='4')
(options, args) = parser.parse_args()
backDays=  int(options.days)
UserPass = os.getenv('b64credentials', '')   # Authorization Basic Base64
if UserPass == '':
    print('Set b64credentials as an environmental variable')
    exit()
SLT_URI  = 'https://asa.alma.cl/webslt/service/api/entries?'
SLTstart = (datetime.datetime.today() - datetime.timedelta(days=backDays)).strftime('%Y-%m-%dT%H:%M:%S')
SLTend   = (datetime.datetime.today() - datetime.timedelta(hours=2)).strftime('%Y-%m-%dT%H:%M:%S')
#queryText = 'curl -H\"Authorization: Basic %s\" \'%sintervalStart=%s&intervalEnd=%s&entryType=SBEX&schedBlockCode=%s&status=success\' > SLT.log' % (UserPass, SLT_URI, SLTstart, SLTend, options.SBcode)
if options.Subject == '':
    queryText = 'curl -H\"Authorization: Basic %s\" \'%sintervalStart=%s&intervalEnd=%s&entryType=SBEX&schedBlockCode=%s&status=success\' > SLT.log' % (UserPass, SLT_URI, SLTstart, SLTend, options.SBcode)
else:
    queryText = 'curl -H\"Authorization: Basic %s\" \'%sintervalStart=%s&intervalEnd=%s&subject=%s&status=success\' > SLT.log' % (UserPass, SLT_URI, SLTstart, SLTend, options.Subject)
#print(queryText)
os.system(queryText)
fp = open('SLT.log', 'r')
SLTline = fp.readlines()
fp.close()
#if len(SLTline) > 0:
SLTentries = SLTline[0].split('\"')
UIDList = [SLTentries[index+2] for index, SLTentry in enumerate(SLTentries) if 'execBlockUid' in SLTentry]
arrayNameList = [SLTentries[index+2] for index, SLTentry in enumerate(SLTentries) if 'arrayName' in SLTentry]
startTimeList = [SLTentries[index+2] for index, SLTentry in enumerate(SLTentries) if 'start' in SLTentry]
endTimeList   = [SLTentries[index+2] for index, SLTentry in enumerate(SLTentries) if 'end' in SLTentry]
for index, UID in enumerate(UIDList): print('new %s %s %s %s' % (UID, arrayNameList[index], startTimeList[index], endTimeList[index]))
