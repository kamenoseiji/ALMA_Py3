import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
exec(open(SCR_DIR + 'Plotters.py').read())
#
FreqList, DxList, DyList = [], [], []
for spw_index, spw in enumerate(spwList):
    antList = np.load('%s-SPW%d-%s.Ant.npy' % (prefix, spw, refantName))
    for ant_index, antName in enumerate(antList):
        Dterm = np.load('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, antName))
        DxList   = DxList   + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList   = DyList   + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
    FreqList = FreqList + [Dterm[0]]
#
pp = PdfPages('D_%s-REF%s-Dspec.pdf' % (prefix, refantName))
plotDSpec(pp, prefix, antList, spwList, FreqList, DxList, DyList)
