import numpy as np
import matplotlib.pyplot as plt
import gaussianFitting as gf


inc = ['20','60','80']
freq = ['06562','12818','21655']
# inc = ["20"]
# freq= ["21655"]
for i in inc:
    for f in freq:
        filename= "halpha_0"+i+"_"+f+".dat"

        fitter = gf.fitGaussian(filename,verbose=False)
        # fitter.calcFWHM()
        # fitter.calcHWZM()
        # fitter.calcCentre()
        fitter.calcPeakFlux()
        if(fitter.nMin > 0):
            print(fitter.peaks[3],f)
        # fitter.calcEqWidths(int(f))
        # print(fitter.eqWidths)
