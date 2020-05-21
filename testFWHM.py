import numpy as np
import matplotlib.pyplot as plt
import gaussianFitting as gf


inc = ['20','60','80']
freq = ['06562','12818','21655']
# inc = ["60"]
# freq= ["12818"]
for i in inc:
    for f in freq:
        filename= "halpha_0"+i+"_"+f+".dat"

        fitter = gf.fitGaussian(filename,verbose=True)
        # fitter.calcFWHM()
        # fitter.calcHWZM()
        # fitter.calcCentre()
        # fitter.calcPeakFlux()
        # fitter.calcEqWidths(int(f))
        # print(fitter.eqWidths)
