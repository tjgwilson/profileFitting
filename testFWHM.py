import numpy as np
import matplotlib.pyplot as plt
import gaussianFitting as gf

fitter = gf.fitGaussian("halpha_020_06562.dat",verbose=True)

fitter.calcFWHM()
print(fitter.FWHM)
