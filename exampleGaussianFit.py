import numpy as np
import matplotlib.pyplot as plt
import line_fitting as gf

filename = "halpha_020_06562.dat"

print("####################Gaussian####################")
ft = gf.fitGaussian(filename=filename,verbose=False)

ft.calcFWHM()
ft.calcHWZM()
ft.calcCentre()
ft.calcEqWidths(wavelength=656.2)
print("FWHM",ft.FWHM[0])
print("HWZM",ft.HWZM[0])
print("centre",ft.centres[0])
print("EqWidth",ft.eqWidths)
