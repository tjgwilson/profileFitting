import numpy as np
import matplotlib.pyplot as plt
import line_fitting as gf

filename = "/Users/twilson/Documents/Torus/Stellar_Wind/finalgrid/Tacc7500_Tsw10000_Macc9_Msw10/halpha_turb_020_10938.dat"
filename = "/Users/twilson/Documents/Torus/Stellar_Wind/finalgrid/Tacc6500_Tsw6000_Macc9_Msw10/halpha_060_21655.dat"




print("####################Gaussian####################")
ft = gf.fitGaussian(filename=filename,verbose=True,blueOnly=True)

ft.calcFWHM()
ft.calcHWZM()
print("FWHM",ft.FWHM[0])
print("HWZM",ft.HWZM[0])
print("EqWidth",ft.eqWidths)
