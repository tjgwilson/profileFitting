import numpy as np
import matplotlib.pyplot as plt


class folhaParam():

    def __init__(self,inputFile):
        self.inputData = np.loadtxt(inputFile)



folhaParam("halpha_020_06562.dat")
