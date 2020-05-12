import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class folhaParam():

    def __init__(self,inputFile):
        self.inputData = np.loadtxt(inputFile)


    def fitParabola(self):

        folhaParam.profileBounds(self)
        x = self.inputData[:,0]
        y = self.inputData[:,1]
        xData = self.inputData[self.lBi:self.uBi,0]
        yData = self.inputData[self.lBi:self.uBi,1]

        popt, pcov = curve_fit(folhaParam.parabola,xData, yData)
        plt.plot(xData,folhaParam.parabola(xData,*popt))
        plt.plot(x,y)
        plt.plot(self.lB,y[self.lBi],'ro')
        plt.plot(self.uB,y[self.uBi],'ro')
        plt.show()

    def profileBounds(self):
        x = self.inputData[:,0]
        y = self.inputData[:,1]
        prec = np.amax(y)/2
        l = len(y)
        for ii in range(0,l,1):
            if((y[ii+1] -1.0) >= prec):
                # if((y[ii+1] - y[ii]) > 0.0):
                self.lB = x[ii]
                self.lBi = ii
                break

        for ii in range(l-1,0,-1):
            if((y[ii-1] - 1.0) >= prec):
                # if((y[ii-1] - y[ii]) > 0.0):
                self.uB = x[ii]
                self.uBi = ii
                break

        # for ii in range(0,l,1):
        #     if(abs(y[ii] - y[ii+1]) >= prec):
        #         self.lB = x[ii]
        #         self.lBi = ii
        #
        #         break
        # for ii in range(l-1,0,-1):
        #     if(abs(y[ii] - y[ii-1]) >= prec):
        #         self.uB = x[ii]
        #         self.uBi = ii
        #         break


    def parabola(x,a0,a1,a2):
        return (a0 + a1*x + a2*x*x)


folhaParam("halpha_020_06562.dat").fitParabola()
