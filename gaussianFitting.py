import numpy as np
from astropy import modeling
from astropy.modeling.models import Gaussian1D
class fitGaussian():
    def __init__(self,filename=None,data=None,col=1,addNoise=False,sensitivity=1,verbose=False):
        self.sensitivity = sensitivity / 100.
        self.verbose = verbose
        if(filename != None):
            input = np.loadtxt(filename)
            self.x = input[:,0]
            self.y = input[:,col]
        elif(data != None):
            self.x = data[:,0]
            self.y = data[:,col]
        else:
            print("ERROR: No imput data")
        if(addNoise):
            noise = np.random.normal(1,0.05*np.amax(self.y),self.x.size)
            self.y += noise
        self.cont = (self.y[0] + self.y[-1])/2.0
        self.y -= self.cont

        self.fit = fitGaussian.createFit(self)

        if(self.verbose):
        ####debugging#####
            import matplotlib.pyplot as plt
            print(self.fit)
            plt.plot(self.x,self.y)
            for ii in range(self.nMin):
                plt.plot(self.xMin[ii],self.yMin[ii],'ro',markersize=5)
            plt.plot(self.xMax,self.yMax,'ko',markersize=5)
            plt.plot(self.x,self.fit(self.x),'g',alpha=0.5)
            plt.show()

    def createFit(self):
        fitGaussian.findMinima(self)
        fitGaussian.findMaxima(self)

        fitter = modeling.fitting.LevMarLSQFitter()
        g1 = Gaussian1D(self.yMax,self.xMax,10)
        g2 = Gaussian1D(self.yMax,self.xMax,5)
        model = g1 + g2
        gaus = []
        for ii in range(self.nMin):
            gaus.append(Gaussian1D(self.yMin[ii],self.xMin[ii],5))
        for g in gaus:
            model += g
        fitted_model = fitter(model, self.x, self.y)
        return fitted_model
        #find 0,1 or 2 sub continuum minima on either side of the 0 point.
    def findMinima(self):
        first = True
        prec = - self.sensitivity
        midii = int(len(self.x)/2.0)
        self.nMin = 0
        self.yMin = []
        self.xMin = []
        tMin = 0.0
        tXMin = 0.0
        for ii in range(len(self.x)):
            if((self.y[ii] <= tMin)):
                tMin = self.y[ii]
                tXMin = self.x[ii]
            if((ii > midii) and first):
                first = not first
                if((tMin <= prec) and (ii != 0) and (tXMin != self.x[midii])):
                    self.yMin.append(tMin)
                    self.xMin.append(tXMin)
                    self.nMin += 1
                    tMin = 0.0
                    tXMin = 0.0
        if((tMin <= prec) and (tXMin < self.x[-1]) and (tXMin > self.x[midii])):
            self.yMin.append(tMin)
            self.xMin.append(tXMin)
            self.nMin += 1
        #find the peak emission flux and velocity of the peak
    def findMaxima(self):
        tMax = 0.0
        tXMax = 0.0
        for ii in range(len(self.x)):
            if(self.y[ii] > tMax):
                tMax = self.y[ii]
                tXMax = self.x[ii]
        self.xMax = tXMax
        self.yMax = tMax
        #Finds the FWHM for the main fit and all the component gaussians: n is the number of points used to create the loop
        #precision is the sensitivity to which it detect the required points on the curve
    def calcFWHM(self,n=10000,precision=0.01):
        import matplotlib.pyplot as plt #remove when debugged
        self.FWHM = []
        x = np.linspace(np.amin(self.x),np.amax(self.x),n,endpoint=True)
        y = self.fit(x)
        halfMax = 0.5*np.amax(y)
        x1 = 0.0
        x2 =0.0
        for ii in range(len(x)): #calculates the FWHM for the main gaussian fit
            if(abs(y[ii] - halfMax) <= precision):
                if(x[ii] < 0.0):
                    x1 = x[ii]
                if(x[ii] >= 0.0):
                    x2 = x[ii]
        self.FWHM.append(x2-x1)

        if(self.verbose):
            plt.plot((x1,x2),(halfMax,halfMax),'r-o',markersize=5)
            plt.plot(x,y)

        for f in self.fit: #calculates the FWHM for the component gaussians
            y = f(x)
            halfMax = 0.5*f.amplitude[0]
            x1 = 0.0
            x2 = 0.0
            first = True
            for ii in range(len(x)):
                if(abs(y[ii] - halfMax) <= precision):
                    if(first):
                        x1 = x[ii]
                        first = False
                    else:
                        x2 = x[ii]
            self.FWHM.append(x2-x1)
            if(self.verbose):
                plt.plot(x,y)
                plt.plot((x1,x2),(halfMax,halfMax),'r-o',markersize=5)


        if(self.verbose):
            plt.plot(self.x,self.y,'k--',alpha=0.5,label="Torus data")
            plt.xlabel("velocity")
            plt.ylabel("flux")
            plt.legend()
            plt.show()



    def lineParameters(self):
        pass
        HWZM = []
        peakFLux = []
        centre = []
        eqWidth = []
