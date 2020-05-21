import numpy as np
from astropy import modeling
from astropy.modeling.models import Gaussian1D
import matplotlib.pyplot as plt

class fitGaussian():
    def __init__(self,filename=None,data=None,col=1,addNoise=False,sensitivity=1,verbose=False):
        self.c = 299792.458
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
        self.filename = filename

        if(self.verbose):
            print(self.fit)
            plt.plot(self.x,self.y,'k--',alpha=1,label="Input Data")
            for ii in range(self.nMin):
                plt.plot(self.xMin[ii],self.yMin[ii],'ro',markersize=5,label="Sub continuum found")
            plt.plot(self.xMax[0],self.yMax[0],'ko',markersize=5,label="Peak Max found")
            plt.plot(self.xMax[1],abs(self.yMax[1]),'co',markersize=5,label="Local Minima found")
            plt.plot(self.x,self.fit(self.x),'b',alpha=0.7,label="Gaussian Fit")
            for f in self.fit:
                plt.plot(self.x,f(self.x),alpha=0.7,label="component Gaussian")
            plt.title(filename)
            plt.legend(fontsize=10)
            plt.show()

    def createFit(self):
        fitGaussian.findMinima(self)
        fitGaussian.findMaxima(self)
        fitter = modeling.fitting.LevMarLSQFitter()

        g1 = Gaussian1D(self.yMax[0],self.xMax[0],10)
        g2 = Gaussian1D(self.yMax[1],self.xMax[1],5)
        if(self.localMinima):
            model = g1 - g2
        else:
            model = g1 + g2

        for ii in range(self.nMin):
            model += Gaussian1D(self.yMin[ii],self.xMin[ii],5)

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
        #find the peak emission flux and velocity of the peak and local mimima - it only reports one local mimia.
        # the smallest there is
    def findMaxima(self):
        prec = self.sensitivity
        self.yMax = []
        self.xMax = []
        tMax = 0.0
        tXMax = 0.0
        peakii = 0
        for ii in range(len(self.x)):
            if(self.y[ii] > tMax):
                tMax = self.y[ii]
                tXMax = self.x[ii]
                peakii = ii
        self.yMax.append(tMax)
        self.xMax.append(tXMax)

        localMinima = False
        tMin = []
        tXMin = []
        for ii in range(1,len(self.x)):
            if(self.y[ii] > prec):
                grad = (self.y[ii] - self.y[ii-1]) / (self.x[ii] - self.x[ii-1])
                if((grad < 0.0) and (ii < peakii)):
                    # plt.plot(self.x[ii],self.y[ii],'rp',alpha=0.5)
                    tMin.append(self.y[ii])
                    tXMin.append(self.x[ii])
                    localMinima = True
                if((grad > 0.0) and (ii > peakii)):
                    # plt.plot(self.x[ii],self.y[ii],'bp',alpha=0.5)
                    tMin.append(self.y[ii])
                    tXMin.append(self.x[ii])
                    localMinima = True
        if(localMinima):
            temp = tMin[0]
            depth = 0.0
            count = 0
            for ii in range(1,len(tMin)):
                if(tMin[ii] < temp):
                    temp = tMin[ii]
                    count += 1
            self.yMax.append(-tMin[count])
            self.xMax.append(tXMin[count])
        else:
            self.yMax.append(tMax)
            self.xMax.append(tXMax)

        # plt.plot(self.x,self.y,'k.')
        # plt.show()
        self.localMinima = localMinima



        #Finds the FWHM for the main fit and all the component gaussians: n is the number of points used to create the loop
        #precision is the sensitivity to which it detect the required points on the curve
    def calcFWHM(self,n=100000,precision=0.001):
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
        #Finds the HWZM for the main fit and all the component gaussians: n is the number of points used to create the loop##[km/s]
        #precision is the sensitivity to which it detect the required points on the curve, detects at 2% peak hight.##[km/s]
    def calcHWZM(self,n=100000,precision=0.001,zeroMag = 0.02):

        self.HWZM = []
        x = np.linspace(np.amin(self.x),np.amax(self.x),n,endpoint=True)
        y = self.fit(x)
        zeroMax = zeroMag*np.amax(y)
        x1 = 0.0
        x2 =0.0
        for ii in range(len(x)): #calculates the FWHM for the main gaussian fit
            if(abs(y[ii] - zeroMax) <= precision):
                if(x[ii] < 0.0):
                    x1 = x[ii]
                if(x[ii] >= 0.0):
                    x2 = x[ii]
        self.HWZM.append((x2-x1)/2.0)

        if(self.verbose):
            plt.plot((x1,x2),(zeroMax,zeroMax),'r-o',markersize=5)
            plt.plot(x,y)

        for f in self.fit: #calculates the FWHM for the component gaussians
            y = f(x)
            zeroMax = zeroMag*f.amplitude[0]
            x1 = 0.0
            x2 = 0.0
            first = True
            for ii in range(len(x)):
                if(abs(y[ii] - zeroMax) <= precision):
                    if(first):
                        x1 = x[ii]
                        first = False
                    else:
                        x2 = x[ii]
            self.HWZM.append((x2-x1)/2.0)
            if(self.verbose):
                plt.plot(x,y)
                plt.plot((x1,x2),(zeroMax,zeroMax),'r-o',markersize=5)


        if(self.verbose):
            plt.plot(self.x,self.y,'k--',alpha=0.5,label="Torus data")
            plt.xlabel("velocity")
            plt.ylabel("flux")
            plt.legend()
            plt.show()
        #Calculate centres of fit, and gaussian components, stored in self.centres ##[km/s]
    def calcCentre(self,n=100000,precision=0.001):
        self.centres = []
        x = np.linspace(np.amin(self.x),np.amax(self.x),n,endpoint=True)
        y = self.fit(x)
        ty = 0.0
        tx = 0.0
        for ii in range(len(x)):
            if(y[ii] >= ty):
                ty = y[ii]
                tx = x[ii]
        self.centres.append(tx)
        for f in self.fit:
            self.centres.append(f.mean[0])

        if(self.verbose):
            plt.plot(self.x,self.y)
            plt.plot(self.x,self.fit(self.x))
            for f in self.fit:
                plt.plot(x,f(x))
            for c in self.centres:
                plt.axvline(c)
            plt.show()
        #Calculates the peak flux for the fit and the gaussian components, stored into self.peaks ##[km/s]
    def calcPeakFlux(self,n=100000,precision=0.001):
        self.peaks = []
        x = np.linspace(np.amin(self.x),np.amax(self.x),n,endpoint=True)
        y = self.fit(x)
        ty = 0.0
        for ii in range(len(x)):
            if(y[ii] >= ty):
                ty = y[ii]
        self.peaks.append(ty)
        for f in self.fit:
            self.peaks.append(f.amplitude[0])

        if(self.verbose):
            plt.plot(self.x,self.y)
            plt.plot(self.x,self.fit(self.x))
            for f in self.fit:
                plt.plot(x,f(x))
            for c in self.peaks:
                plt.axhline(y=c)
            plt.show()
        #calculates the equivilant width of the fitted line and the component gaussians
        #output is in angstroms - converts velocity space to wavelength space
        #prio to integration via trapeziumm method. ###[\AA]
    def calcEqWidths(self,wavelength,bins=100000,dx=0.1):
        from scipy.integrate import trapz
        self.eqWidths = []
        wl= (self.x * wavelength / self.c) + wavelength
        x = np.linspace(np.amin(self.x),np.amax(self.x),bins,endpoint=True)
        w = np.linspace(np.amin(wl),np.amax(wl),bins,endpoint=True)
        y = (1.0 - ((self.fit(x)+self.cont)/self.cont))
        self.eqWidths.append(trapz(y,w,dx=dx))
        for f in self.fit:
            y = (1.0 - ((f(x)+self.cont)/self.cont))
            self.eqWidths.append(trapz(y,w,dx=dx))


        if(self.verbose):
            fix, ax = plt.subplots(2,sharex=True)
            y = self.fit(x)+self.cont
            ax[0].plot(w,(self.fit(x)+self.cont),'r-')
            ax[0].fill_between(w,self.cont,(self.fit(x)+self.cont),color='r',alpha=0.3)
            x1 = wavelength-(self.eqWidths[0]/2.)
            x2 = wavelength+(self.eqWidths[0]/2.)
            ax[0].plot((x1,x1),(0.,self.cont),color='g')
            ax[0].plot((x2,x2),(0.,self.cont),color='g')
            ax[0].fill_between((x1,x2),self.cont,0,color='g',alpha=0.3)

            c = ['c','g','b','r']
            ii = 0
            print(np.sum(np.abs(self.eqWidths)))
            x1 = wavelength-((np.sum(np.abs(self.eqWidths))-abs(self.eqWidths[0]))/2.)
            for f in self.fit:
                y = f(x)+self.cont
                ax[1].plot(w,(f(x)+self.cont),color=c[ii])
                ax[1].fill_between(w,self.cont,(f(x)+self.cont),color=c[ii],alpha=0.3)
                x2 = x1 + abs(self.eqWidths[ii+1])
                ax[1].plot((x1,x1),(0.,self.cont),color=c[ii])
                ax[1].plot((x2,x2),(0.,self.cont),color=c[ii])
                ax[1].fill_between((x1,x2),self.cont,0,color=c[ii],alpha=0.3)
                x1 = x2
                ii += 1
            ax[1].set_xlabel(r"Wavelength [$\AA$]")
            ax[0].set_title("Equivilant widths")
            plt.show()
