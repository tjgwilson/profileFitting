import numpy as np
from astropy import modeling
from astropy.modeling.models import Gaussian1D, Voigt1D
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import trapz
#copy to /usr/local/lib/python3.9/site-packages to allow other code to import

class fitGaussian():
    def __init__(self,filename=None,data=None,col=1,addNoise=False,sensitivity=1,verbose=False,continuum=1.0,blueOnly=False):
        self.c = 299792.458
        self.sensitivity = sensitivity / 100.
        self.verbose = verbose
        self.cont = continuum
        self.blueOnly = blueOnly
        if(filename != None):
            input = np.loadtxt(filename)
            self.x = input[:,0]
            self.y = input[:,col]
            self.filename = filename
        elif(data.any() != None):
            self.x = data[:,0]
            self.y = data[:,col]
        else:
            print("ERROR: No imput data")
        if(addNoise):
            noise = np.random.normal(1,0.05*np.amax(self.y),self.x.size)
            self.y += noise
        self.y = np.subtract(self.y,self.cont)

        self.fit = fitGaussian.createFit(self)
        fitGaussian.calcFWHM(self)
        fitGaussian.calcPeakFlux(self)
        if(self.blueOnly):
            self.fit = fitGaussian.createBlueFit(self)

        if(self.verbose):
            print(self.fit)
            plt.plot(self.x,self.y,'k--',alpha=1,label="Input Data")
            for ii in range(self.nMin):
                plt.plot(self.xMin[ii],self.yMin[ii],'ro',markersize=5,label="Sub continuum found")
            plt.plot(self.xMax[0],self.yMax[0],'ko',markersize=5,label="Peak Max found")
            plt.plot(self.xMax[1],abs(self.yMax[1])/2,'co',markersize=5,label="Local Minima found")
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
        g1.amplitude.fixed = True
        g2 = Gaussian1D(self.yMax[1],self.xMax[1],5)

        model = g1 + g2

        for ii in range(self.nMin):
            model += Gaussian1D(self.yMin[ii],self.xMin[ii],5)

        fitted_model = fitter(model, self.x, self.y)
        return fitted_model


    def createBlueFit(self):
        from scipy.optimize import curve_fit

        def Gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        fitGaussian.findMaxima(self)
        fitter = modeling.fitting.LevMarLSQFitter()

        idx = simpleFit.closestIndex(self.y,np.amax(self.y))+1
        yy1 = self.y[:idx]
        yy2 = yy1[::-1]
        yy = np.append(yy1,yy2)

        origin = self.x[idx]
        xx1 = self.x[:idx] - origin
        xx2 = -1.0*(xx1[::-1])
        xx = np.append(xx1,xx2)
        for ii in range(len(yy)):
            if(yy[ii] < 0.0):
                yy[ii] = 0.0
        self.x = xx
        self.y = yy

        mean = sum(xx * yy) / sum(yy)
        sigma = np.sqrt(sum(yy * (xx - mean)**2) / sum(yy))

        g1 = Gaussian1D(amplitude=np.amax(yy),mean=mean,stddev=sigma)
        g1.amplitude.fixed = True

        fitted_model = fitter(g1, xx, yy)
        self.x = xx
        self.y = yy
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
            self.yMax.append(tMin[count]*2)
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
        n = len(self.x)*100
        precision = np.amax(self.y)/100
        x = np.linspace(np.amin(self.x),np.amax(self.x),n,endpoint=True)
        y = self.fit(x)
        halfMax = 0.5*np.amax(y)

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
            plt.plot((x1,x2),(halfMax,halfMax),'r-o',markersize=5)
            plt.plot(x,y)

        if(not self.blueOnly):
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
        return self.FWHM
        #Finds the HWZM for the main fit and all the component gaussians: n is the number of points used to create the loop##[km/s]
        #precision is the sensitivity to which it detect the required points on the curve, detects at 2% peak hight.##[km/s]
    def calcHWZM(self,n=100000,precision=0.001,zeroMag = 0.1):
        self.HWZM = []
        n = len(self.x)*100
        precision = np.amax(self.y)/100
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
        if(not self.blueOnly):
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
        return self.HWZM
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
            for c in self.centres[0:1]:
                plt.axvline(c)
            plt.show()
        #Calculates the peak flux for the fit and the gaussian components, stored into self.peaks ##
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
            for c in self.peaks[0:1]:
                plt.axhline(y=c)
            plt.show()
        #calculates the equivilant width of the fitted line and the component gaussians
        #output is in angstroms - converts velocity space to wavelength space
        #prio to integration via trapeziumm method. ###[\AA]
    def calcEqWidths(self,wavelength,bins=100000,dx=0.1):
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
        return self.eqWidths

    def closestIndex(array,value):
        array = np.asarray(array)
        idx = (np.abs(array-value)).argmin()
        return idx

class simpleFit():

    def __init__(self,filename=None,data=None,col=1,sensitivity=1,verbose=False,continuum=1.0,blueOnly=False):
        self.c = 299792.458
        self.sensitivity = sensitivity / 100.
        self.verbose = verbose
        self.cont = continuum

        if(filename != None):
            input = np.loadtxt(filename)
            self.x = input[:,0]
            self.y = input[:,col]
            self.filename = filename
        elif(data.any() != None):
            self.x = data[:,0]
            self.y = data[:,col]
        else:
            print("ERROR: No imput data")
        self.y = np.subtract(self.y,self.cont)
        if(blueOnly):
            self.fit = simpleFit.splineFitBlue(self)
        else:
            self.fit = simpleFit.splineFit(self)
        self.peak = simpleFit.calcPeakFlux(self)
        self.min = np.amin(self.y)
        self.centre = simpleFit.calcCentre(self)
        self.fwhm = simpleFit.calcFWHM(self)
        self.hwzm = simpleFit.calcHWZM(self)

        if((self.peak <= 0.0) or (self.fwhm <= 0.0) or (self.hwzm <= 0.0) or (abs(self.peak/self.min) < 0.01)):
            print("No emission detected")
            self.peak = 0.0
            self.centre = 0.0
            self.fwhm = 0.0
            self.hwzm = 0.0
        if(self.verbose):
            xx = np.linspace(np.amin(self.x),np.amax(self.x),100000,endpoint=True)
            plt.plot(self.x,self.y,'k-')
            plt.plot(xx,self.fit(xx),'r--')
            plt.plot(self.centre,self.peak,'g*',markersize=10)
            plt.show()

    def splineFit(self):
        fit = interp1d(self.x,self.y,kind='linear')
        return fit

    def splineFitBlue(self):

        idx = simpleFit.closestIndex(self.y,np.amax(self.y))+1
        yy1 = self.y[:idx]
        yy2 = yy1[::-1]
        yy = np.append(yy1,yy2)

        origin = self.x[idx]
        xx1 = self.x[:idx] - origin
        xx2 = -1.0*(xx1[::-1])
        xx = np.append(xx1,xx2)
        for ii in range(len(yy)):
            if(yy[ii] < 0.0):
                yy[ii] = 0.0
        self.x = xx
        self.y = yy


        fit = interp1d(xx,yy,kind='linear')
        return fit

    def calcPeakFlux(self):
        return np.amax(self.y)

    def calcCentre(self):
        ii = np.argmax(self.y)
        return self.x[ii]

    def closestIndex(array,value):
        array = np.asarray(array)
        idx = (np.abs(array-value)).argmin()
        return idx


    def calcFWHM(self,precision=0.000001):
        if(self.peak > 0.0):
            hp = 0.5*np.amax(self.y)
            matches = []
            found = True

            xnum = len(self.x)*1000
            precision = np.amax(self.y)/1000.
            xx = np.linspace(np.amin(self.x),np.amax(self.x),xnum,endpoint=True)
            yy = self.fit(xx)

            for ii in range(0,len(yy)-1,1):
                if((abs(yy[ii]-hp) <= precision)):
                    # matches.append(ii)
                    dif = abs(yy[ii]-hp)
                    if(abs(yy[ii+1]-hp) < dif):
                        continue
                    else:
                        matches.append(ii)
            x1 = xx[matches[0]]
            x2 = xx[matches[-1]]

            if((matches[0] < int(xnum*0.1)) or (matches[-1] > int(xnum*0.9))):
                x1 = 0.0
                x2 = 0.0
        else:
            x1 = 0.0
            x2 = 0.0

        if(self.verbose):
            print("FWHM Matches",matches)
            plt.plot(self.x,self.y,'k-')
            plt.plot(x1,self.fit(x1),'rp',markersize=5)
            plt.plot(x2,self.fit(x2),'rp',markersize=5)
            plt.show()
        return x2-x1

    def calcHWZM(self,precision=0.000001,zeroMag = 0.1):
        if(self.peak > 0.0):
            hp = zeroMag*np.amax(self.y)
            matches = []
            found = True

            xnum = len(self.x)*1000
            precision = np.amax(self.y)/1000.
            xx = np.linspace(np.amin(self.x),np.amax(self.x),xnum,endpoint=True)
            yy = self.fit(xx)

            for ii in range(0,len(yy)-1,1):
                if((abs(yy[ii]-hp) <= precision)):
                    # matches.append(ii)
                    dif = abs(yy[ii]-hp)
                    if(abs(yy[ii+1]-hp) < dif):
                        continue
                    else:
                        matches.append(ii)
            x1 = xx[matches[0]]
            x2 = xx[matches[-1]]

            if((matches[0] < int(xnum*0.1)) or (matches[-1] > int(xnum*0.9))):
                x1 = 0.0
                x2 = 0.0
        else:
            x1 = 0.0
            x2 = 0.0

        if(self.verbose):
            print("HWZM Matches",matches)
            plt.plot(self.x,self.y,'k-')
            plt.plot(x1,self.fit(x1),'rp',markersize=5)
            plt.plot(x2,self.fit(x2),'rp',markersize=5)
            plt.show()
        return (x2-x1)*0.5

    def calcEqWidths(self,wavelength,dx=0.1):
        if(self.peak > 0.0):
            eqWidths = 0.0
            wl = (self.x * wavelength / self.c) + wavelength
            dx = (abs(np.amax(wl)) + abs(np.amin(wl))) / float(len(wl))
            y = (1.0 - ((self.y+self.cont)/self.cont))
            eqWidths = trapz(y,wl,dx=dx)
        else:
            eqWidths = 0.0
        return eqWidths
