class fitGaussian():
    def __init__(self,filename=None,data=None,col=1,addNoise=False,sensitivity=1,verbose=False):
        import numpy as np
        self.sensitivity = sensitivity / 100.
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

        if(verbose):
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
        from astropy import modeling
        fitGaussian.findMinima(self)
        fitGaussian.findMaxima(self)

        fitter = modeling.fitting.LevMarLSQFitter()
        g1 = modeling.models.Gaussian1D(self.yMax,self.xMax,10)
        g2 = modeling.models.Gaussian1D(self.yMax,self.xMax,5)
        model = g1 + g2
        gaus = []
        for ii in range(self.nMin):
            gaus.append(modeling.models.Gaussian1D(self.yMin[ii],self.xMin[ii],5))
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





#####TEST USE OF FITTER###########
import matplotlib.pyplot as plt
inc = ['20','60','80']
freq = ['06562','12818','21655']
filename = []
for i in inc:
    for f in freq:
        filename.append("halpha_0"+i+"_"+f+".dat")

fig, ax = plt.subplots(len(inc),len(freq))
r=0
c=0
check = 1
for name in filename:
    data = fitGaussian(name)
    x = data.x
    y = data.y
    fit = data.fit
    ax[r,c].plot(x,y,'b--')
    ax[r,c].plot(x,fit(x),'r-')

    c += 1
    if(check%len(inc) == 0):
        r += 1
        c = 0
    check += 1
plt.show()
