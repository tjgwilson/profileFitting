import numpy as np
from astropy import modeling
import matplotlib.pyplot as plt

#####Importing#####
col = 1
inc = ['20','50','80']
freq = ['06562','12818','21655']
filename = []
for i in inc:
    for f in freq:
        filename.append("halpha_0"+i+"_"+f+".dat")

fig, ax = plt.subplots(3,3)
r=0
c=0
check = 1
for name in filename:
    data = np.loadtxt(name)
    # data[:,1] = data[:,1] + np.sqrt(abs(data[:,1])) * np.random.random(data[:,0].size) - 0.5
    data[:,1] = data[:,1] - data[0,1]


    ####Fitting######

    max = 0.0
    xmax = 0.0
    min1 = 0.0
    xmin1 = 0.0
    min2 = 0.0
    xmin2 = 0.0

    for i in range(len(data)):
        if(data[i,col] >= max):
            max = data[i,col]
            xmax = data[i,0]
        if((data[i,col] <= min1) and (data[i,0] < 0.0)):
            min1 = data[i,col]
            xmin1 = data[i,0]
        if((data[i,col] <= min2) and (data[i,0] > 0.0)):
            min2 = data[i,col]
            xmin2 = data[i,0]
    print(min1,xmin1)


    fitter = modeling.fitting.LevMarLSQFitter()
    gaus =[]
    gaus.append(modeling.models.Gaussian1D(max,xmax,10)) #- modeling.models.Gaussian1D()  # Lorentz1D() Voigt1D() Moffatt1D()
    gaus.append(modeling.models.Gaussian1D(max,xmax,5))

    if(min1 != 0):
        gaus.append(modeling.models.Gaussian1D(min1,xmin1,10))
    if(min2 != 0):
        gaus.append(modeling.models.Gaussian1D(min2,xmin2,10))
    model = gaus[0] + gaus[1]
    for i in range(2,len(gaus)):
        model += gaus[i]
    fitted_model = fitter(model, data[:,0], data[:,col])
    # for model in fitted_model:
    #     print(model)


    #####Ploting######

    ax[r,c].plot(data[:,0], data[:,1],'r-')
    ax[r,c].plot(data[:,0], fitted_model(data[:,0]),'b--')
    c += 1
    if(check%3 == 0):
        r += 1
        c = 0
    check += 1
plt.show()
