from line_fitting import fitGaussian, simpleFit

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
    # data = fitGaussian(name,verbose=False,blueOnly=True)
    data = simpleFit(name,verbose=False,blueOnly=True)
    # print(data.calcFWHM()[0])
    x = data.x
    y = data.y
    fits = data.fit

    ax[r,c].plot(x,y,'b--',label="model")
    ax[r,c].plot(x,fits(x),'k',label="fit")


    c += 1
    if(check%len(inc) == 0):
        r += 1
        c = 0
    check += 1
plt.legend()
plt.show()
