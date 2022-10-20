import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sci
from scipy import interpolate as interp
import math
import os
import functions as f




np.set_printoptions(precision=4)

def getMeanData(directory,throttleSettingValue):

    throttleSetting = 'thr' + str(throttleSettingValue)

    fileList = os.listdir(directory)

    meanArray = np.empty([sum(throttleSetting in s for s in fileList),13])

    i = 0
    for file in fileList:

        if throttleSetting in file:
            data = pd.read_csv(os.path.join(directory,file), sep=';', decimal=',')
            meanArray[i,:] = data.mean().to_numpy()
            i += 1

    return meanArray

directory = './data/2022-08-25_Messdaten'

z = getMeanData(directory,100)

plt.close()

x = z[:,7]

y = z[:,9]

xx = x
yy = y

zz = (y*1e3)/(z[:,2]*43e3*300/350)

# plt.plot(x, y)

z = getMeanData(directory,75)

x = z[:,7]

y = z[:,9]

xx = np.append(xx,x,0)
yy = np.append(yy,y,0)
zz = np.append(zz,(y*1e3)/(z[:,2]*43e3*300/350),0)



# plt.plot(x, y)

# plt.axis([3000, 7200, 0, 13])

# plt.show()

points = [xx,yy]
points = np.transpose(points)
values = zz

x_grid,y_grid = np.mgrid[3000:7200:1000j,0:13:1000j]

plt.plot(points[0:10,0],points[0:10,1])
plt.plot(points[11:20,0],points[11:20,1])

grid1 = interp.griddata(points,values,(x_grid,y_grid),method='linear')

Xmax, Ymax = np.unravel_index(np.nanargmax(grid1),np.shape(grid1))

print(x_grid[Xmax,Ymax],y_grid[Xmax,Ymax])

CS = plt.contour(x_grid,y_grid,grid1,9, linestyles = 'dashed', linewidths = 0.5,cmap = 'Reds')
plt.plot(x_grid[Xmax,Ymax],y_grid[Xmax,Ymax],'Dr', markersize =5)
plt.clabel(CS, inline=1, fontsize=10)
plt.text(x_grid[Xmax,Ymax]+30,y_grid[Xmax,Ymax]+0.2,"%.3f" % np.nanmax(grid1))
plt.show()
