from posixpath import join
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sci
from scipy import interpolate as interp
import math
import os
import functions as f
import tikzplotlib as tikz

np.set_printoptions(precision=4)


def getMeanData(directory, throttleSettingValue):

    throttleSetting = 'thr' + str(throttleSettingValue)

    fileList = os.listdir(directory)

    meanArray = np.empty([sum(throttleSetting in s for s in fileList), 21])

    i = 0
    for file in fileList:
        if throttleSetting in file:
            data = pd.read_csv(os.path.join(
                directory, file), sep=';', decimal=',')
            meanArray[i, :] = data.mean().to_numpy()
            i += 1

    return meanArray


directory = './data/2022-11-16_Messdaten'

fig,ax = plt.subplots()
xx = []
yy = []
zz = []

xc = []
yc = []

z = getMeanData(directory, 100)

x = z[:, 9]

y = z[:, 11]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

print(y)
print(np.amax(x), np.amax(y))


ax.plot(x, y,'k')
ax.plot(np.amax(x), np.amax(y), '.r', markersize=10)
ax.text(np.amin(x)-340,np.amin(y)-0.15,'100\\%',fontsize = 14)

print(np.amax(x),np.amax(y),np.amin(z[:,10]))


z = getMeanData(directory, 95)

x = z[:, 9]

y = z[:, 11]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amin(y)-0.15,'95\\%',fontsize = 14)


#95%
z = getMeanData(directory, 90)

x = z[:, 9]

y = z[:, 11]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amin(y)-0.15,'90\\%',fontsize = 14)

#85%
z = getMeanData(directory, 85)

x = z[:, 9]

y = z[:, 11]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amin(y)-0.15,'85\\%',fontsize = 14)

#75%
z = getMeanData(directory, 75)

x = z[:, 9]

y = z[:, 11]

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amin(y)-0.15,'75\\%',fontsize = 14)

#50%
z = getMeanData(directory, 50)

x = z[:, 9]

y = z[:, 11]


xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amin(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (y*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)
print(zz)
ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amin(y)-0.15,'50\\%',fontsize = 14)



points = [xx, yy]
points = np.transpose(points)
values = zz

griddim = 100j
x_grid, y_grid = np.mgrid[1800:7200:griddim, 2:12.7:griddim]

grid1 = interp.griddata(points, values, (x_grid, y_grid), method='cubic',rescale=1)

i = 0
while i < grid1.shape[0]:
    j = 0
    while j < grid1.shape[1]:
    
        if (x_grid[i,j] <= xc[4] or x_grid[i,j] <= xc[5]) and (y_grid[i,j] >= yc[5] or y_grid[i,j] <= yc[4]) and not math.isnan(grid1[i,j]):
            v1 = (xc[4]-xc[5],yc[4]-yc[5])
            vc = (x_grid[i,j]-xc[5],y_grid[i,j]-yc[5])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')
        
        if (x_grid[i,j] <= xc[3] or x_grid[i,j] <= xc[4]) and (y_grid[i,j] >= yc[4] and y_grid[i,j] <= yc[3]) and not math.isnan(grid1[i,j]):
            v1 = (xc[3]-xc[4],yc[3]-yc[4])
            vc = (x_grid[i,j]-xc[4],y_grid[i,j]-yc[4])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')

        if (x_grid[i,j] <= xc[3] or x_grid[i,j] <= xc[0]) and (y_grid[i,j] >= yc[3] and y_grid[i,j] <= yc[0]) and not math.isnan(grid1[i,j]):
            v1 = (xc[0]-xc[3],yc[0]-yc[3])
            vc = (x_grid[i,j]-xc[3],y_grid[i,j]-yc[3])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')

        if (x_grid[i,j] <= xc[0] or x_grid[i,j] <= xc[1]) and (y_grid[i,j] >= yc[0] and y_grid[i,j] <= yc[1]) and not math.isnan(grid1[i,j]):
            v1 = (xc[1]-xc[0],yc[1]-yc[0])
            vc = (x_grid[i,j]-xc[0],y_grid[i,j]-yc[0])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')



        j += 1


    i += 1
    f.printProgressBar(i,np.imag(griddim),'Masking:',)
    
Xmax, Ymax = np.unravel_index(np.nanargmax(grid1), np.shape(grid1))

CS = ax.contour(x_grid, y_grid, grid1, 20,
                 linestyles='dashed', linewidths=1, cmap='Wistia')

#ax.clabel(CS, CS.levels, inline=True, fontsize=10)


#ax.plot(x_grid[Xmax, Ymax], y_grid[Xmax, Ymax], 'xr', markersize=5)
ax.clabel(CS, inline=1,inline_spacing = -7, fontsize=14,fmt="%1.1f")
ax.text(x_grid[Xmax, Ymax]+30, y_grid[Xmax, Ymax] +
         0.2, "%.2f" % np.nanmax(grid1),fontsize=14)

ax.grid(1,'both',ls='--')
ax.set_xlim((1600,7000))
ax.set_ylim((0,13.2))
plt.xlabel('n\\textsubscript{2} [RPM]')
plt.ylabel('P [kW]')


tikz.save("D:/Users/Fabian Zach/OneDrive/Uni/Masterarbeit/Latex/Arbeit/Ressourcen/Bilder/5/kennfeld_eff_2.tex",axis_height = '\\figH', axis_width = '\\figW')

plt.savefig('D:/Users/Fabian Zach/OneDrive/Uni/Masterarbeit/Latex/Arbeit/Ressourcen/Bilder/5/kennfeld_eff_2.pgf', format = 'pgf')
plt.show()
plt.close()

fig,ax = plt.subplots()
xx = []
yy = []
zz = []

xc = []
yc = []

z = getMeanData(directory, 100)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-330,np.amax(y)-0.15,'100\\%',fontsize = 14)


z = getMeanData(directory, 95)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amax(y)-0.15,'95\\%',fontsize = 14)

    


#95%
z = getMeanData(directory, 90)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)


ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amax(y)-0.15,'90\\%',fontsize = 14)



#85%
z = getMeanData(directory, 85)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amax(y)-0.15,'85\\%',fontsize = 14)



#75%
z = getMeanData(directory, 75)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amax(y)-0.15,'75\\%',fontsize = 14)



#50%
z = getMeanData(directory, 50)

x = z[:, 9]

y = z[:, 10]

xc = np.append(xc,np.amin(x))
yc = np.append(yc,np.amax(y))

xx = np.append(xx, x, 0)
yy = np.append(yy, y, 0)
zz = np.append(zz, (z[:,11]*1e3)/(z[:, 2]*43.1e3*300/350)*100, 0)

ax.plot(x, y,'k')
ax.text(np.amin(x)-270,np.amax(y)-0.15,'50\\%',fontsize = 14)

points = [xx, yy]
points = np.transpose(points)
values = zz

griddim = 100j
x_grid, y_grid = np.mgrid[1800:7200:griddim, 2:24:griddim]

grid1 = interp.griddata(points, values, (x_grid, y_grid), method='cubic',rescale=1)
print(xc,yc)
i = 0
while i < grid1.shape[0]:
    j = 0
    while j < grid1.shape[1]:
    
        if (x_grid[i,j] <= xc[4] or x_grid[i,j] <= xc[5]) and (y_grid[i,j] >= yc[5] or y_grid[i,j] <= yc[4]) and not math.isnan(grid1[i,j]):
            v1 = (xc[4]-xc[5],yc[4]-yc[5])
            vc = (x_grid[i,j]-xc[5],y_grid[i,j]-yc[5])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')
        
        if (x_grid[i,j] <= xc[3] or x_grid[i,j] <= xc[4]) and (y_grid[i,j] >= yc[4] and y_grid[i,j] <= yc[3]) and not math.isnan(grid1[i,j]):
            v1 = (xc[3]-xc[4],yc[3]-yc[4])
            vc = (x_grid[i,j]-xc[4],y_grid[i,j]-yc[4])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')

        if (x_grid[i,j] <= xc[3] or x_grid[i,j] <= xc[0]) and (y_grid[i,j] >= yc[3] and y_grid[i,j] <= yc[0]) and not math.isnan(grid1[i,j]):
            v1 = (xc[0]-xc[3],yc[0]-yc[3])
            vc = (x_grid[i,j]-xc[3],y_grid[i,j]-yc[3])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')

        if (x_grid[i,j] <= xc[0] or x_grid[i,j] <= xc[1]) and (y_grid[i,j] >= yc[0] and y_grid[i,j] <= yc[1]) and not math.isnan(grid1[i,j]):
            v1 = (xc[1]-xc[0],yc[1]-yc[0])
            vc = (x_grid[i,j]-xc[0],y_grid[i,j]-yc[0])

            cp = v1[0]*vc[1] - v1[1]*vc[0]
            if cp > 0:
                grid1[i,j] = float('nan')



        j += 1


    i += 1
    f.printProgressBar(i,np.imag(griddim),'Masking:',)
    
Xmax, Ymax = np.unravel_index(np.nanargmax(grid1), np.shape(grid1))

CS = ax.contour(x_grid, y_grid, grid1, 20,
                 linestyles='dashed', linewidths=1, cmap='Wistia')

#ax.clabel(CS, CS.levels, inline=True, fontsize=10)


#ax.plot(x_grid[Xmax, Ymax], y_grid[Xmax, Ymax], 'xr', markersize=5)
ax.clabel(CS, inline=1,inline_spacing = -7, fontsize=14,fmt="%1.1f")
ax.grid(1,'both',ls='--')
ax.set_xlim((1600,7000))
ax.set_ylim((0,25))

plt.xlabel('n\\textsubscript{2} [RPM]')
plt.ylabel('M [Nm]')

tikz.save("D:/Users/Fabian Zach/OneDrive/Uni/Masterarbeit/Latex/Arbeit/Ressourcen/Bilder/5/kennfeld_M_2.tex",axis_height = '\\figH', axis_width = '\\figW')

plt.show()

