import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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

plt.plot(x, y)

z = getMeanData(directory,75)

x = z[:,7]

y = z[:,9]

plt.plot(x, y)

plt.axis([3000, 7200, 0, 13])
plt.show()