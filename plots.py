import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

np.set_printoptions(precision=2)

def getMeanData(directory,nStart,nEnd):

    meanArray = np.zeros((nEnd-nStart+1,13))

    n = nStart
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.startswith(str(n).zfill(2) or str(n)):
            data = pd.read_csv(os.path.join(directory,filename), sep=';', decimal=',')
            meanArray[n-nStart,:] = data.mean().to_numpy()
            n += 1

        if n > nEnd:
            print('done')
            break

    return meanArray

directory = './data/2022-08-25_Messdaten'

z = getMeanData(directory,2,12)

plt.close()

x = z[:,7]

y = z[:,9]

plt.plot(x, y)

z = getMeanData(directory,13,21)

x = z[:,7]

y = z[:,9]

plt.plot(x, y)
plt.axis([3000, 7000, 0, 13])
plt.show()