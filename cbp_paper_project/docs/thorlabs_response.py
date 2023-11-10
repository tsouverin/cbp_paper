import matplotlib.pyplot as plt
plt.interactive(True)

import numpy as np
import pandas as pd

import csv

# Old one
SM05PD1B = pd.read_excel('FDS100_Res_data.xlsx', usecols='C,D',names=['w', 'R'])

# New one
SM05PD3A = pd.read_excel('FD11A_data.xlsx', usecols='C,D',names=['w', 'R'])

plt.figure(1)
plt.clf()
plt.plot(SM05PD1B.w[1:], SM05PD1B.R[1:], label="old/SM05PD1B")
plt.plot(SM05PD3A.w[1:], SM05PD3A.R[1:], label="new/SM05PD3A")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Response (A/W)")
plt.legend()
plt.savefig("thorlabs_response.png")



ga = np.rec.fromarrays([SM05PD3A.w[1:].astype(float), SM05PD3A.R[1:].astype(float)], names=['wavelength', 'Response'])
np.save("thorlabs_SM05PD3A.npy", ga)

bu = np.rec.fromarrays([SM05PD1B.w[1:-1].astype(float), SM05PD1B.R[1:-1].astype(float)], names=['wavelength', 'Response'])
np.save("thorlabs_SM05PD1B.npy", bu)


fh = open('thorlabs_SM05PD3A.csv', 'w')
print(fh, '# Wavelength (nm), Responsivity (A/W)')
w = csv.writer(fh)
w.writerows(ga)
fh.close()
