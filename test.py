from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os
import natsort

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=20)

path = "/home/mmaschke/BA_Code/Data/out/SpecificHeats/forFit/SpecHeatN10J0_448980.txt"
file = open(path)

betas = []
Cs = []
lines = file.readlines()
for line in lines:
    data = line.split(" ")

    if float(data[0]) < 0:
        continue
    if data[1] == "nan\n" or data[1] == "-nan\n" or data[1] == "inf\n" or data[1] == "-inf\n":
        data[1] = "0"
    betas.append(float(data[0]))
    Cs.append(float(data[1].replace("\n", "")))

plt.plot(betas, Cs)
plt.show()
