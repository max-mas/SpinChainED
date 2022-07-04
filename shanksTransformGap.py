import scipy.interpolate
from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os
import natsort
import copy

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=24)
plt.rc('axes', titlesize=30)
plt.rcParams["figure.figsize"] = (12, 9)


def shanks(a0, a1, a2):
    return a2 - (a2 - a1)**2/((a2 - a1) - (a1 - a0))


def weird_transform(Js, Vals):
    A = copy.deepcopy(Vals)

    for y in range(len(Js)):
        A[y] *= 1 / (1 + Js[y])
    return A


nMin = 6
nMax = 18
nNum = 4#int((nMax - nMin) / 2) + 1
dataPointNum = 50

gaps = []

for N in [6, 10, 14, 18]:
    path = "/home/mmaschke/BA_Code/Data/out/GapFit/spin/gapsHighJ" + str(int(N)) + ".txt"
    file = open(path, "r")
    lines = file.readlines()
    Js = []
    gapsN = []
    for line in lines:
        data = line.split(" ")
        Js.append(float(data[0]))
        gapsN.append(float(data[1].replace("\n", "")))
    gaps.append(gapsN)

gaps_shanks = []
gaps_plot = []
for j in range(dataPointNum):
    gaps_shanksJ = []
    for i in range(nNum):
        if i < 1 or i >= nNum - 1:
            continue
        gaps_shanksJ.append(shanks(gaps[i-1][j], gaps[i][j], gaps[i+1][j]))
    gaps_shanks.append(gaps_shanksJ)
    gaps_plot.append(gaps_shanksJ[-1])
"""
gaps_shanks2 = []
for j in range(dataPointNum):
    gaps_shanksJ = []
    for i in range(nNum):
        if i < 1 or i >= nNum - 3:
            continue
        gaps_shanksJ.append(shanks(gaps_shanks[j][i-1], gaps_shanks[j][i], gaps_shanks[j][i+1]))
    gaps_shanks2.append(gaps_shanksJ)

gaps_shanks3 = []
for j in range(dataPointNum):
    gaps_shanksJ = []
    for i in range(nNum):
        if i < 1 or i >= nNum - 5:
            continue
        gaps_shanksJ.append(shanks(gaps_shanks2[j][i-1], gaps_shanks2[j][i], gaps_shanks2[j][i+1]))
    gaps_shanks3.append(gaps_shanksJ)
    gaps_plot.append(gaps_shanksJ[-1])
"""
fig, ax = plt.subplots()
ax.scatter(Js, weird_transform(Js, gaps_plot))
ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Extrapolated Spin Gap Energy $\\Delta/(J_1+J_2)$", title="QT Data")
#ax.set_ylim(0, 0.8)
ax.set_xlim(0, 2)
plt.show()


