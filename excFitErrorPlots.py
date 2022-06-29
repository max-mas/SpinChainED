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

path = "/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/GapN"

fig, ax = plt.subplots()
for N in np.linspace(6, 16, 6):
    Js = []
    Gaps5 = []
    Gaps1 = []
    path5 = path + str(int(N)) + "It5"
    file5 = open(path5)
    lines5 = file5.readlines()
    for line in lines5:
        data = line.split(" ")
        Js.append(float(data[0]))
        Gaps5.append(float(data[1].replace("\n", "")))
    path1 = path + str(int(N)) + "It1"
    file1 = open(path1)
    lines1 = file1.readlines()
    for line in lines1:
        data = line.split(" ")
        Gaps1.append(float(data[1].replace("\n", "")))
    lab = "$N=$ " + str(int(N))
    ax.plot(Js, np.abs((np.asarray(Gaps5) - np.asarray(Gaps1))/np.asarray(Gaps5)), ".-", label=lab)
ax.semilogy()
ax.legend()
ax.set(xlabel="$J_1/J_2$", ylabel="$|\\left(\\Delta_{n=1}-\\Delta_{n=5}\\right)/\\Delta_{n=5}|$", title="Relative Differenz der angepassten Spinlücken \n für $n=1$ und $n=5$")
#fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/RelDiffs.png")
plt.show()

"""
Ns =[]
diffs1 = []
diffs5 = []
path = "/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/MeanRelDiffsIt1.txt"
file = open(path)
lines = file.readlines()
for line in lines:
    Data = line.split(" ")
    Ns.append(int(Data[0]))
    diffs1.append(float(Data[1].replace("\n", "")))

path = "/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/MeanRelDiffsIt5.txt"
file = open(path)
lines = file.readlines()
for line in lines:
    Data = line.split(" ")
    diffs5.append(float(Data[1].replace("\n", "")))

fig, ax = plt.subplots()
ax.set(xlabel="$N$", ylabel="$|\\overline\\delta_{n=1}-\\overline\\delta_{n=5}|$", title="Differenz der mittleren relativen Residuen der\n angepassten Spinlücke für $n=1$ und $n=5$")
ax.plot(Ns, np.abs(np.asarray(diffs1) - np.asarray(diffs5)), ".-")
fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/MeanDiffs.png")
#plt.show()
"""

