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
    return a2 - (a2 - a1) ** 2 / ((a2 - a1) - (a1 - a0))


def epsilon(a):
    N = len(a)
    eps = np.zeros((N, N+1))
    for l in range(N):
        eps[l][1] = a[l]
    for k in range(N+2):
        if k < 2:
            continue
        for l in range(N):
            if l >= N-(k-1):
                continue
            eps[l][k] = eps[l + 1][k - 2] + 1 / (eps[l + 1][k - 1] - eps[l][k - 1])
    if N % 2 == 0:
        return eps[0][N-1]
    else:
        return eps[0][N]


def euler(a1, a2):
    return ()


def weird_transform(Js, Vals):
    A = copy.deepcopy(Vals)

    for y in range(len(Js)):
        A[y] *= 1 / (1 + Js[y])
    return A


def lin(x, a, b):
    return a*x + b


nMin = 6
nMax = 18
nNum = 7  # int((nMax - nMin) / 2) + 1
nNumLow = 4
dataPointNum = 200

fullGaps = []
gapsLow = []
for N in [6, 8, 10, 12, 14, 16, 18]:
    path1 = "/home/mmaschke/BA_Code/Data/out/ExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
    file1 = open(path1, "r")
    lines1 = file1.readlines()
    fullJs = []
    JsLow = []
    gapsNLow = []
    fullGapsN = []
    for line1 in lines1:
        data1 = line1.split(" ")
        fullJs.append(float(data1[0]))
        fullGapsN.append(float(data1[1].replace("\n", "")))
        if float(data1[0]) < 0.65:
            JsLow.append(float(data1[0]))
            gapsNLow.append(float(data1[1].replace("\n", "")))
    gapsLow.append(gapsNLow)
    fullGaps.append(fullGapsN)

gapsHigh = []
for N in [6, 10, 14, 18]:
    path2 = "/home/mmaschke/BA_Code/Data/out/ExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
    if N == 22:
        path2 = "/home/mmaschke/BA_Code/Data/out/GapFit/spin/gapsLowJ" + str(int(N)) + ".txt"
    file2 = open(path2, "r")
    lines2 = file2.readlines()
    JsHigh = []
    gapsNHigh = []
    for line2 in lines2:
        data2 = line2.split(" ")
        if float(data2[0]) >= 0.65:
            data2 = line2.split(" ")
            JsHigh.append(float(data2[0]))
            gapsNHigh.append(float(data2[1].replace("\n", "")))
    gapsHigh.append(gapsNHigh)

gaps_extrapLow = []
gaps_extrapHigh = []
for j in range(len(gapsLow[0])):
    gapsNLow = []
    for i in range(nNum):
        gapsNLow.append(gapsLow[i][j])
    gaps_extrapLow.append(epsilon(gapsNLow))
for j in range(len(gapsHigh[0])):
    gapsNHigh = []
    for i in range(nNum):
        if i < nNumLow:
            gapsNHigh.append(gapsHigh[i][j])
    gaps_extrapHigh.append(epsilon(gapsNHigh))

fig, ax = plt.subplots()
N = 6
for gap in fullGaps:
    ax.plot(fullJs, weird_transform(fullJs, gap), label="ED, $N=$" + str(N))
    N += 2
ax.plot(JsLow, weird_transform(JsLow, gaps_extrapLow), "r--", label="QT Extrapolation Low")
ax.plot(JsHigh, weird_transform(JsHigh, gaps_extrapHigh), "g--", label="QT Extrapolation High")
ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Spin Gap Energy $\\Delta/(J_1+J_2)$", title="QT Data")
ax.set_ylim(0, 0.8)
ax.set_xlim(0, 2)
ax.legend()
#plt.show()
plt.close(fig)


offsets = []
gaps = gapsLow #+ gapsHigh
Ns = np.linspace(nMin, nMax, nNum)
NsHigh = [6, 10, 14, 18]
RecipNsPlot = np.linspace(0.001, 0.5, 200)
for j in range(dataPointNum):
    fig, ax = plt.subplots()
    jGaps = []
    for gap in fullGaps:
        jGaps.append(gap[j] / (1 + fullJs[j]))
    jGapsHigh = []
    for gap in gapsHigh:
        jGapsHigh.append(gap[j-len(gapsLow[0])] / (1 + fullJs[j]))
    if fullJs[j] >= 0.65:
        parameters, covariance = scipy.optimize.curve_fit(lin, 1/np.asarray(NsHigh), jGapsHigh)
    else:
        parameters, covariance = scipy.optimize.curve_fit(lin, 1/np.asarray(Ns), jGaps)
    offsets.append(parameters[1])
    ax.plot(1/np.asarray(Ns), jGaps, ".-")
    ax.plot(RecipNsPlot, lin(RecipNsPlot, parameters[0], parameters[1]), "--")
    ax.set(xlabel="$1/N$", ylabel="Reduced Spin Gap Energy $\\Delta/(J_1+J_2)$", title="ED Fit, $J_1/J_2=$" + str(fullJs[j]))
    ax.set_ylim(0, 0.8)
    ax.set_xlim(0, 0.5)
    #ax.semilogy()
    #ax.set_xlim(0, 2)
    fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Extrap/Single_pointsED/J" + str(fullJs[j]).replace(".", "_") + ".png")
    #plt.show()
    plt.close(fig)

fig, ax = plt.subplots()
N = 6
for gap in fullGaps:
    ax.plot(fullJs, weird_transform(fullJs, gap), label="ED, $N=$" + str(N))
    N += 2
ax.plot(fullJs, weird_transform(fullJs, offsets), "r--", label="ED Fit-Extrapolation")
ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Spin Gap Energy $\\Delta/(J_1+J_2)$", title="ED Data")
#ax.set_ylim(0, 0.8)
ax.set_xlim(0, 2)
ax.legend()
plt.show()
