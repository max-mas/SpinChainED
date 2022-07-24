import scipy.interpolate
from matplotlib import pyplot as plt
import matplotlib.lines as lines
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os
import natsort
import copy

plt.rcParams['text.usetex'] = True
plt.rc('font', family='serif')
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=20)


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


def quad(x, a, b, c):
    return a*x**2 + b*x + c


nMin = 6
nMax = 18
nNum = 7  # int((nMax - nMin) / 2) + 1
nNumLow = 4
dataPointNum = 200
ED = True
lowCutoff = 0.65
fitRangeLow = 0.4
fitRangeHigh = 0.8

fullGaps = []
fullErrs = []
gapsLow = []
for N in [6, 8, 10, 12, 14, 16, 18]:
    path1 = "D:/Code/C++/spinChainData/ba_data_2207/Data/out/ExcitationErgs/Excergs" + str(int(N)) + ".txt"
    file1 = open(path1, "r")
    lines1 = file1.readlines()
    fullJs = []
    JsLow = []
    gapsNLow = []
    fullGapsN = []
    fullErrsN = []
    for line1 in lines1:
        data1 = line1.split(" ")
        fullJs.append(float(data1[0]))
        fullGapsN.append(float(data1[1].replace("\n", "")))
        if ED:
            fullErrsN.append(0)
        else:
            fullErrsN.append(float(data1[2].replace("\n", "")))
        if float(data1[0]) < lowCutoff:
            JsLow.append(float(data1[0]))
            gapsNLow.append(float(data1[1].replace("\n", "")))
    gapsLow.append(gapsNLow)
    fullGaps.append(fullGapsN)
    fullErrs.append(fullErrsN)

gapsHigh = []
gapErrsHigh = []
for N in [6, 10, 14, 18]:
    path2 = "D:/Code/C++/spinChainData/ba_data_2207/Data/out/ExcitationErgs/Excergs" + str(int(N)) + ".txt"
    #if N == 22:
    #    path2 = "/home/mmaschke/BA_Code/Data/out/GapFit/spin/gapsIt20lowJ" + str(int(N)) + ".txt"
    file2 = open(path2, "r")
    lines2 = file2.readlines()
    JsHigh = []
    gapsNHigh = []
    gapErrsNHigh = []
    for line2 in lines2:
        data2 = line2.split(" ")
        if float(data2[0]) >= lowCutoff:
            data2 = line2.split(" ")
            JsHigh.append(float(data2[0]))
            gapsNHigh.append(float(data2[1].replace("\n", "")))
        if ED:
            gapErrsNHigh.append(0)
        else:
            gapErrsNHigh.append(float(data1[2].replace("\n", "")))
    gapsHigh.append(gapsNHigh)
    gapErrsHigh.append(gapErrsNHigh)

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
ax.plot(JsLow, weird_transform(JsLow, gaps_extrapLow), "r--", label="ED Extrapolation Low")
ax.plot(JsHigh, weird_transform(JsHigh, gaps_extrapHigh), "g--", label="ED Extrapolation High")
ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Spin Gap Energy $\\Delta/(J_1+J_2)$", title="ED Data")
ax.set_ylim(0, 0.8)
ax.set_xlim(0, 2)
ax.legend()
#plt.show()
plt.close(fig)


offsets = []
offsetErrs = []
gaps = gapsLow #+ gapsHigh
Ns = np.linspace(nMin, nMax, nNum)
NsHigh = [6, 10, 14, 18]
RecipNsPlot = np.linspace(0.001, 0.5, 200)
for j in range(dataPointNum):
    fig, ax = plt.subplots()
    jGaps = []
    jErrs = []
    for gapErr in zip(fullGaps, fullErrs):
        jGaps.append(gapErr[0][j])
        jErrs.append(gapErr[1][j])
    jGapsHigh = []
    jErrsHigh = []
    for gapErr in zip(gapsHigh, gapErrsHigh):
        jGapsHigh.append(gapErr[0][j-len(gapsLow[0])])
        jErrsHigh.append(gapErr[1][j])
    if ED:
        if fullJs[j] >= lowCutoff:
            parameters, covariance = scipy.optimize.curve_fit(quad, 1/np.asarray(NsHigh), jGapsHigh, bounds=[[-np.inf, -np.inf, -np.inf], [1, np.inf, np.inf]])
        else:
            parameters, covariance = scipy.optimize.curve_fit(quad, 1/np.asarray(Ns), jGaps, bounds=[[-np.inf, -np.inf, -np.inf], [0.5, np.inf, np.inf]])
    else:
        if fullJs[j] >= lowCutoff:
            parameters, covariance = scipy.optimize.curve_fit(quad, 1/np.asarray(NsHigh), jGapsHigh, sigma=jErrsHigh, absolute_sigma=False)
        else:
            parameters, covariance = scipy.optimize.curve_fit(quad, 1/np.asarray(Ns), jGaps, sigma=jErrs, absolute_sigma=False)
    offsets.append(parameters[2])
    offsetErrs.append(np.sqrt(covariance[2][2]))
    ax.errorbar(1/np.asarray(Ns), np.asarray(jGaps)/(1+fullJs[j]), fmt=".-", xerr=None, yerr=np.asarray(jErrs)/(1+fullJs[j]), capsize=2)
    ax.plot(RecipNsPlot, np.asarray(quad(RecipNsPlot, parameters[0], parameters[1], parameters[2]))/(1+fullJs[j]), "--")
    ax.set(xlabel="$1/N$", ylabel="Reduzierte Spinlücke $\\Delta/(J_1+J_2)$")
    ax.set_ylim(0, 0.8)
    ax.set_xlim(0, 0.2)
    #ax.semilogy()
    #ax.set_xlim(0, 2)
    #fig.savefig("D:/Code/C++/spinChainData/ba_data_2207/Data/plots/GapFit/spin/Extrap/Single_pointsED/J" + str(fullJs[j]).replace(".", "_") + ".pdf")
    #plt.show()
    plt.close(fig)

fig, ax = plt.subplots()
N = 6
for gapErr in zip(fullGaps, fullErrs):
    if ED:
        ax.plot(fullJs, weird_transform(fullJs, gapErr[0]), label="$N=$" + str(N))
    else:
        ax.errorbar(fullJs, weird_transform(fullJs, gapErr[0]), yerr=weird_transform(fullJs, gapErr[1]), xerr=None, label="$N=$" + str(N),
                    fmt=".-", capsize=2, lw=1)
    N += 2
ax.plot(fullJs, weird_transform(fullJs, offsets), "r--", label="Fit-Extrapolation", lw=1)
ax.fill_between(fullJs,  weird_transform(fullJs, np.asarray(offsets) - np.asarray(offsetErrs)),  weird_transform(fullJs, np.asarray(offsets) + np.asarray(offsetErrs)), color="r", alpha=0.1)
ax.set(xlabel="$j$", ylabel="Reduzierte Spinlücke $\\Delta/(J_1+J_2)$")
ax.set_ylim(-0.05, 0.75)
ax.set_xlim(0, 1.25)

ax.axhline(0, 0, 1, color="grey", ls="--", alpha=0.3, lw=0.8)

ax.legend(fontsize=5)
plt.show()
