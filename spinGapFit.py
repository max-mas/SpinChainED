from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=20)


def read_vec_from_file(arg_path):
    vals = []
    arg_file = open(arg_path)
    arg_lines = arg_file.readlines()
    for arg_line in arg_lines:
        if float.is_integer(float(arg_line)):
            stri = arg_line.replace("\n", "") + "_000000"
        else:
            stri = arg_line.replace(".", "_")
            stri = stri.replace("\n", "") + "00000"
        vals.append(stri)
    return vals


def exp_fit_fn(beta, A, B):
    return beta**2 * A * np.exp(-beta * B)

nMin = 6
nMax = 12
nNum = int((nMax-nMin)/2) + 1
"""
gaps = []

path = "D:/Code/C++/spinChainData/out/SpecificHeats_DQT/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gaps.append([])
    gaps[i].append([])
    gaps[i].append([])
    for filePath in os.listdir(path):
        file = open(path+filePath)
        name1 = filePath.split("J")
        name2 = name1[1].split("I")
        name = name2[0]
        J_str = name.replace("_", ".")
        N_str = name1[0].split("N")
        N_test = int(N_str[1])
        if N_test != N:
            continue
        J = float(J_str)
        #print(J)

        betas = []
        Cs = []
        lines = file.readlines()
        for line in lines:
            data = line.split(" ")
            if float(data[0]) < 25:
                continue
            betas.append(float(data[0]))
            Cs.append(float(data[1].replace("\n", "")))

        try:
            parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
        except RuntimeError:
            parameters = [0, 0]
            print("Fit Error for N = " + str(N) + " and J = " + str(J))
        gaps[i][0].append(J)
        gaps[i][1].append(parameters[1])
    i += 1

fig, ax = plt.subplots()
N = nMin
for arr in gaps:
    lab = "$N$ = " + str(int(N))
    ax.scatter(arr[0], arr[1], label=lab)
    N += 2

ax.set(xlabel="$J_1/J_2$", ylabel="Spin Gap Energy $\\Delta E$ from DQT Fit ($J_2$)")
ax.legend()
ax.set_ylim(0, 2)
#ax.set_xlim(0, 2)
plt.show()
"""
gaps = []

path = "D:/Code/C++/spinChainData/out/SpecificHeats/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gaps.append([])
    gaps[i].append([])
    gaps[i].append([])
    for filePath in os.listdir(path):
        file = open(path+filePath)
        name1 = filePath.split("J")
        name2 = name1[1].split(".")
        name = name2[0]
        J_str = name.replace("_", ".")
        N_str = name1[0].split("N")
        N_test = int(N_str[1])
        if N_test != N:
            continue
        J = float(J_str)
        #print(J)

        betas = []
        Cs = []
        lines = file.readlines()
        for line in lines:
            data = line.split(" ")
            if float(data[0]) < 25:
                continue
            betas.append(float(data[0]))
            Cs.append(float(data[1].replace("\n", "")))

        try:
            parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
        except RuntimeError:
            parameters = [0, 0]
            print("Fit Error for N = " + str(N) + " and J = " + str(J))
        gaps[i][0].append(J)
        gaps[i][1].append(parameters[1])
    i += 1

fig, ax = plt.subplots()
N = nMin
for arr in gaps:
    lab = "$N$ = " + str(int(N))
    ax.plot(arr[0], arr[1], label=lab)
    N += 2

ax.set(xlabel="$J_1/J_2$", ylabel="Spin Gap Energy $\\Delta E$ from ED Fit ($J_2$)")
ax.legend()
ax.set_ylim(0, 2)
#ax.set_xlim(0, 2)
plt.show()




"""
file = open(path)
lines = file.readlines()
betas = []
Cs = []
for line in lines:
    data = line.split(" ")
    if float(data[0]) < 33:
        continue
    betas.append(float(data[0]))
    Cs.append(float(data[1].replace("\n", "")))

parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
print(parameters)
print(covariance)
"""