from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os
import natsort

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
    return beta* A * np.exp(-beta * B)


nMin = 6
nMax = 14
nNum = int((nMax - nMin) / 2) + 1

gapsQT = []

path = "/home/mmaschke/BA_Code/Data/out/Susceptibilities_DQT/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gapsQT.append([])
    gapsQT[i].append([])
    gapsQT[i].append([])
    for filePath in natsort.natsorted(os.listdir(path)):
        file = open(path + filePath)
        name1 = filePath.split("J")
        name2 = name1[1].split("I")
        name = name2[0]
        J_str = name.replace("_", ".")
        N_str = name1[0].split("N")
        N_test = int(N_str[1])
        if N_test != N:
            continue
        J = float(J_str)
        # print(J)

        betas = []
        Cs = []
        lines = file.readlines()
        for line in lines:
            data = line.split(" ")
            if float(data[0]) < 20:
                continue
            if data[1] == "nan\n" or data[1] == "-nan\n" or data[1] == "inf\n" or data[1] == "-inf\n":
                data[1] = "0"
            betas.append(float(data[0]))
            Cs.append(float(data[1].replace("\n", "")))

        try:
            parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
        except RuntimeError:
            parameters = [0, 0]
            covariance = [0, 0]
            print("Fit Error for N = " + str(N) + " and J = " + str(J))
        perr = np.sqrt(np.diag(covariance))
        gapsQT[i][0].append(J)
        gapsQT[i][1].append(parameters[1])
        print(str(N) + " " + J_str + " " + str(perr[1]))
    i += 1

fig, ax = plt.subplots()
N = nMin
for arr in gapsQT:
    lab = "$N$ = " + str(int(N)) + " QT Fit"
    ax.plot(arr[0], arr[1], ".-", label=lab, alpha=0.5)
    N += 2

# ax.set_xlim(0, 2)

gapsEQ = []

path = "/home/mmaschke/BA_Code/Data/out/Susceptibilities/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gapsEQ.append([])
    gapsEQ[i].append([])
    gapsEQ[i].append([])
    for filePath in natsort.natsorted(os.listdir(path)):
        file = open(path + filePath)
        name1 = filePath.split("J")
        name2 = name1[1].split(".")
        name = name2[0]
        J_str = name.replace("_", ".")
        N_str = name1[0].split("N")
        N_test = int(N_str[1])
        if N_test != N:
            continue
        J = float(J_str)

        betas = []
        Cs = []
        lines = file.readlines()
        for line in lines:
            data = line.split(" ")

            if float(data[0]) < 15:
                continue
            if data[1] == "nan\n" or data[1] == "-nan\n" or data[1] == "inf\n" or data[1] == "-inf\n":
                data[1] = "0"
            betas.append(float(data[0]))
            Cs.append(float(data[1].replace("\n", "")))

        try:
            parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
        except RuntimeError:
            parameters = [0, 0]
            covariance = [0, 0]
            print("Fit Error for N = " + str(N) + " and J = " + str(J))
            # plt.plot(betas, Cs)
            # plt.show()
        perr = np.sqrt(np.diag(covariance))
        gapsEQ[i][0].append(J)
        gapsEQ[i][1].append(parameters[1])
        print(str(N) + " " + J_str + " " + str(perr[1]))
    i += 1

N = nMin
plt.gca().set_prop_cycle(None)
for arr in gapsEQ:
    lab = "$N$ = " + str(int(N)) + " ED Fit"
    ax.plot(arr[0], arr[1], ".-", label=lab, alpha=0.2)
    N += 2

plt.gca().set_prop_cycle(None)
for N in np.linspace(nMin, nMax, nNum):
    path = "/home/mmaschke/BA_Code/Data/out/ExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
    file = open(path)
    lines = file.readlines()
    Ts = []
    excErg = []
    for line in lines:
        data = line.split(" ")
        Ts.append(float(data[0]))
        excErg.append(float(data[1].replace("\n", "")))
    lab = "$N$ = " + str(int(N)) + " exact"
    ax.plot(Ts, excErg, label=lab)

ax.set(xlabel="$J_1/J_2$", ylabel="Spin Gap Energy $\\Delta$ ($J_2$)")
ax.legend()
ax.set_ylim(0, 2)
ax.set_xlim(0, 2)
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
