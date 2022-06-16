from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
from matplotlib.ticker import MaxNLocator
import os
import natsort

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=24)
plt.rc('axes', titlesize=30)


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
nNum = int((nMax - nMin) / 2) + 1
numOfRuns = 1

gapsQT = []
gapsQTavg = []
gapsQTdev = []

path = "/home/mmaschke/BA_Code/Data/out/SpecificHeats_DQT/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gapsQT.append([])
    j = 0
    for k in range(1, numOfRuns + 1):

        runPath = path + str(k) + "/"

        gapsQT[i].append([])
        gapsQT[i][k-1].append([])
        gapsQT[i][k-1].append([])
        for filePath in natsort.natsorted(os.listdir(runPath)):
            file = open(runPath + filePath)
            name1 = filePath.split("J")
            name2 = name1[1].split("I")
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
            gapsQT[i][k-1][0].append(J)
            gapsQT[i][k-1][1].append(parameters[1])
            print(str(N) + " " + J_str + " " + str(perr[1]))

        j += 1
    gapsQTavg.append([])
    gapsQTavg[i].append(gapsQT[i][0][0])
    gapsQTavg[i].append([])

    gapsQTdev.append([])
    gapsQTdev[i].append(gapsQT[i][0][0])
    gapsQTdev[i].append([])

    for n in range(len(gapsQT[i][0][0])):
        avg = 0.0
        dev = 0.0
        for l in range(numOfRuns):
            avg += gapsQT[i][l][1][n]/numOfRuns
        for l in range(numOfRuns):
            dev += (gapsQT[i][l][1][n] - avg)**2 / numOfRuns
        dev = np.sqrt(dev)
        gapsQTavg[i][1].append(avg)
        gapsQTdev[i][1].append(dev)

    i += 1

fig, ax = plt.subplots()
N = nMin
i = 0
for arr in gapsQTavg:
    lab = "$N$ = " + str(int(N)) + " QT Fit"
    ax.errorbar(arr[0], arr[1], yerr=gapsQTdev[i][1], xerr=None, fmt=".-",  label=lab, alpha=0.5, capsize=2)
    N += 2
    i += 1
"""
gapsEQ = []

path = "/home/mmaschke/BA_Code/Data/out/SpecificHeats/forFit/"
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
"""
plt.gca().set_prop_cycle(None)
for N in np.linspace(nMin, nMax, nNum):
    path = "/home/mmaschke/BA_Code/Data/out/ActualExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
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

ax.set(xlabel="$J_1/J_2$", ylabel="Excitation Energy $\\Delta$ ($J_2$)", title="$n = $ " + str(numOfRuns))
ax.legend(prop={'size': 16})
ax.set_ylim(0, 2)
ax.set_xlim(0, 2)
fig.set_figwidth(12)
fig.set_figheight(9)

fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/exc/ExcErg"+str(numOfRuns)+".pdf")
fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/exc/ExcErg"+str(numOfRuns)+".png")
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
