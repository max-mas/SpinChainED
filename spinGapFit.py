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
    return beta * A * np.exp(-beta * B)


def alt_fit_fn(beta, A, B):
    return A + 2 * np.log(beta) - B * beta


def fd_up_fit_fn(beta, A, C, B):
    return beta**2 * A * ( B**2 * (np.exp(-2*beta*B) + np.exp(-beta*(B+C))) + C**2 * (np.exp(-2*beta*C) + np.exp(-beta*(B+C))) )


def weird_transform(Js, Vals):
    A = copy.deepcopy(Vals)

    for y in range(len(Js)):
        A[y] *= 1 / (1 + Js[y])
    return A


nMin = 20
nMax = 20
nNum = int((nMax - nMin) / 2) + 1
numOfRuns = 1
resids = True
diffs = False

gapsQT = []
gapsQTavg = []
gapsQTdev = []

path = "/home/mmaschke/BA_Code/remoteData/out/Susceptibilities_DQT/forFit/"
i = 0
for N in np.linspace(nMin, nMax, nNum):
    gapsQT.append([])
    j = 0
    for k in range(1, numOfRuns + 1):
        if N == 18 and k > 1:
            continue

        runPath = path + str(k) + "/"

        gapsQT[i].append([])
        gapsQT[i][k - 1].append([])
        gapsQT[i][k - 1].append([])
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
            fullBetas = []
            Cs = []
            fullCs = []
            lines = file.readlines()
            cutoff = 20
            upperCutoff = 50
            #if J > 1.5:
            #    cutoff = 15
            #elif 0.4 < J < 0.75:
            #    upperCutoff = 30

            for line in lines:
                data = line.split(" ")
                fullBetas.append(float(data[0]))
                if data[1] == "nan\n" or data[1] == "-nan\n" or data[1] == "inf\n" or data[1] == "-inf\n":
                    data[1] = "0"
                fullCs.append(float(data[1].replace("\n", "")))
                if float(data[0]) < cutoff or float(data[0]) > upperCutoff:
                    continue
                betas.append(float(data[0]))
                Cs.append(float(data[1].replace("\n", "")))

            try:
                #if 0.4 < J < 0.75:
                #    for m in range(len(Cs)):
                #        if Cs[m] < 1e-12:
                #            Cs[m] = 1e-12
                #    lnCs = np.log(Cs)
                #    parameters, covariance = opt.curve_fit(alt_fit_fn, betas, lnCs)
                #else:
                    parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
            except RuntimeError:
                parameters = [0, 0]
                covariance = [0, 0]
                print("Fit Error for N = " + str(N) + " and J = " + str(J))
            perr = np.sqrt(np.diag(covariance))
            gapsQT[i][k - 1][0].append(J)
            gapsQT[i][k - 1][1].append(parameters[1])
            if k == 1 and resids:
                fullFitted = []
                fitted = []
                for beta in fullBetas:
                    try:
                        fullFitted.append(exp_fit_fn(beta, parameters[0], parameters[1]))
                        if beta > cutoff:
                            fitted.append(exp_fit_fn(beta, parameters[0], parameters[1]))
                    except:
                        fullFitted.append(0)
                        if beta > cutoff:
                            fitted.append(0)
                fig2, ax2 = plt.subplots()
                ax2.plot(fullBetas, fullCs, label="QT", alpha=0.4)
                ax2.plot(fullBetas, fullFitted, label="Fit")
                ax2.legend()
                ax2.semilogy()
                ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility per Spin $\\chi/N$", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
                fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/QT/fitN"+str(N)+"J"+J_str+".png")
                plt.close(fig2)
                fig2, ax2 = plt.subplots()
                ax2.plot(fullBetas, np.abs(np.asarray(fullFitted)-np.asarray(fullCs)))
                ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility fit residual", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
                ax2.semilogy()
                fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/QTResid/residN"+str(N)+"J"+J_str+".png")
                plt.close(fig2)
                fig2, ax2 = plt.subplots()
                ax2.plot(fullBetas, np.abs((np.asarray(fullFitted)-np.asarray(fullCs))/np.asarray(fullCs)))
                ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility fit relative residual", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
                ax2.semilogy()
                fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/QTResidRel/relresidN"+str(N)+"J"+J_str+".png")
                plt.close(fig2)
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
            if N == 18 and l == 0:
                avg += gapsQT[i][l][1][n]
                break
            avg += gapsQT[i][l][1][n] / numOfRuns
        for l in range(numOfRuns):
            if N == 18 and l == 0:
                dev += 0
                break
            dev += (gapsQT[i][l][1][n] - avg) ** 2 / numOfRuns
        dev = np.sqrt(dev)/np.sqrt(numOfRuns)
        gapsQTavg[i][1].append(avg)
        gapsQTdev[i][1].append(dev)

    i += 1

fig, ax = plt.subplots()
N = nMin
i = 0
for arr in gapsQTavg:
    lab = "$N$ = " + str(int(N)) + " QT Fit"
    ax.errorbar(arr[0], weird_transform(arr[0], arr[1]), yerr=weird_transform(arr[0], gapsQTdev[i][1]), xerr=None,
            fmt=".-", label=lab, alpha=0.5, capsize=2)
    N += 2
    i += 1
"""
gapsED = []
path = "/home/mmaschke/BA_Code/remoteData/out/Susceptibilities/forFit/"
i = 0
for N in np.linspace(6, 16, 6):
    gapsED.append([])
    gapsED[i].append([])
    gapsED[i].append([])
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
        fullBetas = []
        Cs = []
        fullCs = []
        lines = file.readlines()
        cutoff = 20
        upperCutoff = 50
        #if J > 1.5:
        #    cutoff = 15
        #elif 0.4 < J < 0.75:
        #    upperCutoff = 30

        for line in lines:
            data = line.split(" ")
            fullBetas.append(float(data[0]))
            if data[1] == "nan\n" or data[1] == "-nan\n" or data[1] == "inf\n" or data[1] == "-inf\n":
                data[1] = "0"
            fullCs.append(float(data[1].replace("\n", "")))
            if float(data[0]) < cutoff or float(data[0]) > upperCutoff:
                continue
            betas.append(float(data[0]))
            Cs.append(float(data[1].replace("\n", "")))

        try:
            #if 0.4 < J < 0.75:
            #    for m in range(len(Cs)):
            #        if Cs[m] < 1e-12:
            #            Cs[m] = 1e-12
            #    lnCs = np.log(Cs)
            #    parameters, covariance = opt.curve_fit(alt_fit_fn, betas, lnCs)
            #else:
                parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
        except RuntimeError:
            parameters = [0, 0]
            covariance = [0, 0]
            print("Fit Error for N = " + str(N) + " and J = " + str(J))
        perr = np.sqrt(np.diag(covariance))
        gapsED[i][0].append(J)
        gapsED[i][1].append(parameters[1])
        if resids:
            fullFitted = []
            fitted = []
            for beta in fullBetas:
                try:
                    fullFitted.append(exp_fit_fn(beta, parameters[0], parameters[1]))
                    if beta > cutoff:
                        fitted.append(exp_fit_fn(beta, parameters[0], parameters[1]))
                except:
                    fullFitted.append(0)
                    if beta > cutoff:
                        fitted.append(0)
            fig2, ax2 = plt.subplots()
            ax2.plot(fullBetas, fullCs, label="ED", alpha=0.4)
            ax2.plot(fullBetas, fullFitted, label="Fit")
            ax2.legend()
            ax2.semilogy()
            ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility per Spin $\\chi/N$", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
            fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/ED/fitN"+str(N)+"J"+J_str+".png")
            plt.close(fig2)
            fig2, ax2 = plt.subplots()
            ax2.plot(fullBetas, np.abs(np.asarray(fullFitted)-np.asarray(fullCs)))
            ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility fit residual", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
            ax2.semilogy()
            fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/EDResid/residN"+str(N)+"J"+J_str+".png")
            plt.close(fig2)
            fig2, ax2 = plt.subplots()
            ax2.plot(fullBetas, np.abs((np.asarray(fullFitted)-np.asarray(fullCs))/np.asarray(fullCs)))
            ax2.set(xlabel="$\\beta$ $1/J_2$", ylabel="Susceptibility fit relative residual", title="$J=$ "+ J_str + ", $N=$ " + str(int(N)))
            ax2.semilogy()
            fig2.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/EDResidRel/relresidN"+str(N)+"J"+J_str+".png")
            plt.close(fig2)
        print(str(N) + " " + J_str + " " + str(perr[1]))
    i += 1
path = "/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/GapN"
N = nMin
plt.gca().set_prop_cycle(None)
for arr in gapsQTavg:
    lab = "$N$ = " + str(int(N)) + " ED Fit"
    ax.plot(arr[0], weird_transform(arr[0], arr[1]), ".-", label=lab, alpha=0.2)
    filePath = path + str(N) + "It" + str(numOfRuns)
    file = open(filePath, "w")
    for i in range(len(arr[0])):
        file.write(str(arr[0][i]) + " " + str(arr[1][i]) + "\n")
    file.close()
    N += 2
"""
if diffs:
    meanRelDiffsN = [[], []]
    N = nMin
    path = "/home/mmaschke/BA_Code/Data/plots/GapFit/spin/Fits/MeanDiffs/MeanRelDiffsIt" + str(numOfRuns) + ".txt"
    file = open(path, "w")
    plt.gca().set_prop_cycle(None)
    for i in range(len(gapsQTavg)):
        diffs = np.abs((np.asarray(gapsQTavg[i][1]) - np.asarray(gapsED[i][1]))/np.asarray(gapsED[i][1]))
        meanRelDiffsN[0].append(N)
        meanRelDiffsN[1].append(np.mean(diffs))
        file.write(str(N) + " " + str(meanRelDiffsN[1][-1]) + "\n")
        N += 2
    file.close()

plt.gca().set_prop_cycle(None)
for N in np.linspace(6, 18, 7):
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
    ax.plot(Ts, weird_transform(Ts, excErg), label=lab)

ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Spin Gap Energy $\\Delta/(J_1+J_2)$", title="$n = $ " + str(numOfRuns))
ax.legend(prop={'size': 6})
ax.set_ylim(0, 0.8)
ax.set_xlim(0, 2)

#fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/ExcErg" + str(numOfRuns) + ".pdf")
#fig.savefig("/home/mmaschke/BA_Code/Data/plots/GapFit/spin/ExcErg" + str(numOfRuns) + ".png")
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
