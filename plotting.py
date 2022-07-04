from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import sys
import copy

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=20)
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


def weird_transform(Js, Vals):
    A = copy.deepcopy(Vals)

    for y in range(len(Js)):
        A[y] *= 1 / (1 + Js[y])
    return A


if len(sys.argv) != 8:
    print("Incorrect amount of arguments.")
    sys.exit(69)

nMin = int(sys.argv[1])
nMax = int(sys.argv[2])
nNum = int((nMax - nMin) / 2) + 1
datapointNum = int(sys.argv[3])
J_ratios = read_vec_from_file(sys.argv[4])
T_values = read_vec_from_file(sys.argv[5])
flags_str = sys.argv[6]
saveToPath = sys.argv[7]
flags = []

reps = [1, 2, 3]

for i in range(0, len(flags_str)):
    flag = int(flags_str[i])
    flags.append(bool(flag))

# Excitation Erg
if flags[0]:
    fig, ax = plt.subplots()
    for N in np.linspace(nMin, nMax, nNum):
        path = saveToPath + "/out/ActualExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
        file = open(path)
        lines = file.readlines()
        Js = []
        excErg = []
        for line in lines:
            data = line.split(" ")
            Js.append(float(data[0]))
            excErg.append(float(data[1].replace("\n", "")))
        lab = "$N$ = " + str(int(N))
        ax.plot(Js, weird_transform(Js, excErg), label=lab)
        ax.legend()
        ax.set(xlabel="$J_1/J_2$", ylabel="Reduced Excitation Energy $\\Delta/(J_1+J_2)$ ($J_2$)")
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 0.75)
    #fig.savefig(saveToPath + "/plots/ActualExcitationErgs/Excergs.pdf")
    #fig.savefig(saveToPath + "/plots/ActualExcitationErgs/Excergs.png")
    plt.show()
    plt.close()

# Ground state erg
if flags[1]:
    fig, ax = plt.subplots()
    for N in np.linspace(nMin, nMax, nNum):
        path = saveToPath + "/out/GroundStateErgs/GSErgs" + str(int(N)) + ".txt"
        file = open(path)
        lines = file.readlines()
        Ts = []
        excErg = []
        for line in lines:
            data = line.split(" ")
            Ts.append(float(data[0]))
            excErg.append(float(data[1].replace("\n", "")))
        lab = "$N$ = " + str(int(N))
        ax.plot(Ts, excErg, label=lab)
        ax.legend()
        ax.set(xlabel="$J_1/J_2$", ylabel="Ground state energy per spin $E_0 / N$ ($J_2$)")
    fig.savefig(saveToPath + "/plots/GroundStateErgs/GSErgs.pdf")
    fig.savefig(saveToPath + "/plots/GroundStateErgs/GSErgs.png")

# Specific Heat (T)
if flags[2]:
    for J_ratio in J_ratios:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            if np.abs(float(J_ratio.replace("_", ".")) - 0.5) < 1e-5 and N == 18:
                continue
            path = saveToPath + "/out/SpecificHeats/SpecHeatN" + str(int(N)) + "J" + J_ratio + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            specHeat = []
            for line in lines:
                data = line.split(" ")
                if float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                if N == 18:
                    Ts[-1] = 1/Ts[-1]
                specHeat.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, specHeat, label=lab)
        b = 1/np.asarray(Ts)
        ax.plot(Ts, 3*np.exp(b)/(np.exp(b)+3)**2 * b**2, "r-", label="High-T")
        ax.plot(Ts, 3/13/np.asarray(Ts)**2, "b-", label="High-T 2nd Order")
        ax.plot(Ts, 3/16/np.asarray(Ts)**2, "g-", label="High-T 2nd Order")
        ax.legend()
        ax.set_ylim(0, 0.4)
        ax.set_xlim(0, 4)
        J_ratio = J_ratio.replace("_", ".")
        ax.set(xlabel="$T$ ($J_2$)", ylabel="Specific heat per Spin $C/N$", title="$J_1/J_2=\\,$ " + J_ratio)
        J_ratio = J_ratio.replace(".", "_")
        #fig.savefig(saveToPath + "/plots/SpecificHeats/SpecHeatJ" + J_ratio + ".pdf")
        #fig.savefig(saveToPath + "/plots/SpecificHeats/SpecHeatJ" + J_ratio + ".png")
        plt.show()

# Specific Heat (J)
if flags[3]:
    for T in T_values:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/SpecificHeatsForJ/SpecHeatN" + str(int(N)) + "T" + T + ".txt"
            file = open(path)
            lines = file.readlines()
            Js = []
            specHeat = []
            for line in lines:
                data = line.split(" ")
                Js.append(float(data[0]))
                specHeat.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Js, specHeat, label=lab)
        ax.legend()
        T = T.replace("_", ".")
        ax.set(xlabel="$J_1/J_2$", ylabel="Specific heat per Spin $C/N$", title="$T =\\,$" + T)
        T = T.replace(".", "_")
        fig.savefig(saveToPath + "/plots/SpecificHeatsForJ/SpecHeatT" + T + ".pdf")
        fig.savefig(saveToPath + "/plots/SpecificHeatsForJ/SpecHeatT" + T + ".png")

# Susceptibility (T)
if flags[4]:
    for J_ratio in J_ratios:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/Susceptibilities/SuscN" + str(int(N)) + "J" + J_ratio + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            susc = []
            for line in lines:
                data = line.split(" ")
                if data[0] == "-nan" or data[0] == "nan" or float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                if N == 16:
                    Ts[-1] = 1/Ts[-1]
                susc.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, susc, label=lab)
        ax.legend()
        ax.set_xlim(0, 3)
        J_ratio = J_ratio.replace("_", ".")
        ax.set(xlabel="$T$ ($J_2$)", ylabel="Susceptibility per Spin $\\chi / N$", title="$J_1/J_2 =\\,$" + J_ratio)
        J_ratio = J_ratio.replace(".", "_")
        fig.savefig(saveToPath + "/plots/Susceptibilities/SuscJ" + J_ratio + ".pdf")
        fig.savefig(saveToPath + "/plots/Susceptibilities/SuscJ" + J_ratio + ".png")

# Susceptibility (J)
if flags[5]:
    for T in T_values:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/SusceptibilitiesForJ/SuscN" + str(int(N)) + "T" + T + ".txt"
            file = open(path)
            lines = file.readlines()
            Js = []
            chi = []
            for line in lines:
                data = line.split(" ")
                Js.append(float(data[0]))
                chi.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Js, chi, label=lab)
        ax.legend()
        T = T.replace("_", ".")
        ax.set(xlabel="$J_1/J_2$", ylabel="Susceptibility per Spin $\\chi / N$ ", title="$T =\\,$" + T)
        T = T.replace(".", "_")
        fig.savefig(saveToPath + "/plots/SusceptibilitiesForJ/SuscT" + T + ".pdf")
        fig.savefig(saveToPath + "/plots/SusceptibilitiesForJ/SuscT" + T + ".png")

# Dispersion
if flags[6]:
    for J_ratio in J_ratios:
        for N in np.linspace(nMin, nMax, nNum):
            fig, ax = plt.subplots()
            path = saveToPath + "/out/Dispersion/DispN" + str(int(N)) + "J" + J_ratio + ".txt"
            file = open(path)
            lines = file.readlines()
            ms = []
            ks = []
            ergs = []
            for line in lines:
                data = line.split(" ")
                ms.append(float(data[0]))
                ks.append(float(data[1]))
                ergs.append(float(data[2].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            sc = ax.scatter(ks, ergs, s=400, c=ms, cmap="hsv", marker="_", linewidths=0.5)
            ax.legend()
            J_ratioNum = J_ratio.replace("_", ".")
            ax.set(xlabel="Momentum quantum number $k$", ylabel="State Energy $E$ ($J_2$)",
                   title="$J_1/J_2 =\\,$" + J_ratioNum + ", $N= $ " + str(N))
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            fig.colorbar(sc, cmap="hsv", values=ms, label="Magnetization $m$")
            fig.savefig(saveToPath + "/plots/Dispersion/DispN" + str(N) + "J" + J_ratio + ".pdf")
            fig.savefig(saveToPath + "/plots/Dispersion/DispN" + str(N) + "J" + J_ratio + ".png")

"""
# QT Error stats for C (WIP)
reps = [1, 3, 10]
for J_ratio in J_ratios:
    for rep in reps:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/QTErrorStats/SpecHeatDiffs/DiffN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            deltas = []
            for line in lines:
                data = line.split(" ")
                if data[0] == "-nan" or data[0] == "nan" or float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                deltas.append(float(data[1]))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, deltas, label=lab)
        dBeta = 1/Ts[0]
        ax.legend()
        J_ratioNum = J_ratio.replace("_", ".")
        ax.set(xlabel="Temperature $T$ ($J_2$)", ylabel="Relative Error of DQT spec. Heat $\\delta C$",
                   title="$J_1/J_2 =\\,$" + J_ratioNum + ", $n =$ " + str(rep) + ", $d\\beta$ =" + str(dBeta))
        ax.set_xlim(0, 10)
        ax.semilogy()
        fig.savefig(saveToPath + "/plots/QTErrorStats/SpecHeatDiffs/Diff" + "J" + J_ratio + "It" + str(rep) + ".pdf")
        fig.savefig(saveToPath + "/plots/QTErrorStats/SpecHeatDiffs/Diff" + "J" + J_ratio + "It" + str(rep) + ".png")
        plt.close()

# QT Error stats for X (WIP)
for J_ratio in J_ratios:
    for rep in reps:
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/QTErrorStats/SuscDiffs/DiffN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            deltas = []
            for line in lines:
                data = line.split(" ")
                if data[0] == "-nan" or data[0] == "nan" or float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                deltas.append(float(data[1]))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, deltas, label=lab)
        dBeta = 1/Ts[0]
        ax.legend()
        J_ratioNum = J_ratio.replace("_", ".")
        ax.set(xlabel="Temperature $T$ ($J_2$)", ylabel="Realative Error of DQT Susceptibility $\\delta\\chi$",
               title="$J_1/J_2 =\\,$" + J_ratioNum + ", $n =$ " + str(rep) + ", $d\\beta$ =" + str(dBeta))
        ax.set_xlim(0, 10)
        ax.semilogy()
        fig.savefig(saveToPath + "/plots/QTErrorStats/SuscDiffs/Diff" + "J" + J_ratio + "It" + str(rep) + ".pdf")
        fig.savefig(saveToPath + "/plots/QTErrorStats/SuscDiffs/Diff" + "J" + J_ratio + "It" + str(rep) + ".png")
        plt.close()

maxErrs = []
# Specific Heat DQT (T)
for J_ratio in J_ratios:
    maxErrs.append([])
    for rep in reps:
        maxErrs[-1].append([])
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/SpecificHeats_DQT/SpecHeatDQTN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".txt"
            if N > 16:
                path = "/home/mmaschke/BA_Code/remoteData/out/SpecificHeats_DQT/SpecHeatDQTN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".txt"
            if N > 16 and rep == 10:
                continue
            file = open(path)
            lines = file.readlines()
            Ts = []
            specHeat = []
            errs = []
            for line in lines:
                data = line.split(" ")
                if float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                specHeat.append(float(data[1]))
                if N < 18:
                    errs.append(float(data[2].replace("\n", ""))/np.sqrt(rep))
                else:
                    errs.append(0)
            lab = "$N$ = " + str(int(N)) + ", QT"
            dBeta = 1/Ts[0]
            Ts = np.asarray(Ts)
            specHeat = np.asarray(specHeat)
            errs = np.asarray(errs)
            relErrs = []
            for u in range(len(errs)):
                if specHeat[u] > 0 and Ts[u] > 0.25:
                    relErrs.append(errs[u]/np.abs(specHeat[u]))
            ax.plot(Ts, specHeat, label=lab)
            ax.fill_between(Ts, specHeat - errs, specHeat + errs, alpha=0.1)
            maxErrs[-1][-1].append(np.mean(relErrs))

        #plt.gca().set_prop_cycle(None)
        #for N in np.linspace(nMin, nMax, nNum):
        #    path = saveToPath + "/out/SpecificHeats/SpecHeatN" + str(int(N)) + "J" + J_ratio + ".txt"
        #    file = open(path)
        #    lines = file.readlines()
        #    Ts = []
        #    specHeat = []
        #    for line in lines:
        #        data = line.split(" ")
        #        if float(data[0]) == 0:
        #            continue
        #        Ts.append(1/float(data[0]))
        #        specHeat.append(float(data[1].replace("\n", "")))
        #    lab = "$N$ = " + str(int(N)) + ", ED"
        #    ax.plot(Ts, specHeat, "--", label=lab, alpha = 0.4)
        #    dBeta = 1/Ts[0]

        ax.legend()
        #ax.set_ylim(0, 0.4)
        ax.set_xlim(0, 2)
        J_ratio = J_ratio.replace("_", ".")
        ax.set(xlabel="$\\beta$ ($1/J_2$)", ylabel="Specific heat per Spin $C/N$", title="$J_1/J_2 =\\,$" + J_ratio + ", $n=$ " + str(rep) + ", $d\\beta$ =" + str(dBeta))
        J_ratio = J_ratio.replace(".", "_")
        fig.savefig("/home/mmaschke/BA_Code/Data/plots/SpecificHeats_DQT/SpecHeatJ" + J_ratio + "It" + str(rep) + ".pdf")
        fig.savefig("/home/mmaschke/BA_Code/Data/plots/SpecificHeats_DQT/SpecHeatJ" + J_ratio + "It" + str(rep) + ".png")
        #plt.show()
        plt.close()

i = 0
N_list = np.linspace(nMin, nMax, nNum)
for J_ratio in J_ratios:
    j = 0
    fig, ax = plt.subplots()
    for rep in reps:
        lab = "$n=$ " + str(rep)
        ax.plot(N_list, maxErrs[i][j], ".-", label=lab)
        j += 1
    ax.legend()
    ax.legend(prop={'size': 6})
    J_ratio = J_ratio.replace("_", ".")
    ax.set_xlabel("$N$", fontsize = 24)
    ax.set_ylabel("Mittlerer relativer Standardfehler des Mittelwerts $\\overline\\sigma^{-}_{C / N}$", fontsize = 24)
    ax.set_title("$J_1/J_2 =$ " + J_ratio + ", $d\\beta =$ " + str(dBeta) + ", $T>0.25$", fontsize = 30)
    fig.set_figwidth(12)
    fig.set_figheight(9)
    J_ratio = J_ratio.replace(".", "_")
    ax.semilogy()
    fig.savefig(saveToPath + "/plots/QTErrorStats/SpecHeatMeanSEOM/SpecHeatMeanSEOM" + J_ratio + "rel" + ".pdf")
    fig.savefig(saveToPath + "/plots/QTErrorStats/SpecHeatMeanSEOM/SpecHeatMeanSEOM" + J_ratio + "rel" + ".png")
    #plt.show()
    plt.close()
    i += 1

# Susceptibility DQT (T)
maxErrs = []
for J_ratio in J_ratios:
    maxErrs.append([])
    for rep in reps:
        maxErrs[-1].append([])
        fig, ax = plt.subplots()
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/Susceptibilities_DQT/SuscDQTN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            susc = []
            errs = []
            for line in lines:
                data = line.split(" ")
                if data[0] == "-nan" or data[0] == "nan" or float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                susc.append(float(data[1]))
                errs.append(float(data[2].replace("\n", ""))/np.sqrt(rep))
            Ts = np.asarray(Ts)
            susc = np.asarray(susc)
            errs = np.asarray(errs)
            relErrs = []
            for u in range(len(errs)):
                if susc[u] > 0 and Ts[u] > 0.25:
                    relErrs.append(errs[u]/np.abs(susc[u]))
            lab = "$N$ = " + str(int(N)) + ", QT"
            ax.plot(Ts, susc, label=lab)
            ax.fill_between(Ts, susc - errs, susc + errs, alpha=0.1)
            dBeta = 1/Ts[0]
            maxErrs[-1][-1].append(np.mean(relErrs))

        plt.gca().set_prop_cycle(None)
        for N in np.linspace(nMin, nMax, nNum):
            path = saveToPath + "/out/Susceptibilities/SuscN" + str(int(N)) + "J" + J_ratio + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            susc = []
            for line in lines:
                data = line.split(" ")
                if data[0] == "-nan" or data[0] == "nan" or float(data[0]) == 0:
                    continue
                Ts.append(1/float(data[0]))
                susc.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N)) + ", ED"
            ax.plot(Ts, susc, "--", label=lab, alpha=0.4)

        ax.legend()
        ax.set_xlim(0, 1)
        J_ratio = J_ratio.replace("_", ".")
        ax.set(xlabel="$T$ ($J_2$)", ylabel="Susceptibility per Spin $\\chi / N$", title="$J_1/J_2 =\\,$" + J_ratio + ", $n=$ " + str(rep) + ", $d\\beta$ =" + str(dBeta))
        J_ratio = J_ratio.replace(".", "_")
        fig.savefig(saveToPath + "/plots/QTErrorStats/SuscQTandED/SuscN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".pdf")
        fig.savefig(saveToPath + "/plots/QTErrorStats/SuscQTandED/SuscN" + str(int(N)) + "J" + J_ratio + "It" + str(rep) + ".png")
        #plt.show()
        plt.close()

i = 0
N_list = np.linspace(nMin, nMax, nNum)
for J_ratio in J_ratios:
    j = 0
    fig, ax = plt.subplots()
    for rep in reps:
        lab = "$n=$ " + str(rep)
        ax.plot(N_list, maxErrs[i][j], ".-", label=lab)
        j += 1
    ax.legend()
    ax.legend(prop={'size': 6})
    J_ratio = J_ratio.replace("_", ".")
    ax.set_xlabel("$N$", fontsize = 24)
    ax.set_ylabel("Mittlerer relativer Standardfehler des Mittelwerts $\\overline\\sigma^{-}_{\\chi / N}$", fontsize = 24)
    ax.set_title("$J_1/J_2 =$ " + J_ratio + ", $d\\beta =$ " + str(dBeta), fontsize = 30)
    fig.set_figwidth(12)
    fig.set_figheight(9)
    J_ratio = J_ratio.replace(".", "_")
    ax.semilogy()
    fig.savefig(saveToPath + "/plots/QTErrorStats/SuscMeanSEOM/SuscMeanSEOM" + J_ratio + ".pdf")
    fig.savefig(saveToPath + "/plots/QTErrorStats/SuscMeanSEOM/SuscMeanSEOM" + J_ratio + ".png")
    #plt.show()
    plt.close()
    i += 1
"""