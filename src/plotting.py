from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import sys

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

for i in range(0, len(flags_str)):
    flag = int(flags_str[i])
    flags.append(bool(flag))

# Excitation Erg
if flags[0]:
    fig, ax = plt.subplots()
    for N in np.linspace(nMin, nMax, nNum):
        path = saveToPath + "/out/ExcitationErgs/ExcErgs" + str(int(N)) + ".txt"
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
        ax.set(xlabel="$J_1/J_2$", ylabel="Spin Gap Energy $\\Delta E$ ($J_2$)")
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 0.75)
    fig.savefig(saveToPath + "/plots/ExcitationErgs/Excergs.pdf")
    fig.savefig(saveToPath + "/plots/ExcitationErgs/Excergs.png")

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
            path = saveToPath + "/out/SpecificHeats/SpecHeatN" + str(int(N)) + "J" + J_ratio + ".txt"
            file = open(path)
            lines = file.readlines()
            Ts = []
            specHeat = []
            for line in lines:
                data = line.split(" ")
                Ts.append(float(data[0]))
                specHeat.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, specHeat, label=lab)
        ax.legend()
        J_ratio = J_ratio.replace("_", ".")
        ax.set(xlabel="$T$ ($J_2$)", ylabel="Specific heat per Spin $C/N$", title="$J_1/J_2 =\\,$" + J_ratio)
        J_ratio = J_ratio.replace(".", "_")
        fig.savefig(saveToPath + "/plots/SpecificHeats/SpecHeatJ" + J_ratio + ".pdf")
        fig.savefig(saveToPath + "/plots/SpecificHeats/SpecHeatJ" + J_ratio + ".png")

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
                if data[0] == "-nan" or data[0] == "nan":
                    continue
                Ts.append(float(data[0]) * 8.5)
                susc.append(float(data[1].replace("\n", "")))
            lab = "$N$ = " + str(int(N))
            ax.plot(Ts, susc, label=lab)
        ax.legend()
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
