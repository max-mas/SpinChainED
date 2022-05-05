from matplotlib import pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=20)

J_ratios = ["0_100000", "1_000000", "2_000000", "5_000000"]
Ts = ["0_200000", "0_500000", "1_000000"]

for T in Ts:
    fig, ax = plt.subplots()
    
    for N in np.linspace(6, 12, 4):
        path = "/home/mmaschke/BA_Code/Data/SpecificHeatsForJ/SpecHeatN"+str(int(N))+"T"+T+".txt"
        file = open(path)
        Js = []
        specHeat = []
        for i in range(200) :
            stri = file.readline()
            data = stri.split(" ")
            Js.append( float( data[0]) )
            specHeat.append( float( data[1].replace("\n", "") ) )
        lab = "$N$ = "+ str( int(N) )
        ax.plot(Js, specHeat, label=lab)
    ax.legend()
    T = T.replace("_", ".")
    ax.set(xlabel="$J_1/J_2$", ylabel="Specific heat per Spin $C/N$", title="$T =\\,$"+T)
    plt.show()

"""
for J_ratio in J_ratios:
    fig, ax = plt.subplots()
    
    for N in np.linspace(6, 12, 4):
        path = "/home/mmaschke/BA_Code/Data/SpecificHeats/SpecHeatN"+str(int(N))+"J"+J_ratio+".txt"
        file = open(path)
        Ts = []
        specHeat = []
        for i in range(200) :
            stri = file.readline()
            data = stri.split(" ")
            Ts.append( float( data[0]) )
            specHeat.append( float( data[1].replace("\n", "") ) )
        lab = "$N$ = "+ str( int(N) )
        ax.plot(Ts, specHeat, label=lab)
    ax.legend()
    J_ratio = J_ratio.replace("_", ".")
    ax.set(xlabel="$T$ ($J_2$)", ylabel="Specific heat per Spin $C/N$", title="$J_1/J_2 =\\,$"+J_ratio)
    plt.show()


for J_ratio in J_ratios:
    fig, ax = plt.subplots()
    
    for N in np.linspace(6, 12, 4):
        path = "/home/mmaschke/BA_Code/Data/Susceptibilities/SuscN"+str(int(N))+"J"+J_ratio+".txt"
        file = open(path)
        Ts = []
        susc = []
        for i in range(200) :
            stri = file.readline()
            data = stri.split(" ")
            Ts.append( float( data[0]) )
            susc.append( float( data[1].replace("\n", "") ) )
        lab = "$N$ = "+ str( int(N) )
        ax.plot(Ts, susc, label=lab)
    ax.legend()
    J_ratio = J_ratio.replace("_", ".")
    ax.set(xlabel="$T$ ($J_2$)", ylabel="Susceptibility per Spin $C$ ($1/J_2$)", title="$J_1/J_2 =\\,$"+J_ratio)
    plt.show()

fig, ax = plt.subplots()
for N in np.linspace(6, 12, 4):
    path = "/home/mmaschke/BA_Code/Data/ExcitationErgs/ExcErgs"+str(int(N))+".txt"
    file = open(path)
    Ts = []
    excErg = []
    for i in range(200) :
        stri = file.readline()
        data = stri.split(" ")
        Ts.append( float( data[0]) )
        excErg.append( float( data[1].replace("\n", "") ) )
    lab = "$N$ = "+ str( int(N) )
    ax.plot(Ts, excErg, label=lab)
    ax.legend()
    ax.set(xlabel="$J_1/J_2$", ylabel="1st Excitation Energy $\\Delta E$ ($J_2$)")
plt.show()
"""

