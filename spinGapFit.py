from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
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


def exp_fit_fn(beta, A, B):
    return A * np.exp(-beta * B)


path = "D:/Code/C++/spinChainData/out/Susceptibilities_DQT/SuscDQTN12J0_500000It10.txt"
file = open(path)
lines = file.readlines()
betas = []
Cs = []
for line in lines:
    data = line.split(" ")
    if float(data[0]) < 25:
        continue
    betas.append(float(data[0]))
    Cs.append(float(data[1].replace("\n", "")))

parameters, covariance = opt.curve_fit(exp_fit_fn, betas, Cs)
print(parameters)
