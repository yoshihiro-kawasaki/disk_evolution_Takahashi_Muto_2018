# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from matplotlib import gridspec
import re
import os
import linecache
import scipy.constants as sc

# 物理定数
pi           = np.pi                    # PI
k_b          = sc.k*1e7                 # Boltzmann constant in erg/K
m_p          = sc.proton_mass*1e3       # proton mass in g
Grav         = sc.G*1e3                 # gravitational constant in cm^3 g^-1 s^-2
AU           = sc.au*1e2                # astronomical unit in cm
yr         = sc.Julian_year           # year in s
mu           = 2.3e0                    # mean molecular mass in proton masses
M_sun        = 1.9891e+33               # mass of the sun in g
R_sun        = 69550800000.0            # radius of the sun in cm
L_sun        = 3.828000000000000255e+33
sigmaSB      = 5.670500000000000003e-05

# dir
dir_name = "../output/St1em2_evap_only"
savename = "St1em2_evap_only_mdot.mp4"

plot_sigma = True # if True, plot Sigma, else plot T

## data ID
ID_RADIUS = 0   # radius
ID_SIGMAG = 1   # gas surface density
ID_TEMP   = 2   # temperature
ID_CS     = 3   # sound speed
ID_OMEGA  = 4   # angular velocity
ID_HG     = 5   # gas scale height
ID_QT     = 6   # Toomre's Q parameter
ID_MR     = 7   # enclosed mass
ID_ALPHA  = 8   # alpha parameter
ID_VGR    = 9   # gas radial velocity
ID_VGRVIS = 10  # viscous velocity
ID_VGRETA = 11  # gas pressure gradient velocity
ID_VGRSRC = 12  # velocty by source term
ID_CVGVIS = 13  # gas viscous velocity coefficient by gas and dust interaction
ID_CVGETA = 14  # gas pressure gradient velocity coefficient by gas and dust interaction
ID_QVIS   = 15  # viscous heating
ID_QIRR   = 16  # irradiation heating
ID_QINF   = 17  # infalling material heating
ID_QRAD   = 18  # radiation cooling
ID_SIGMAD = 19  # dust surface density
ID_VDR    = 20  # dust radial velocity
ID_CVDVIS = 21  # dust viscous velocity coefficient by gas and dust interaction
ID_CVDETA = 22  # gas pressure gradient velocity coefficient by gas and dust interaction
ID_FV     = 23
ID_MAX    = 24  # number of data ID

# read function
def read_input(filename):
    file = open(filename, 'r', encoding='UTF-8')
    d = file.readlines()
    nlen = len(d)
    input_prams = {}
    for i in range(nlen):
        line = d[i].split()
        if len(line) == 0:
            continue
        if (line[0] == '\n') or (line[0] == '#'):
            continue
        input_prams[line[0]] = line[2]
    return input_prams
        

def read_data(filename, nr):
    with open(filename) as f:
        data = np.fromfile(f, dtype=float)
        time = data[0]
        data = data[13:].reshape([ID_MAX, nr])
    return time, data

## input file
inputfile = dir_name + "/input"
input_params = read_input(inputfile)
nr = int(input_params['nr'])

## logfile
logfile = dir_name + "/log"
log_data = np.loadtxt(logfile, dtype=float).T
count = log_data[0]
time  = log_data[1]
mstar = log_data[2]
Ls = log_data[3]
Lacc = log_data[4]
gas_total_mass = log_data[5]
dust_total_mass = log_data[6]
total_mass = log_data[7]
mdot_acc = log_data[8]
mdot_inf = log_data[9]
total_infall_mass = log_data[10]
mdot_wind = log_data[11]
total_wind_loss_mass = log_data[12]

mdot_acc = mdot_acc * yr / M_sun
mdot_inf = mdot_inf * yr / M_sun
time1e3yr = time / (1e3*yr)

ndata = int(len(count)/2)
print(ndata)

filename = dir_name + "/disk" + str(int(count[0]))
t, data = read_data(filename, nr)
rAU = data[ID_RADIUS] / AU
del data

fig, axs = plt.subplots(1, 2, figsize=(26, 10), tight_layout=True)

if plot_sigma:
    ax01, = axs[0].plot([], [], lw=4, color="black", label="Gas")
    ax02, = axs[0].plot([], [], linestyle="dashed", lw=4, color="black", label="Dust")
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_xlim(0.1, 1e3)
    axs[0].set_ylim(1.0e-4, 1.0e7)
    axs[0].set_xlabel("r [au]", fontsize=30)
    axs[0].set_ylabel("Sigma [g cm^-2]", fontsize=30)
    axs[0].tick_params(labelsize=30)
    axs[0].legend(fontsize=30, frameon=False)
else:
    ax01, = axs[0].plot([], [], lw=4, color="black")
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_xlim(0.5, 1e3)
    axs[0].set_ylim(5.0, 1.0e4)
    axs[0].set_xlabel("r [au]", fontsize=30)
    axs[0].set_ylabel("T [K]", fontsize=30)
    axs[0].tick_params(labelsize=30)

# axs[1].plot(time1e3yr, mdot_acc, color="blue")
# ax11, = axs[1].plot(0, 0, marker="o", markersize=10, color="blue")
# axs[1].plot(time1e3yr, mdot_inf, color="orange")
# ax12, = axs[1].plot(0, 0, marker="o", markersize=10, color="orange")
# # axs[1].set_xscale("log")
# axs[1].set_yscale("log")
# axs[1].set_xlim(time1e3yr[0], time1e3yr[-1])
# axs[1].set_ylim(1.0e-7, 1.0e-4)
# axs[1].set_xlabel("time [1000 yr]", fontsize=30)
# axs[1].set_ylabel("Mdot", fontsize=30)
# axs[1].tick_params(labelsize=30)

axs[1].plot(time1e3yr, mdot_acc/1e-5, color="blue", label="Mdot_acc")
ax11, = axs[1].plot(0, 0, marker="o", markersize=10, color="blue")
axs[1].plot(time1e3yr, mdot_inf/1e-5, color="orange", label="Mdot_inf")
ax12, = axs[1].plot(0, 0, marker="o", markersize=10, color="orange")
# axs[1].set_xscale("log")
# axs[1].set_yscale("log")
axs[1].set_xlim(time1e3yr[0], time1e3yr[-1])
axs[1].set_ylim(0.0, 4.0)
axs[1].set_xlabel("time [1000 yr]", fontsize=30)
axs[1].set_ylabel("Mdot", fontsize=30)
axs[1].tick_params(labelsize=30)
axs[1].legend(fontsize=30, frameon=False)

def plot(idx):
    n = 2*idx
    filename = dir_name + "/disk" + str(int(count[n]))
    is_file = os.path.isfile(filename)
    if not is_file:
        return
    t, data = read_data(filename, nr)
    title = "{:.2e}".format(t/yr) + "year"
    # # 面密度
    if plot_sigma:
        ax01.set_data(rAU, data[ID_SIGMAG])
        ax02.set_data(rAU, data[ID_SIGMAD])
        axs[0].set_title(title, fontsize=30)
    else:
        #温度
        ax01.set_data(rAU, data[ID_TEMP])
        axs[0].set_title(title, fontsize=30)
        
    #mdot
    t1e3yr = t /(1e3 * yr)
    ax11.set_data(t1e3yr, mdot_acc[n]/1e-5)
    ax12.set_data(t1e3yr, mdot_inf[n]/1e-5)

ani = animation.FuncAnimation(fig, plot, interval=100, frames=ndata)
# plt.show()
# ani.save(savename, writer = "pillow")
ani.save(savename, writer = "ffmpeg")
# ani.save(savename, writer = "imagemagick")