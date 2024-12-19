import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc
import linecache
from read_file import *

## eps形式で保存して、texに貼り付けたときに画像がずれたり切れたりしないようにするおまじない。
plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.rcParams['font.family'] = "Arial" #'Times New Roman'
# plt.rcParams["mathtext.fontset"] = "stixsans" #"stix" 
plt.rcParams['text.latex.preamble'] = r'\usepackage{sfmath}'
#[r'\usepackage[notext]{stix2}', r"\renewcommand{\rmdefault}{phv}", r"\renewcommand{\sfdefault}{phv}"]



dirname = "../output/"
model   = "St1em2"

outdir = "../tex/"

savefigname = outdir + model + "_mdot.eps"

# read input file
inputfile = dirname + model + "/input"
input_params = read_input(inputfile)
nr = int(input_params['nr'])

# logfile 
logfile = dirname + model + "/log"
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

mask = (mdot_inf > 0.0)


time1e3yr = time[mask] / (1e3*yr)
mdot_accMs_yr = mdot_acc[mask] * yr / M_sun
mdot_infMS_yr = mdot_inf[mask] * yr / M_sun

fig = plt.figure(figsize=(14, 22))
# plt.suptitle(r"$\mathrm{St} = 10^{-2}$", fontsize=30)

ax1 = fig.add_subplot(2, 1, 1)
# ax1.set_title(r"$\mathrm{St} = 10^{-3}$", fontsize=30)
ax1.plot(time1e3yr, mdot_accMs_yr, label=r"$\dot{M}_{\mathrm{acc}}$")
ax1.plot(time1e3yr, mdot_infMS_yr, label=r"$\dot{M}_{\mathrm{inf}}$")
# ax1.set_ylim(0.0, 4.0)
ax1.set_xlim(time1e3yr[0], time1e3yr[-1])
ax1.set_ylim(1.0e-8, 1.0e-4)
ax1.set_yscale("log")
ax1.set_xlabel(r"time $[10^{3} \ \mathrm{yr}]$", fontsize=30)
ax1.set_ylabel(r"$\dot{M} \ [10^{-5} \mathrm{M_{\odot} / yr}]$", fontsize=30)
ax1.tick_params(labelsize=30)
ax1.legend(fontsize=30, frameon=False)

mstarMsun = mstar[mask] / M_sun
mdiskMsun = (gas_total_mass[mask] + dust_total_mass[mask])/M_sun

ax2 = fig.add_subplot(2, 1, 2)
# ax2.set_title(r"$\mathrm{St} = 10^{-2}$", fontsize=30)
ax2.plot(time1e3yr, mstarMsun, label=r"$M_{\mathrm{star}}$")
ax2.plot(time1e3yr, mdiskMsun, label=r"$M_{\mathrm{disk}}$")
ax2.set_xlabel(r"time $[10^{3} \ \mathrm{yr}]$", fontsize=30)
ax2.set_ylabel(r"$M \ [M_{\odot}]$", fontsize=30)
ax2.set_xlim(time1e3yr[0], time1e3yr[-1])
ax2.tick_params(labelsize=30)
ax2.legend(fontsize=30, frameon=False)


# ax3 = fig.add_subplot(2, 2, 3)
# # ax3.set_title(r"$\mathrm{St} = 10^{-1}$", fontsize=30)
# ax3.plot(time/1000/yr, Ls/L_sun, label=r"$L_{s}$")
# ax3.plot(time/1000/yr, Lacc/L_sun, label=r"$L_{\mathrm{acc}}$")
# #ax3.set_ylim(0.0, 4.0)
# ax3.set_xlabel("time [1000yr]", fontsize=30)
# ax3.set_ylabel(r"$L \ [L_{\odot}]$", fontsize=30)
# ax3.tick_params(labelsize=30)
# ax3.legend(fontsize=30, frameon=False)


# ax4 = fig.add_subplot(2, 2, 4)
# ax4.set_title(r"$\mathrm{St} = 1$", fontsize=30)
# ax4.plot(time/1000/yr, mdot_acc*yr/M_sun/1e-5, label=r"$\dot{M}_{\mathrm{acc}}$")
# ax4.plot(time/1000/yr, mdot_inf*yr/M_sun/1e-5, label=r"$\dot{M}_{\mathrm{inf}}$")
# ax4.set_ylim(0.0, 4.0)
# ax4.set_xlabel("time [1000yr]", fontsize=30)
# ax4.set_ylabel(r"$\dot{M} \ [10^{-5} \mathrm{M_{\odot} / yr}]$", fontsize=30)
# ax4.tick_params(labelsize=30)
# ax4.legend(fontsize=30, frameon=False)

# plt.show()
fig.tight_layout()
fig.savefig(savefigname)
