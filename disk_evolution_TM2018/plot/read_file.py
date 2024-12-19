import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc
import linecache

## 物理定数
pi           = np.pi                    # PI
k_b          = sc.k*1e7                 # Boltzmann constant in erg/K
m_p          = sc.proton_mass*1e3       # proton mass in g
Grav         = sc.G*1e3                 # gravitational constant in cm^3 g^-1 s^-2
AU           = sc.au*1e2                # astronomical unit in cm
yr           = sc.Julian_year             # year in s
mu           = 2.34e0                    # mean molecular mass in proton masses
M_sun        = 1.9891e+33               # mass of the sun in g
R_sun        = 69550800000.0            # radius of the sun in cm
L_sun        = 3.828000000000000255e+33
sigmaSB      = 5.670500000000000003e-05

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
ID_MAX    = 23  # number of data ID

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
        data = data[12:].reshape([ID_MAX, nr])
    return time, data