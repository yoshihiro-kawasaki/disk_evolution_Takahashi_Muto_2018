"""
Takahashi et al 2013, 2018
tinf = (2/pi)*tff*(factor by density enhancement)
This script calculates the factor by density enhancement
"""

import numpy as np
from scipy import integrate

f = 1.4

def func(x):
    return 1.0 / np.sqrt(1.0/f*np.log(x) + 1.0/x - 1.0)

ans, err = integrate.quad(func, 0.0, 1.0)

print(ans)