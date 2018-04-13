# -*- coding: utf-8 -*-
"""
Ising main
"""
import numpy as np
from labutil.src.plugins.ising import *
import matplotlib.pyplot as plt
# Variable initialization

T = 2 # Temperature range to scan (kb = 1)
N = 12                      # Spin lattice size

n_eq = 75          # Average number of equillibriation steps (flips) per site
n_mc = 100          # Average number of Monte Carlo steps (flips) per site

#### RUN MC CODE OVER ####
import time
start_time = time.time()
E, M, E_eq, M_eq = mc_run(N,n_eq,n_mc,T)
print("--- %s seconds ---" % (time.time() - start_time))

x_eq = np.arange(n_eq*N*N)
x_mc = np.arange(n_mc*N*N)+(n_eq*N*N)

fig, (ax_energy, ax_mag) = plt.subplots(2,1, sharex=True, sharey=False)

ax_energy.plot(x_eq,E_eq,label='Equilibration')
ax_energy.plot(x_mc,E,label='Production')
ax_energy.axvline(n_eq*(N**2),color = 'black')
ax_energy.legend()
ax_energy.set_ylabel('Energy per site/ J')

ax_mag.plot(x_eq,M_eq,label='Equilibration')
ax_mag.plot(x_mc,M,label='Production')
ax_mag.axvline(n_eq*(N**2),color = 'black')
ax_mag.legend()
ax_mag.set_ylabel('Magnetization per site')
ax_mag.set_xlabel('Number of flip attempts')


plt.show()