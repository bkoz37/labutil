# -*- coding: utf-8 -*-
"""
Predefined functions for AP275 lab 4 -
Ising model on 2D square lattice
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def initiate_state(N):
    """
    Given an integer N, outputs a randomly generated square spin lattice
    of size N x N
    """
    N = int(N)
    state = 2*np.random.randint(2,size=(N,N)) - 1   # Each site randomly -1 or 1
    return state


def state_hamiltonian(state):
    """
    Given a spin state, sum through each nearest neighbor pair to find the
    energy of the spin hamiltonian.
    
    A factor of 1/2 is used to adjust for double counting
    """
    N = state.shape[0]
    energy = 0
    for ii in range(N):
        for jj in range(N):
            
            # Sum through sites
            site = state[ii, jj]
            
            #Sum nearest neighbor bonds, using modulo N to enforce periodic boundaries
            nn = state[(ii + 1) % N, jj] + state[(ii - 1) % N, jj] + \
                 state[ii, (jj + 1) % N] + state[ii, (jj - 1) % N]
                            
            # Update energy, eliminate double counting
            energy += -0.5*site*nn
    return energy


def state_magnetization(state):
    """
    Given a spin state, returns the total magnetization
    """
    return np.sum(state)


def spin_flip(state, beta):
    """
    Given a spin state, randomly select a site and flip it.
    
    If deltaE < 0, accept
    else if rand < exp(-deltE*beta), accept
    else reject
    """
    N = state.shape[0]
    
    # Randomly select site
    ii = np.random.randint(0, N)
    jj = np.random.randint(0, N)
    site_spin =  state[ii, jj]
    
    # Sum site's nearest neighbors
    nn = state[(ii + 1) % N, jj] + state[(ii - 1) % N, jj] + state[ii, (jj + 1) % N] + state[ii, (jj - 1) % N]
    
    deltaE = 2*site_spin*nn
    
    if np.random.rand() < np.exp(-deltaE*beta):
        site_spin *= -1
        state[ii, jj] = site_spin
        deltaM = 2 * site_spin
    else:
        deltaE = 0
        deltaM = 0

    return state, deltaE, deltaM

def animator(imagelist):
    fig_animation = plt.figure() # make figure
    im = plt.imshow(imagelist[0], cmap=plt.get_cmap('gray'))
    def updatefig(j):
        # set the data in the axesimage object
        im.set_array(imagelist[j])
        # return the artists set
        return [im]
    
    ani = animation.FuncAnimation(fig_animation, updatefig, frames=range(len(imagelist)), interval=25, blit=True, repeat = False)
    plt.show()
    return ani

def mc_run(N,n_eq,n_mc,T):
    """
    Inputs:
        N - lattice size (NxN)
        n_eq - average number of spin flip attempts per site to come to equilibrium
        n_mc - average number of spin flip attemots per site to fluctuate
        T - temperature
        
    Outputs:
        E - energy per fluctuation
        M - magnetization per fluctuation
        E_eq - energy per equilibration flip
        M_eq - magnetization per equilibration flip
    """
    beta = 1.0/T
    state = initiate_state(N)      # Create random state
    energy = state_hamiltonian(state)
    magnetization = state_magnetization(state)
    E = []
    M = []
    E_eq = []
    M_eq = []    
    imagelist = []
    
    # Come to equilibrium
    for i in range(n_eq*(N**2)):
        state, deltaE, deltaM = spin_flip(state,beta)       # Spin flip
        energy += deltaE
        magnetization += deltaM
        E_eq.append(energy / N**2)                          # Record energy
        M_eq.append(magnetization / N**2)                   # Record magnetization
        if (i%(N**2)==0):
            imagelist.append(state.copy())

    
    # Collect statistics
    for i in range(n_mc*(N**2)):
        state, deltaE, deltaM = spin_flip(state,beta)     # Spin flip
        energy += deltaE
        magnetization += deltaM
        E.append(energy / N**2)                             # Record energy
        M.append(magnetization / N**2)                      # Record magnetization
        if (i%(N**2)==0):
            imagelist.append(state.copy())   

    
    return E,M,E_eq,M_eq,imagelist
    