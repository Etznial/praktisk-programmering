import numpy as np

def ME(A, Z):
    if not 'MEs' in dir(): # load data file only once
        Ns, Zs, As, MEs, BAs = np.loadtxt('masses.dat', unpack=True)
    return MEs[(Zs == Z) & (As == A)][0]

print("Mass excess of alpha particle:", ME(4, 2), "keV")

def QB(A, Z):
    return ME(A, Z) - ME(A, Z + 1)

print("3H -> 3He + e + ν. Q-value:", QB(3, 1), "keV")

print("235U -> 87Br + 145La + 3n. Q-value:", (ME(235, 92) - (ME(87, 35) + ME(145, 57) + 3*ME(1, 0)))/1e3, "MeV")
