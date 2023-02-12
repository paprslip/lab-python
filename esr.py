import numpy as np
from matplotlib import pyplot as plt
from scipy import constants

N = 320 # N is number of turns in the coil
R = 7.25e-2 # radius of the coil in cm

class Coil:
    def __init__(self, frequency, current):
        self.frequency = frequency*10**6                    # freq are inputed in MHz
        self.current = current/2
        self.bfield = (4/5)**(3/2)*(4e-7*np.pi)*N*self.current/R # B field is aligned with Z axis, thus magnitude of B = Bz
        self.gamma = self.frequency*2*np.pi/self.bfield          # gyromagnetic ratio (ratio of magnetic moment to angular momentum)
        self.lande = 2*constants.m_e*self.gamma/constants.e # landé g factor
        self.energy = constants.h*self.frequency
        self.energyp = 0.5*self.gamma*constants.hbar*self.bfield
        self.energyn = -0.5*self.gamma*constants.hbar*self.bfield

big = Coil(np.array([20.429,23.057,27.983,29.629,32.642]), np.array([0.37,0.424,0.504,0.555,0.622]))              # [0.185,0.212,0.252,0.2775,0.311]
med = Coil(np.array([30.321,41.557,44.315,47.449,52.781,56.047,78.148]), np.array([0.531,0.718,0.791,0.873,0.955,1.006,1.319]))
small = Coil(np.array([49.216,50.521,50.864,51.761,52.638]), np.array([0.925,0.931,0.944,0.939,0.969]))

# Landé g factors
print(big.lande)
print(med.lande)
print(small.lande)

# Plot energy vs field
plt.title('Energy vs. Magnetic Field Strength')
plt.xlabel('Magnetic Field Strength')
plt.ylabel('Potential energy Energy of Free Electron')
plt.plot(big.bfield, big.energyp)
plt.plot(big.bfield, big.energyn)
plt.plot(big.bfield, big.energy)
plt.show()

