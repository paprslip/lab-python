import numpy as np
from matplotlib import pyplot as plt
from scipy import constants

N = 320 # N is number of turns in the coil
R = 7.25e-2 # radius of the coil in cm

class Coil:
    def __init__(self, frequency, current):
        self.frequency = frequency*10**6                            # freq are inputed in MHz
        self.current = current/2
        self.bfield = (4/5)**(3/2)*(4e-7*np.pi)*N*self.current/R    # B field is aligned with Z axis, thus magnitude of B = Bz
        self.gamma = self.frequency*2*np.pi/self.bfield             # gyromagnetic ratio (ratio of magnetic moment to angular momentum)
        self.lande = 2*constants.m_e*self.gamma/constants.e         # landé g factor
        self.energy = constants.h*self.frequency
        self.energyp = 0.5*self.gamma*constants.hbar*self.bfield
        self.energyn = -0.5*self.gamma*constants.hbar*self.bfield

    def std(self):
        return np.std(self.lande)

    def plot(self):
        plt.title('Angular Frequency vs. Magnetic Field Strength')
        plt.xlabel('Magnetic Field Strength [T]')
        plt.ylabel('Angular Frequency [rad/s]')
        plt.legend(loc="lower right")
        plt.show()

    def plot_nu(self):
        plt.plot(self.bfield, 2*np.pi*self.frequency, 'ob', label='Data')
        plt.errorbar(self.bfield, 2*np.pi*self.frequency, xerr=np.std(self.bfield), yerr=np.std(2*np.pi*self.frequency), fmt="o", color="tab:blue")

        p, cov = np.polyfit(self.bfield, 2*np.pi*self.frequency, 1, cov=True)
        np.polyval(p,self.bfield)
        plt.plot(self.bfield, p[0]*self.bfield + p[1], '--', color='tab:orange', label='Linear Fit')
        #print('Slope: ', a, 'Intercept: ', b)
        print('Slope + Intercept: ', p)
        print('Cov:', np.sqrt(np.diag(cov)))
        plt.annotate('Slope: ' + str(p[0]) + ' +/- ' + str(np.sqrt(np.diag(cov))[0]), xy=(0.05, 0.9), xycoords='axes fraction')
        plt.annotate('Intercept: ' + str(p[1]) + ' +/- ' + str(np.sqrt(np.diag(cov))[1]), xy=(0.05, 0.85), xycoords='axes fraction')

    def chi_square(self):
        p, cov = np.polyfit(self.bfield, 2*np.pi*self.frequency, 1, cov=True)
        np.polyval(p,self.bfield)

        chi_square = 0
        for i in range(len(self.bfield)):
            chi_square += (2*np.pi*self.frequency[i] - (p[0]*self.bfield[i] + p[1]))**2/(np.std(2*np.pi*self.frequency))**2
        print('Chi Square: ', chi_square)

    def r_squared(self):
        p, cov = np.polyfit(self.bfield, 2*np.pi*self.frequency, 1, cov=True)
        np.polyval(p,self.bfield)

        y_bar = np.average(2*np.pi*self.frequency)
        ss_tot = 0
        ss_res = 0
        for i in range(len(self.bfield)):
            ss_tot += (2*np.pi*self.frequency[i] - y_bar)**2
            ss_res += (2*np.pi*self.frequency[i] - (p[0]*self.bfield[i] + p[1]))**2
        r_squared = 1 - ss_res/ss_tot
        print('R Squared: ', r_squared)
        plt.annotate('R Squared: ' + str(r_squared), xy=(0.05, 0.95), xycoords='axes fraction')

    def residuals(self):
        p, cov = np.polyfit(self.bfield, 2*np.pi*self.frequency, 1, cov=True)
        np.polyval(p,self.bfield)

        plt.plot(self.bfield, 2*np.pi*self.frequency - (p[0]*self.bfield + p[1]), 'ob')
        plt.errorbar(self.bfield, 2*np.pi*self.frequency - (p[0]*self.bfield + p[1]), xerr=np.std(self.bfield), yerr=np.std(2*np.pi*self.frequency), fmt="o", color="tab:blue")
        plt.title('Residuals')
        plt.xlabel('Magnetic Field Strength [T]')
        plt.ylabel('Residuals [rad/s]')
        plt.axhline(y=0, color='tab:orange', linestyle='-')
        plt.show()

big = Coil(np.array([20.429,23.057,27.983,29.629,32.642]), np.array([0.37,0.424,0.504,0.555,0.622]))              # [0.185,0.212,0.252,0.2775,0.311]
med = Coil(np.array([41.557,44.315,47.449,52.781,56.047]), np.array([0.718,0.791,0.873,0.955,1.006]))
small = Coil(np.array([49.216,50.521,50.864,51.761,52.638]), np.array([0.925,0.931,0.944,0.939,0.969]))

print('Spin Gyromagnetic Ratio: ', np.average(small.gamma), '+/-', np.std(small.gamma))
small.plot_nu()
small.chi_square()
small.r_squared()
print('Lande g-factor', np.average(small.lande), '+/-', np.std(small.lande))
small.plot()
small.residuals()
print(small.bfield)
## Landé g factors
#print(big.lande)
#print(med.lande)
#print(small.lande)
#
## Plot energy vs field
#plt.title('Energy vs. Magnetic Field Strength')
#plt.xlabel('Magnetic Field Strength')
#plt.ylabel('Potential energy Energy of Free Electron')
#plt.plot(big.bfield, big.energyp, 'or')
#plt.plot(big.bfield, big.energyn, 'og')
#plt.plot(big.bfield, big.energy, 'ob')
#
## Lines of best fit
#a1, b1 = np.polyfit(big.bfield, big.energyp, 1)
#c1, d1 = np.polyfit(big.bfield, big.energyn, 1)
#a2, b2 = np.polyfit(med.bfield, med.energyp, 1)
#c2, d2 = np.polyfit(med.bfield, med.energyn, 1)
#a3, b3 = np.polyfit(small.bfield, small.energyp, 1)
#c3, d3 = np.polyfit(small.bfield, small.energyn, 1)
#
#plt.plot(big.bfield, a1*big.bfield + b1, '--r')
#plt.plot(big.bfield, c1*big.bfield + d1, '--r')

#plt.show()

