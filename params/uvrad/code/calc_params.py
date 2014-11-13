import numpy as np

h = 6.626068e-27
c = 2.99792458e10
k = 1.3806503e-16

Mb = 1.0e-18
eV = 1.60217646e-12

# limits of integration
lambdaHminus = h * c / 0.755 / eV #1.6e4 Angstroms
lambdaH2 = h * c / 11.2 / eV      #1.1e3
lambdaH = h * c / 13.6 / eV       #9.1e2
lambdaHe = h * c / 24.59 / eV     #5.0e2
lambdaHeplus = h * c / 54.42 / eV #2.3e2
lambdaMax = h * c / 50000 / eV    #2.5e-2

E_H2 = 11.2 * eV
E_H = 13.6 * eV
E_He = 24.59 * eV
E_Heplus = 54.42 * eV

