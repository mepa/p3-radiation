#!/usr/bin/python

import numpy as np
import sys
import getopt
import math

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

def sq(x): return x*x

# photoionization cross sections from Verner et al. 1996
def sigma_H(wavelength):
    if(wavelength > lambdaH): return 0.0
    E = h * c / wavelength
    sigma0 = 5.475e4 * Mb
    E0 = 0.4298 * eV
    yw = 0
    ya = 32.88
    P = 2.963
    y0 = 0
    y1 = 0
    x = E / E0 - y0
    y = np.sqrt(sq(x) + sq(y1))
    F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + np.sqrt(y / ya), -P)
    return sigma0 * F

def sigma_He(wavelength):
    if(wavelength > lambdaHe): return 0.0
    E = h * c / wavelength
    sigma0 = 9.492e2 * Mb
    E0 = 13.61 * eV
    yw = 2.039
    ya = 1.469
    P = 3.188
    y0 = 0.4434
    y1 = 2.136
    x = E / E0 - y0
    y = np.sqrt(sq(x) + sq(y1))
    F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + np.sqrt(y / ya), -P)
    return sigma0 * F

def sigma_Heplus(wavelength):
    if(wavelength > lambdaHeplus): return 0.0
    E = h * c / wavelength
    sigma0 = 1.369e4 * Mb
    E0 = 1.72 * eV
    yw = 0
    ya = 32.88
    P = 2.963
    y0 = 0
    y1 = 0
    x = E / E0 - y0
    y = np.sqrt(sq(x) + sq(y1))
    F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + np.sqrt(y / ya), -P)
    return sigma0 * F

def sigma_Hminus(wavelength):
    if(wavelength > lambdaHminus): return 0.0
    E = h * c / wavelength
    E_in_eV = E / eV
    return 2.1e-16 * pow(E_in_eV - 0.75, 1.5) / pow(E_in_eV, 3.11)

def sigma_H2(wavelength):
    if(wavelength > lambdaH2): return 0.0
    return 3.78e-18

def readHeader(f, age, wavelength_points):
    age = f.readline().split(":")[1].rstrip()
    f.readline()
    f.readline()
    wavelength_points = f.readline().split(":")[1].rstrip()
    f.readline()
    f.readline()
    f.readline()
    return age, wavelength_points
    

def main(argv):
    infile = ""
    outfile = ""
    try:
        opts, args = getopt.getopt(argv,'ha:i:o:',['tage=','ifile=','ofile='])
    except getopt.GetoptError:
        print "usage: calc_params.py -a <age_in_Myr> -i <infile> -o <outfile>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "calc_params.py -a <age_in_Myr> -i <infile> -o <outfile>"
            sys.exit()
        elif opt in ('-a', '--tage'):
            target_age = int(arg)
        elif opt in ('-i', '--ifile'):
            infile = arg
        elif opt in ('-o', '--ofile'):
            outfile = arg

    age = 0.0
    wavelength_points = 0
    with open(infile) as f:
        f.readline()
        f.readline()

        age, wavelength_points = readHeader(f, age, wavelength_points)
        age = float(age)
        wavelength_points = int(wavelength_points)

        wavelength = np.zeros(wavelength_points)
        flux = np.zeros(wavelength_points)
        spectrum_age = -1.0
        while(age < target_age):
            for i in range(wavelength_points):
                line = f.readline()
                wavelength[i] = float(line.split()[0]) * 1.0e-8
                flux[i] = float(line.split()[1]) * 1.0e8

            spectrum_age = age
      
            f.readline()
            age, wavelength_points = readHeader(f, age, wavelength_points)
            age = float(age)
            wavelength_points = int(wavelength_points)

    print ""

    print "target_age = ", target_age, " Myr"
    print "spectrum_age = ", spectrum_age, " Myr"
 
    bins = 5

    bin_max = [lambdaHminus, lambdaH2, lambdaH, lambdaHe, lambdaHeplus] 
    bin_min = [lambdaH2, lambdaH, lambdaHe, lambdaHeplus, lambdaMax]
    
    N_dot_Hminus = [0.0]*bins
    N_dot_H = [0.0]*bins
    N_dot_He = [0.0]*bins
    N_dot_Heplus = [0.0]*bins
    N_dot_H2 = [0.0]*bins
    
    Gamma_Hminus = [0.0]*bins
    Gamma_H = [0.0]*bins
    Gamma_He = [0.0]*bins
    Gamma_Heplus = [0.0]*bins
    
    GammaE_H = [0.0]*bins
    GammaE_He = [0.0]*bins
    GammaE_Heplus = [0.0]*bins
    
    for i in range(1,wavelength_points):
        delta_wavelength = wavelength[i] - wavelength[i - 1]
        wavelength_1 = wavelength[i - 1]
        wavelength_2 = wavelength[i]
        E_1 = h * c / wavelength_1
        E_2 = h * c / wavelength_2
        flux_1 = flux[i - 1]
        flux_2 = flux[i]
        N_dot_1 = flux_1 / E_1
        N_dot_2 = flux_2 / E_2

        for bin in range(bins):
            if wavelength[i - 1] >= bin_min[bin] and wavelength[i] <= bin_max[bin]:
 
                N_dot_Hminus[bin] += delta_wavelength * 0.5 * (N_dot_1 + N_dot_2)
                Gamma_Hminus[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_Hminus(wavelength_1) + N_dot_2 * sigma_H(wavelength_2))

                N_dot_H2[bin] += delta_wavelength * 0.5 * (N_dot_1 + N_dot_2) #LW photons

                N_dot_H[bin] += delta_wavelength * 0.5 * (N_dot_1 + N_dot_2)
                Gamma_H[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_H(wavelength_1) + N_dot_2 * sigma_H(wavelength_2))	 
                GammaE_H[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_H(wavelength_1) * (E_1 - E_H) + N_dot_2 * sigma_H(wavelength_2) * (E_2 - E_H))	
                
                N_dot_He[bin] += delta_wavelength * 0.5 * (N_dot_1 + N_dot_2)
                Gamma_He[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_He(wavelength_1) + N_dot_2 * sigma_He(wavelength_2))
                GammaE_He[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_He(wavelength_1) * (E_1 - E_He) + N_dot_2 * sigma_He(wavelength_2) * (E_2 - E_He))  
                
                N_dot_Heplus[bin] += delta_wavelength * 0.5 * (N_dot_1 + N_dot_2)
                Gamma_Heplus[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_Heplus(wavelength_1) + N_dot_2 * sigma_Heplus(wavelength_2)) 
                GammaE_Heplus[bin] += delta_wavelength * 0.5 * (N_dot_1 * sigma_Heplus(wavelength_1) * (E_1 - E_Heplus) + N_dot_2 * sigma_Heplus(wavelength_2) * (E_2 - E_Heplus))

    print ""

    print "N_dot_Hminus = ", N_dot_Hminus[0], N_dot_Hminus[1], N_dot_Hminus[2], N_dot_Hminus[3], N_dot_Hminus[4]
    print "N_dot_H2 = ", N_dot_H2[0], N_dot_H2[1], N_dot_H2[2], N_dot_H2[3], N_dot_H2[4]
    print "N_dot_H = ", N_dot_H[0], N_dot_H[1], N_dot_H[2], N_dot_H[3], N_dot_H[4]
    print "N_dot_He = ", N_dot_He[0], N_dot_He[1], N_dot_He[2], N_dot_He[3], N_dot_He[4]
    print "N_dot_Heplus = ", N_dot_Heplus[0], N_dot_Heplus[1], N_dot_Heplus[2], N_dot_Heplus[3], N_dot_Heplus[4]
    
    sigma_Hminus_average = [0.0]*bins
    sigma_H_average = [0.0]*bins
    sigma_He_average = [0.0]*bins
    sigma_Heplus_average = [0.0]*bins
    
    sigmaE_H_average = [0.0]*bins
    sigmaE_He_average = [0.0]*bins
    sigmaE_Heplus_average = [0.0]*bins
    
    DeltaE_H_average = [0.0]*bins
    DeltaE_He_average = [0.0]*bins
    DeltaE_Heplus_average = [0.0]*bins
    
    for bin in range(bins):
        sigma_Hminus_average[bin] = Gamma_Hminus[bin] / N_dot_Hminus[bin] if  N_dot_Hminus[bin] > 0.0 else 0.0
        sigma_H_average[bin] = Gamma_H[bin] / N_dot_H[bin] if N_dot_H[bin] > 0.0 else 0.0
        sigma_He_average[bin] = Gamma_He[bin] / N_dot_He[bin] if N_dot_He[bin] > 0.0 else 0.0
        sigma_Heplus_average[bin] = Gamma_Heplus[bin] / N_dot_Heplus[bin] if N_dot_Heplus[bin] > 0.0 else 0.0
        
        sigmaE_H_average[bin] = GammaE_H[bin] / N_dot_H[bin] if N_dot_H[bin] > 0.0 else 0.0
        sigmaE_He_average[bin] = GammaE_He[bin] / N_dot_He[bin] if N_dot_He[bin] > 0.0 else 0.0
        sigmaE_Heplus_average[bin] = GammaE_Heplus[bin] / N_dot_Heplus[bin] if N_dot_Heplus[bin] > 0.0 else 0.0
        
        DeltaE_H_average[bin] = sigmaE_H_average[bin] / sigma_H_average[bin] if sigma_H_average[bin] > 0.0 else 0.0
        DeltaE_He_average[bin] = sigmaE_He_average[bin] / sigma_He_average[bin] if sigma_He_average[bin] > 0.0 else 0.0
        DeltaE_Heplus_average[bin] = sigmaE_Heplus_average[bin] / sigma_Heplus_average[bin] if sigma_Heplus_average[bin] > 0.0 else 0.0
        
        DeltaE_H_average[bin] /= eV
        DeltaE_He_average[bin] /= eV
        DeltaE_Heplus_average[bin] /= eV

    #Gamma_H_reference = 1.0e-12
    Gamma_H_reference = 0.3e-12

    Gamma_H_total_old = Gamma_H[2] + Gamma_H[3] + Gamma_H[4]
    
    factor = Gamma_H_reference / Gamma_H_total_old
    
    N_dot_H_total = (N_dot_H[2] + N_dot_H[3] + N_dot_H[4]) * factor
    N_dot_He_total = (N_dot_He[3] + N_dot_He[4]) * factor
    N_dot_Heplus_total = (N_dot_Heplus[4]) * factor
    N_dot_Hminus_total = (N_dot_Hminus[0] + N_dot_Hminus[1]) * factor
    N_dot_H2_total = (N_dot_H2[1]) * factor
    
    Gamma_H_total = (Gamma_H[2] + Gamma_H[3] + Gamma_H[4]) * factor
    Gamma_He_total = (Gamma_He[3] + Gamma_He[4]) * factor
    Gamma_Heplus_total = (Gamma_Heplus[4]) * factor
    
    print ""
    
    print "sigma_Hminus = ", sigma_Hminus_average[0], sigma_Hminus_average[1], sigma_Hminus_average[2], sigma_Hminus_average[3], sigma_Hminus_average[4]
    print "sigma_H = ", sigma_H_average[0], sigma_H_average[1], sigma_H_average[2], sigma_H_average[3], sigma_H_average[4]
    print "sigma_He = ", sigma_He_average[0], sigma_He_average[1], sigma_He_average[2], sigma_He_average[3], sigma_He_average[4]
    print "sigma_Heplus = ", sigma_Heplus_average[0], sigma_Heplus_average[1], sigma_Heplus_average[2], sigma_Heplus_average[3], sigma_Heplus_average[4]
    
    print ""

    print "sigmaE_H = ", sigmaE_H_average[0], sigmaE_H_average[1], sigmaE_H_average[2], sigmaE_H_average[3], sigmaE_H_average[4]
    print "sigmaE_He = ", sigmaE_He_average[0], sigmaE_He_average[1], sigmaE_He_average[2], sigmaE_He_average[3], sigmaE_He_average[4]
    print "sigmaE_Heplus = ", sigmaE_Heplus_average[0], sigmaE_Heplus_average[1], sigmaE_Heplus_average[2], sigmaE_Heplus_average[3], sigmaE_Heplus_average[4]
    
    print ""

    print "DeltaE_H = ", DeltaE_H_average[0], DeltaE_H_average[1], DeltaE_H_average[2], DeltaE_H_average[3], DeltaE_H_average[4]
    print "DeltaE_He = ", DeltaE_He_average[0], DeltaE_He_average[1], DeltaE_He_average[2], DeltaE_He_average[3], DeltaE_He_average[4]
    print "DeltaE_Heplus = ", DeltaE_Heplus_average[0], DeltaE_Heplus_average[1], DeltaE_Heplus_average[2], DeltaE_Heplus_average[3], DeltaE_Heplus_average[4]
    
    print ""

    print "N_dot_Hminus_total = ", N_dot_Hminus_total
    print "N_dot_H_total = ", N_dot_H_total
    print "N_dot_He_total = ", N_dot_He_total
    print "N_dot_Heplus_total = ", N_dot_Heplus_total
    print "N_dot_H2 = ", N_dot_H2_total
    
    print ""
    
    print "J21_LW = ", 1.0e21 * N_dot_H2_total * 0.5 * (E_H2 + E_H) / (4.0 * np.pi * (E_H - E_H2) / h)
    
    print ""

    print "Gamma_H_total = ", Gamma_H_total
    print "Gamma_He_total = ", Gamma_He_total
    print "Gamma_Heplus_total = ", Gamma_Heplus_total
    
    print ""

    #print N_dot_Heplus[0], N_dot_Heplus[1], N_dot_Heplus[2], N_dot_Heplus[3], N_dot_Heplus[4]

    # ostrstream filename
    # filename, "spec_", spectrum_age, ".out", ends 
    # ofstream spec_out(filename.str())
    # for(size_t i = 0 i < wavelength.size() i ++)
    # {
    #   spec_out, wavelength[i], flux[i] * factor
    # }


if __name__ == "__main__":
    main(sys.argv[1:])


