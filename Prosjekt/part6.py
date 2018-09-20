from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, Boltzmann, Avogadro
plt.style.use("ggplot")

Amu_to_kg = 1.66054e-27

wavelengths, flux = np.transpose(np.load("spectrum_data.npy"))
_, flux_standard_deviation = np.transpose(np.load("sigma_noise.npy"))
measured_flux = flux + flux_standard_deviation



gases = {'O2':{'lines': [630, 690, 760], 'weight':32 },
         'H2O':{'lines': [720, 820, 940], 'weight':18 },
         'CO2':{'lines': [1400, 1600], 'weight': 44},
         'CH4':{'lines': [1660, 2200], 'weight':16 },
         'CO':{'lines': [2340], 'weight': 28},
         'N2O':{'lines': [2870], 'weight': 30}}

def get_flux_min_range():
    return np.linspace(0.7, 1.0, 30) # [lines, 30]

def get_lambda_center_range(name_of_gas):
    v_max = 10000           #[m/s]
    lambda_0 = gases[name_of_gas]["lines"]
    lambda_0_range = np.zeros((len(lambda_0), 300))
    for line in xrange(len(lambda_0)):
        max_doppler_shift = lambda_0[line] * v_max / speed_of_light
        wavelengths_within_max_shift_mask = np.abs(wavelengths - lambda_0[line]) < max_doppler_shift
        wavelengths_within_max_shift = wavelengths[wavelengths_within_max_shift_mask]
        lambda_0_range[line] = np.linspace(wavelengths_within_max_shift[0], wavelengths_within_max_shift[-1], 300)

    return lambda_0_range #[lines, 300]

def get_sigma_range(name_of_gas, T_min=150, T_max=450):
    T = np.linspace(T_min, T_max, 30)
    k = Boltzmann
    c = speed_of_light
    lambda_0 = gases[name_of_gas]["lines"]
    m = gases[name_of_gas]["weight"] * Amu_to_kg
    sigma_range = np.zeros((len(lambda_0), 30))
    for line in xrange(len(lambda_0)):
        sigma_range[line] = lambda_0[line] * np.sqrt(k*T/m) / c

    return sigma_range #[lines, 30]
"""
def get_mask(line, lambda_center_range):
    v_max = 10000
    max_doppler_shift = lines[line] * v_max / speed_of_light
    max_wl = np.max(lambda_center_range) + 4 * max_doppler_shift
    min_wl = np.min(lambda_center_range) - 4 * max_doppler_shift
    mask1 = wavelengths < max_wl
    mask2 = wavelengths > min_wl
    return np.logical_and(mask1, mask2)
    """
def get_mask(lambda_center_range):
    mask1 = wavelengths < lambda_center_range[-1]
    mask2 = wavelengths > lambda_center_range[0]
    return np.logical_and(mask1, mask2)

def get_mask2(lambdsaas, min, mac):
    mask1 = lambdsaas < mac
    mask2 = lambdsaas > min
    return np.logical_and(mask1, mask2)

def flux_model(lambda_values, sigma, lambda_center, flux_min):
    f_max = 1
    f_model = f_max + (flux_min - f_max)*np.exp(-(lambda_values-lambda_center)**2/(2*sigma**2))
    return f_model

def chi_squared(sigma_range, lambda_center_range, flux_min_range):
    chi = np.zeros((len(sigma_range), len(sigma_range[0]), len(lambda_center_range[0]), len(flux_min_range)))
    for line in xrange(len(sigma_range)):
        print line
        mask = get_mask(lambda_center_range[line])
        lambda_values = wavelengths[mask]
        for sigma in xrange(len(sigma_range[line])):
            print "sigma = ", sigma
            for lambda_center in xrange(len(lambda_center_range[line])):
                if lambda_center%20 == 0:
                    print "lambda_center = ", lambda_center
                for f_min in xrange(len(flux_min_range)):
                    f_model = flux_model(lambda_values, sigma_range[line, sigma], lambda_center_range[line, lambda_center], flux_min_range[f_min])
                    val = (measured_flux[mask] - f_model)**2/flux_standard_deviation[mask]**2
                    chi[line, sigma, lambda_center, f_min] = np.sum(val)

    return chi


def plot_chi_square():
    chi = np.load("chi_O2.npy")
    for line in range(len(chi)):
            print np.unravel_index(np.argmin(chi[line]), chi[line].shape)

    l_c_range = get_lambda_center_range("O2")
    sigma_range = get_sigma_range("O2")
    f_min_range = get_flux_min_range()

    sigma = sigma_range[2, 0]
    l_c = l_c_range[2, 0]
    f_min = f_min_range[26]

    mask = get_mask2(wavelengths, l_c - 8*sigma, l_c + 8*sigma)
    lambda_values = wavelengths[mask]
    f_model = flux_model(lambda_values, sigma, l_c, f_min)
    val = (measured_flux[mask] - f_model)**2/flux_standard_deviation[mask]**2

    ax = plt.subplot()
    legend = ["model", "obs", "$\lambda$_center = %f"%l_c]

    ax.plot(lambda_values, f_model)
    ax.plot(lambda_values, measured_flux[mask])

    ax.scatter(lambda_values, measured_flux[mask])

    plt.show()

    #print sigma_range[0, 29]
    #print l_c_range[0, 11]
    #print f_min_range[25]
"""
def plaatt():

    plt.plot(meget_godt_nedkortet_versjon_av_wavelengths, meget_godt_nedkortet_versjon_av_flux)
    plt.vlines(2869.91143779, 0, 1)
    plt.show()
    #print sigma, lambda_center, f_min

"""
#x = np.linspace(-10,10,1000)
#y = (3*x**2 + 4*x) / (3*x**4 - x**2)
#plt.plot(x, y)
#plt.show()

plot_chi_square()

#plot_chi_square()
#chi = chi_squared(get_sigma_range("H2O"), get_lambda_center_range("H2O"), get_flux_min_range())
#np.save("chi_H2O", chi)
