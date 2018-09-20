from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, Boltzmann, Avogadro
plt.style.use("ggplot")

Amu_to_kg = 1.66054e-27

wavelengths, flux = np.transpose(np.load("spectrum_data.npy"))
_, flux_standard_deviation = np.transpose(np.load("sigma_noise.npy"))


gases = {'O2':{'lines': [630, 690, 760], 'weight':32 },
         'H20':{'lines': [720, 820, 940], 'weight':18 },
         'CO2':{'lines': [1400, 1600], 'weight': 44},
         'CH4':{'lines': [1660, 2200], 'weight':16 },
         'CO':{'lines': [2340], 'weight': 28},
         'N2O':{'lines': [2870], 'weight': 30}}

def sigma_range(name_of_gas, T_min=150, T_max=450):
    T = np.linspace(T_min, T_max, 30)
    k = Boltzmann
    c = speed_of_light
    lambda_0 = np.reshape(gases[name_of_gas]["lines"], (3, 1))
    m = gases[name_of_gas]["weight"] * Amu_to_kg

    sigma = lambda_0 * np.sqrt(k*T/m) / c
    return sigma

def flux_min_range():
    return np.linspace(0.7, 1.0, 30)

def lambda_center_range(name_of_gas):
    v_max = 10000           #[m/s]
    lambda_0 = gases[name_of_gas]["lines"]
    lambda_0_range = np.zeros((len(lambda_0), 300))
    for line in range(len(lambda_0)):
        max_doppler_shift = lambda_0[line] * v_max / speed_of_light
        wavelengths_within_max_shift_mask = np.abs(wavelengths - lambda_0[line]) < max_doppler_shift
        wavelengths_within_max_shift = wavelengths[wavelengths_within_max_shift_mask]
        lambda_0_range[line] = np.linspace(wavelengths_within_max_shift[0], wavelengths_within_max_shift[-1], 300)

    return lambda_0_range

def spectral_line_model():
    F_max = 1

def xhi_squared(sigma_range, flux_min_range, lambda_center_range):
    xhi = np.zeros((len(lambda_center_range), len(lambda_center_range[0]), len(flux_min_range), len(sigma_range[0])))
    for line in range(len(lambda_center_range[:])):
        max_shift = np.max(lambda_center_range[line]*10000/speed_of_light)
        max_wl = np.max(lambda_center_range[line]) + 4 * max_shift
        min_wl = np.min(lambda_center_range[line]) - 4 * max_shift
        mask1 = wavelengths < max_wl
        mask2 = wavelengths > min_wl
        mask = np.logical_and(mask1, mask2)
        lambda_values = wavelengths[mask]
        for sigma in range(len(sigma_range[0])):
            for f_min in range(len(flux_min_range)):
                for l_c in range(len(lambda_center_range[0])):
                    max_wl = np.max(lambda_center_range[line, l_c]) + 4 * max_shift
                    min_wl = np.min(lambda_center_range[line, l_c]) - 4 * max_shift
                    mask1 = lambda_values < max_wl
                    mask2 = lambda_values > min_wl
                    mask = np.logical_and(mask1, mask2)
                    l_values = lambda_values[mask]
                    flux_max = 1
                    flux_model = flux_max + (flux_max - flux_min_range[f_min])*np.exp(-(l_values - lambda_center_range[line, l_c])**2/(2*sigma_range[line, sigma]**2))
                    val = (flux[mask] - flux_model)**2 / flux_standard_deviation[mask]**2
                    xhi[line, l_c, f_min, sigma] = np.sum(val)

    print xhi
xhi_squared(sigma_range("O2"), flux_min_range(), lambda_center_range("O2"))