import glob,matplotlib as mpl, numpy as np,pandas as pd
import os
from Photonic_fit_funtction_lib_V2021_10 import *
from Photonic_instruments_lib_V2023_07 import *
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

"""
Calculates peaks and fittes them
Input:
    x: e.g. wavelength
    y: e.g. Intensity
    expected_peaks: Number of expected peaks to be found -> Returns that many peaks
    fit_type: Type of fit that shall be used
        If None: No fit, returns detected peaks
        Alternatives: Currently voigtlin_same_cen_const
    plot: 
        If True: Plots data and fit function
        Else: No plot
    
Return: Detected peaks

Example: 
p0, p1, p2 = calculate_peaks(np.asarray(df.loc[:, "Wavelength (nm)"]), np.asarray(df.loc[:, "Detector 1"]), expected_peaks = 3, fit_type = voigtlin_same_cen_const, plot = False)
"""
def calculate_peaks(x, y, expected_peaks = 3, fit_type = None, plot = False):
    dx = abs(x[0]-x[101])/100
    dy = abs(y[0]-np.max(y))/1000
    peaks, properties = find_peaks(y, height = 0.5, distance = 1/dx, width=1/dy)

    if plot == True:
        load_default_plot_font(fontsize = 20)

    if fit_type != None:
        width = properties["widths"]
        amp = properties["peak_heights"]
        cen = x[peaks]
        left_base = properties['left_ips']
        right_base = properties['right_ips']

        peaks_fit = []
        for i in range(0, len(peaks)):
            if fit_type == voigtlin_same_cen_const:
                try: 
                    x_fit = x[int(left_base[i]-width[i]*2):int(right_base[i]+width[i]*2)]
                    y_fit = y[int(left_base[i]-width[i]*2):int(right_base[i]+width[i]*2)]
                    p0 = [amp[i], 1, cen[i], amp[i], width[i]]
                    popt, pcov = curve_fit(fit_type, x_fit, y_fit, p0=p0)
                    peak_fit, properties = find_peaks(fit_type(x, *popt), height = 0.75, distance = 1/dx, width=1/dy)
                    if plot == True:
                        plt.plot(x, y, x_fit, fit_type(x_fit, *popt))
                        plt.title('Fit')
                    if abs(x[peaks[i]]-x[peak_fit[0]]) < 0.5:
                        peaks_fit.append(peak_fit[0])
                    else:
                        peaks_fit.append(peaks[i])
                        print('No Fit')
                except:
                    peaks_fit.append(peaks[i])
                    print('No Fit')

    if plot == True and fit_type != None:
        plt.xlabel('Wavelength / nm')
        plt.ylabel('Intensity')
        plt.title('Fit of Peaks')
        plt.show()
    if fit_type != None:
        peaks = peaks_fit

    if expected_peaks == 3:
        if len(peaks) == 3:
            peak2 = x[peaks[0]]
            peak1 = x[peaks[1]]
            peak0 = x[peaks[2]]
        elif len(peaks) == 2:
            peak2 = x[peaks[0]]
            peak1 = x[peaks[1]]
            peak0 = 0
        elif len(peaks) == 1:
            peak2 = x[peaks[0]]
            peak1 = 0
            peak0 = 0
        else:
            print(f'AMOUNT OF PEAKS NOT MATCH AMOUNT OF FBGs! {len(peaks)}')
            return 0,0,0
        return peak0, peak1, peak2
    
    elif expected_peaks == 1:
        if len(peaks) == 1:
            return x[peaks[0]]
        else:
            print(f'AMOUNT OF PEAKS NOT MATCH AMOUNT OF FBGs! {len(peaks)}')
            return 0


"""
Calculates wavelength to temperature using calculated parameters
Input:
    peaks: Calculated peaks
    parms: Parameters from polynomial fit function
    reference_value: e.g. ice point value if fit function is normed
    polynomial_deg: degree of polynominal function used
        Second and third are currently included
    plot: 
        If True: Plots temperature over wavelength difference (to reference value)
        Else: No plot
    
Return: Temperature values

Example: 
parms = [7.29485246, 90.4193956, -1.78066693]
t0 = wavelength_to_temperature(np.array(peaks0), parms, 1553.5102583333335, polynomial_deg = 2, plot = True)
"""
def wavelength_to_temperature(peaks, parms, reference_value = 0, polynomial_deg = 2, plot = False):
    def poly2_fit(x, a, b, c):
        return a + b*x + c*x**2
    def poly3_fit(x, a, b, c, d):
        return a + b*x + c*x**2 + d*x**3

    if polynomial_deg == 2:
        temperature = poly2_fit(peaks-reference_value, parms[0], parms[1], parms[2])
    if polynomial_deg == 3:
        temperature = poly3_fit(peaks-reference_value, parms[0], parms[1], parms[2], parms[3])

    if plot == True:
        plt.plot(peaks-reference_value, temperature, label = f'Fit for {reference_value:.2f}')
        plt.title(f'FBG fits')
        plt.xlabel('Wavelength difference / nm')
        plt.ylabel('Temperature / °C')
        plt.legend(loc='best')

    return temperature
