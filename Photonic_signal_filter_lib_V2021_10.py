# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 18:11:10 2020

@author: eiserm01
"""

import numpy as np, matplotlib.pyplot as plt
import os
from scipy.signal import butter, sosfiltfilt, sosfreqz

def plot_filtered_spektrum(wavelength,signal,signal_filt,fs,cutoff,savepath=False,savename=False,x_min=False,x_max=False):
    plt.figure(figsize=(9,6),dpi=300)
    plt.plot(wavelength,signal,"--",color="black",label='signal')
    plt.plot(wavelength,signal_filt,":",color="red",label='signal filtered (fc='+str(cutoff)+'Hz)')
    plt.xlabel("wavelength / nm")
    plt.ylabel("signal amplitude / V")
    plt.ylim(0.0,max(signal))
    if not x_min==False:
        plt.xlim(x_min,x_max)
    plt.legend(loc='best')
    plt.tight_layout()
    if os.path.exists(savepath):
        #plt.title(savename+' fs='+str(fs)+'Hz')
        plt.tight_layout()
        #plt.savefig(savepath+savename+"_signal_filt_V3",dpi=300)
    else:
        print ("path don't exist")

def plot_filter(fs=1,lowcut=1,highcut=100,filtertyp="low"):
        plt.figure(1)
        plt.clf()
        for order in [3, 6, 9]:
            if filtertyp=="low":
                sos = butter_lowpass(lowcut, fs, order=order)
            elif filtertyp=="high":
                sos = butter_highpass(highcut, fs, order=order)
            elif filtertyp=="band":
                sos = butter_bandpass(lowcut,highcut, fs, order=order)
            else:
                plt.close()
                return "Please chose valid filter: low, high, band"
            w, h = sosfreqz(sos, worN=2000)
            plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)
        
        plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
                 '--', label='sqrt(0.5)')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Gain')
        plt.grid(True)
        plt.legend(loc='best')  

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5, axis=-1):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y

def butter_lowpass(lowcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        sos = butter(order, [low], analog=False, btype='lowpass', output='sos')
        return sos

def butter_lowpass_filter(data, lowcut, fs, order=5, axis=-1):
        sos = butter_lowpass(lowcut, fs, order=order)
        y = sosfiltfilt(sos, data)
        return y
    
def butter_highpass(highcut, fs, order=5):
        nyq = 0.5 * fs
        high = highcut / nyq
        sos = butter(order, [high], analog=False, btype='highpass', output='sos')
        return sos

def butter_highpass_filter(data, highcut, fs, order=5, axis=-1):
        sos = butter_highpass(highcut, fs, order=order)
        y = sosfiltfilt(sos, data, axis)
        return y