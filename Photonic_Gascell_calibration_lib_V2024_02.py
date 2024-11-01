# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:40:15 2021

@author: eiserm01
"""
import os, matplotlib.pyplot as plt,numpy as np,pandas as pd
from Photonic_fit_funtction_lib_V2021_10 import *
#from Photonic_instruments_lib_V2023_07 import *
from Photonic_signal_filter_lib_V2021_10 import *
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

#load_default_plot_font(fontsize = 20)

def get_thresh_lorentzianlin(thresh,fit_params,x1,x2):
    fit1thresh = root_scalar(lorentzianlin,\
        args=(fit_params[0],fit_params[1]-thresh,fit_params[2],fit_params[3],fit_params[4],fit_params[5]),\
        bracket=[x1,x2])
    return fit1thresh.root

def peak_fitting_Voigt_2center(Signal,peaks_CH,peak_properties,fit_range_faktor=0.5,plot=False,plot_save=False):
        Peaks_pos_Fit           = np.zeros((len(peaks_CH)))
        Peaks_amp_Fit           = np.zeros((len(peaks_CH)))
                 
        test = np.zeros((len(peaks_CH)))
        
        if plot == True:
            fig, ax = plt.subplots(2, 1,sharex=True,figsize=(15,11)) 
            ax[0].plot(Signal,'k:',linewidth=3,label='Transmission spectra')
            ax[0].plot(peaks_CH,Signal[peaks_CH],'bo',linewidth=3,label='Simple peak detection')
            ax[0].set_title('Voigt with Gaus and Lorentz centre')
    
        for j,peak in enumerate(peaks_CH):
            FWHM              = np.round(peak_properties['right_ips'][j])-np.round(peak_properties['left_ips'][j])
            
            pos1            = int(np.round(peak_properties['left_ips'][j]-(FWHM*fit_range_faktor)))
            pos2            = int(np.round(peak_properties['right_ips'][j]+(FWHM*fit_range_faktor)))
    
            Pos_range       = np.arange(pos1,pos2)
            
            guess_voigtlin                         = [ peak,peak_properties['peak_heights'][j],(FWHM)/2,peak,peak_properties['peak_heights'][j],(FWHM)/2,0,0]   # Center Gaus Lorentz kann unterschiedlich sein
            popt_voigtlin, pcov_voigtlin           = curve_fit(voigtlin, xdata=Pos_range,ydata=-Signal[pos1:pos2]+1,p0=guess_voigtlin,bounds=(0,np.inf),maxfev=100000)      
            
            y_voigtlin= -voigtlin(Pos_range, *popt_voigtlin)+1
                                
            presicion = 0.1 # Subsampling des Peaks
            Pos_range_fine  = np.arange(pos1,pos2,presicion)
            y_fine  = -voigtlin(Pos_range_fine, *popt_voigtlin)+1
            y_peak_voigtlin = y_fine[np.argmin(y_fine)]

            Peaks_pos_Fit[j] = Pos_range_fine[np.argmin(y_fine)]
            Peaks_amp_Fit[j] = y_peak_voigtlin
                    
            if plot == True:
                ax[0].plot(Pos_range_fine,y_fine,'g:')#,label='voigtlin fit')
                ax[1].plot(Pos_range,y_voigtlin-Signal[pos1:pos2],'g:')#,label='voigtlin')
                
        if plot == True:
            ax[0].plot(Peaks_pos_Fit,Peaks_amp_Fit,'ro',label='Peaks of Voigt+Lin.')
            ax[1].plot(Pos_range,y_voigtlin-Signal[pos1:pos2],'g:',label='Voigt fit with Lin. background')
            ax[0].legend()
            plt.tight_layout()
            if plot_save==True:
                plt.savefig(save_path+'Peak_fit_HCN_ana_'+version+str(name)+'.png',dpi=300)
        return Peaks_pos_Fit,Peaks_amp_Fit
    
def peak_fitting_Voigt_1center(Signal,peaks_CH,peak_properties,fit_range_faktor=0.5,plot=False,plot_save=False):
    Peaks_pos_Fit           = np.zeros((len(peaks_CH)))
    Peaks_amp_Fit           = np.zeros((len(peaks_CH)))
             
    test = np.zeros((len(peaks_CH)))
    
    if plot == True:
        fig, ax = plt.subplots(2, 1,sharex=True,figsize=(15,11)) 
        ax[0].plot(Signal,'k:',linewidth=3,label='Transmission spectra')
        ax[0].plot(peaks_CH,Signal[peaks_CH],'bo',linewidth=3,label='Simple peak detection')
        ax[0].set_title('Voigt one centre for Gauss and Lorentz')

    for j,peak in enumerate(peaks_CH):
        FWHM              = np.round(peak_properties['right_ips'][j])-np.round(peak_properties['left_ips'][j])
        
        pos1            = int(np.round(peak_properties['left_ips'][j]-(FWHM*fit_range_faktor)))
        pos2            = int(np.round(peak_properties['right_ips'][j]+(FWHM*fit_range_faktor)))

        Pos_range       = np.arange(pos1,pos2)
        
        guess_voigtlin                         = [ peak,peak_properties['peak_heights'][j],(FWHM)/2,peak,peak_properties['peak_heights'][j],(FWHM)/2,0,0]   # Center Gaus & Lorentz gleich
        
        popt_voigtlin_2, pcov_voigtlin_2       = curve_fit(voigtlin_same_cen, xdata=Pos_range,ydata=-Signal[pos1:pos2]+1,p0=guess_voigtlin[1:],bounds=(0,np.inf),maxfev=400000)           

        y_voigtlin= -voigtlin_same_cen(Pos_range, *popt_voigtlin_2)+1
                
        presicion = 0.1 # Subsampling des Peaks
        Pos_range_fine  = np.arange(pos1,pos2,presicion)

        y_fine  = -voigtlin_same_cen(Pos_range_fine, *popt_voigtlin_2)+1
        y_peak_voigtlin = -voigtlin_same_cen(popt_voigtlin_2[2], *popt_voigtlin_2)+1
        
        Peaks_pos_Fit[j] = popt_voigtlin_2[2]
        Peaks_amp_Fit[j] = y_peak_voigtlin
                
        if plot == True:
            ax[0].plot(Pos_range_fine,y_fine,'g:')#,label='voigtlin fit')
            ax[1].plot(Pos_range,y_voigtlin-Signal[pos1:pos2],'g:')#,label='voigtlin')
    if plot == True:
        ax[0].plot(Peaks_pos_Fit,Peaks_amp_Fit,'ro',label='Peaks of Voigt+Lin.')
        ax[1].plot(Pos_range,y_voigtlin-Signal[pos1:pos2],'g:',label='Voigt fit with Lin. background')
        
        ax[0].legend()
        plt.tight_layout()
        if plot_save==True:
            plt.savefig(save_path+'Peak_fit_HCN_ana_'+version+str(name)+'.png',dpi=300)
    return Peaks_pos_Fit,Peaks_amp_Fit

def Gascell_peaks_check(Transmission,df_gascell_ref,C2H2 =False,HCN=False,CO12=False,CO13=False,plot=False,save=False,plot_save=False):
        """ Ermittlung des größten Wellenlängen Sprungs = R0 
        Bitte beachten: Es muss der Übergang !R0 zu P1 vorhanden sein!
        R-Branch Parms Anzahl um 1 kleiner, da Array Index bei 0 beginnt.
        Für C2H2 wird der Peak P26 von der Fit Routine nicht erkannt / später ignoriert    """
        if C2H2 ==True:
            gascell_parms = {'height': 0.101,'distance':10000,'R':28,'P':26,'name' : 'C2H2'} 
            
        if HCN ==True:
            #gascell_parms = {'height': 0.05,'distance':5000,'R':26,'P':27,'name' : 'HCN'}
            gascell_parms = {'height': 0.015,'distance':5000,'R':26,'P':27,'name' : 'HCN'}
            
        peaks_Transmission, properties_Transmission = gascell_peak_detection(Transmission,height = gascell_parms['height'],distance=gascell_parms['distance'],plot = plot)

        print ('Peaks detected:',str(len(peaks_Transmission)))
        x1 = np.argmax(np.diff(peaks_Transmission))      
        x2 = len(peaks_Transmission[x1+1:])
        print ("Peaks R-Branch:",x1+1," | P-Branch:" , x2)

        start   = df_gascell_ref.iloc[0,0]
        stop    = df_gascell_ref.iloc[-1,0]
        
        if x1==gascell_parms['R']:
            print ('R branch fits')
        elif x1<gascell_parms['R']:
            print ('R branch peak missing: ',str(gascell_parms['R']-x1))
            start   = df_gascell_ref.iloc[gascell_parms['R']-x1,0]
        else:
            print ('R branch peak to much: '+str(x1-gascell_parms['R']))
            peaks_Transmission      = peaks_Transmission[x1-gascell_parms['R']:]
            for key in properties_Transmission:
                properties_Transmission[key] = properties_Transmission[key][x1-gascell_parms['R']:]
            x1 -= x1-gascell_parms['R']
            
        if x2==gascell_parms['P']:
            print ('P branch fits')

        elif x2<gascell_parms['P']:
            print ('P branch peak missing: ',str(gascell_parms['P']-x2))
            stop    = df_gascell_ref.iloc[(x2-gascell_parms['P']),0]

        else:
            peak_count = x2-gascell_parms['P']
            print ('P branch peak to much:',str(peak_count))
            peaks_Transmission      = peaks_Transmission[:-(peak_count)]
            for key in properties_Transmission:
                properties_Transmission[key] = properties_Transmission[key][:-(peak_count)]
            x2 -=peak_count
            
        print ('WL cal. range: ',start,'nm -',stop,'nm')

        if plot==True:                        
            plt.figure(figsize=(15,11))
            plt.title('Peaks detected:'+str(len(peaks_Transmission))+'\n Peaks used: '+str(x1+1)+' R-branch | '+str(x2)+' P-branch')
            plt.plot(Transmission,'r--',label='Gascell transmittance')
            plt.plot(peaks_Transmission,Transmission[peaks_Transmission],'ko',label=str(len(peaks_Transmission))+' Peaks used')
            plt.plot(peaks_Transmission[x1],Transmission[peaks_Transmission[x1]],'bx',label='Pos. R0')
            plt.plot(peaks_Transmission[0],Transmission[peaks_Transmission[0]],'gs',label='R Branch left')
            plt.plot(peaks_Transmission[-1],Transmission[peaks_Transmission[-1]],'gs',label='P Branch right')
            plt.ylabel('Amplitude arb. unit')
            plt.xlabel('Samples #')
            plt.legend()   
            plt.tight_layout() 
            if plot_save == True:
                plt.savefig(save_path+'Peak_check_'+gascell_parms['name']+version+str(name)+'.png',dpi=300)
            
        return peaks_Transmission, properties_Transmission,(x1,x2)

def gascell_peak_detection(Transmission,width       = 10,height      = 0.05,distance    = 1000, plot =False):
    """ Input: Transmission, Peak parms """
    
    peaks, properties                       = find_peaks((Transmission*-1)+1, height=height, distance=distance, width=width)   
    ####### Amplitude of Signal is inverted so that the valleys in transmission become peaks
    if plot == True:    
        plt.figure()
        plt.plot(Transmission,'r--',label='Gascell Transmittance')
        plt.plot(peaks,Transmission[peaks],'x',label=str(len(peaks))+' Peaks detected')
        plt.ylabel('Amplitude arb. unit')
        plt.ylabel('Amplitude arb. unit')
        plt.xlabel('Samples #')
        plt.legend()   
        plt.tight_layout()

    return peaks, properties

def baseline_cor_2steps(Amp,scale=0.95,fs=200000,lowcut_1st = 5,lowcut_2nd = 15,order=7):
    """
    2 stufiger Filter zur Erkennung der Transmission ohne RR_peaks oder Gasabsorptionen
    Ablauf:
    I.  Basislinie 1 = Tiefpass mit 5 Hz zur Erkennung der Baseline
    II.  Ausschneiden der Transmissions Täler aus Datensatz
    III. Basislinie 2 = Tiefpass mit 15 Hz aus reduzierten Datensatz (ohne Transmissionstäler)
    IV. Rohsignal wird durch Baselinie 2 geteilt
    
    """
    Amp_baseline                           = butter_lowpass_filter(Amp,lowcut_1st,fs,order=order)
    Amp_baseline_red                       = np.where(Amp>Amp_baseline,Amp,Amp_baseline)
    #Amp_baseline_red[Amp>scale*Amp_baseline]     = Amp_baseline[Amp>scale*Amp_baseline]
    Amp_baseline2                          = butter_lowpass_filter(Amp_baseline_red,lowcut_2nd,fs,order=order)
    return (Amp/Amp_baseline2)

def C2H2_gascell_peaks():
    dic_C2H2 ={"R29": [1511.73022, 0.0003],
    "R28": [1512.08822, 0.0003],
    "R27": [1512.45256, 0.00012],
    "R26": [1512.82303, 0.0003],
    "R25": [1513.19984, 0.0003],
    "R24": [1513.58304, 0.0003],
    "R23": [1513.97244, 0.0003],
    "R22": [1514.36815, 0.0003],
    "R21": [1514.77015, 0.0003],
    "R20": [1515.17846, 0.0003],
    "R19": [1515.59306, 0.0003],
    "R18": [1516.01397, 0.0003],
    "R17": [1516.44117, 0.00011],
    "R16": [1516.87458, 0.0003],
    "R15": [1517.31438, 0.0003],
    "R14": [1517.76049, 0.0003],
    "R13": [1518.213, 0.0003],
    "R12": [1518.6717, 0.0003],
    "R11": [1519.13677, 0.00011],
    "R10": [1519.60821, 0.0003],
    "R9": [1520.08592, 0.0003],
    "R8": [1520.56992, 0.0003],
    "R7": [1521.06033, 0.0001],
    "R6": [1521.55714, 0.0003],
    "R5": [1522.06024, 0.0003],
    "R4": [1522.56965, 0.0003],
    "R3": [1523.08546, 0.0003],
    "R2": [1523.60766, 0.0003],
    "R1": [1524.13606, 0.0001],
    "P1": [1525.75984, 0.0006],
    "P2": [1526.31394, 0.0003],
    "P3": [1526.87428, 0.0001],
    "P4": [1527.44107, 0.0001],
    "P5": [1528.01425, 0.0001],
    "P6": [1528.59383, 0.0001],
    "P7": [1529.17983, 0.0003],
    "P8": [1529.77223, 0.0003],
    "P9": [1530.37103, 0.0003],
    "P10": [1530.9762, 0.0001],
    "P11": [1531.58783, 0.0003],
    "P12": [1532.20594, 0.0003],
    "P13": [1532.83039, 0.0001],
    "P14": [1533.46129, 0.0001],
    "P15": [1534.09863, 0.0003],
    "P16": [1534.74243, 0.0003],
    "P17": [1535.39273, 0.0003],
    "P18": [1536.04943, 0.0006],
    "P19": [1536.71253, 0.0003],
    "P20": [1537.38212, 0.0003],
    "P21": [1538.05822, 0.0003],
    "P22": [1538.74082, 0.0003],
    "P23": [1539.42983, 0.00011],
    "P24": [1540.12534, 0.00011],
    "P25": [1540.82734, 0.00011],
    "P26": [1541.53579, 0.0003],
    "P27": [1542.25068, 0.0003],}
    return pd.DataFrame.from_dict(dic_C2H2,orient='index',columns=[ "wavelength/nm", "expanded uncertainties/nm"])

def HCN_gascell_peaks():
    dic_HCN =    {"R26": [1527.63342, 0.00012],
    "R25": [1528.05474, 0.00015],
    "R24": [1528.48574,	9e-05],
    "R23": [1528.92643,	6e-05],
    "R22": [1529.37681,	7e-05],
    "R21": [1529.83688,	6e-05],
    "R20": [1530.30666, 8e-05],
    "R19": [1530.78615, 8e-05],
    "R18": [1531.27537, 7e-05],
    "R17": [1531.7743, 8e-05],
    "R16": [1532.28298, 8e-05],
    "R15": [1532.80139, 7e-05],
    "R14": [1533.32954, 8e-05],
    "R13": [1533.86745, 7e-05],
    "R12": [1534.41514, 6e-05],
    "R11": [1534.97258, 6e-05],
    "R10": [1535.53981, 5e-05],
    "R09": [1536.11683, 4e-05],
    "R08": [1536.70364, 5e-05],
    "R07": [1537.30029, 6e-05],
    "R06": [1537.90675, 0.00013],
    "R05": [1538.52305, 7e-05],
    "R04": [1539.14921, 0.00012],
    "R03": [1539.78523, 9e-05],
    "R02": [1540.4312, 0.0001],
    "R01": [1541.08703, 0.0001],
    "R00": [1541.7528, 6e-05],
    "P01": [1543.11423, 5e-05],
    "P02": [1543.80967, 0.00018],
    "P03": [1544.51503, 8e-05],
    "P04": [1545.23033, 7e-05],
    "P05": [1545.95549, 7e-05],
    "P06": [1546.69055, 8e-05],
    "P07": [1547.43558, 0.00024],
    "P08": [1548.19057, 7e-05],
    "P09": [1548.95555, 4e-05],
    "P10": [1549.73051, 4e-05],
    "P11": [1550.51546, 5e-05],
    "P12": [1551.31045, 9e-05],
    "P13": [1552.11546, 0.0001],
    "P14": [1552.93051, 9e-05],
    "P15": [1553.75562, 0.00012],
    "P16": [1554.59079, 0.0001],
    "P17": [1555.43605, 0.00011],
    "P18": [1556.29141, 0.00015],
    "P19": [1557.15686, 0.00015],
    "P20": [1558.0324, 0.00015],
    "P21": [1558.91808, 0.00014],
    "P22": [1559.81389, 0.00014],
    "P23": [1560.71983, 0.0001],
    "P24": [1561.63593, 9e-05],
    "P25": [1562.56218, 0.00013],
    "P26": [1563.49859, 0.00016],
    "P27": [1564.44519,	2.1E-4],}
    return pd.DataFrame.from_dict(dic_HCN,orient='index',columns=[ "wavelength/nm", "expanded uncertainties/nm"])



def CO12_gascell_peaks():
    dic_CO12 = {"R21":  [1560.5025, 0.006],
    "R20":  [1560.868, 0.007],
    "R19":  [1561.26, 0.007],
    "R18":  [1561.6786, 0.007],
    "R17":  [1562.1237, 0.007],
    "R16":  [1562.5953, 0.005],
    "R15":  [1563.0935, 0.005],
    "R14":  [1563.6183, 0.005],
    "R13":  [1564.1697, 0.005],
    "R12":  [1564.7477, 0.004],
    "R11":  [1565.3523, 0.004],
    "R10":  [1565.9835, 0.004],
    "R9":  [1566.6414, 0.004],
    "R8":  [1567.3261, 0.004],
    "R7":  [1568.0375, 0.004],
    "R6":  [1568.7756, 0.005],
    "R5":  [1569.5405, 0.004],
    "R4":  [1570.3323, 0.004],
    "R3":  [1571.1509, 0.005],
    "R2":  [1571.9965, 0.005],
    "R1":  [1572.8691, 0.004],
    "R0":  [1573.7687, 0.004],
    "P1":  [1575.6498, 0.004],
    "P2":  [1576.6311, 0.004],
    "P3":  [1577.6397, 0.006],
    "P4":  [1578.6758, 0.004],
    "P5":  [1579.7392, 0.005],
    "P6":  [1580.83, 0.005],
    "P7":  [1581.9485, 0.004],
    "P8":  [1583.0945, 0.005],
    "P9":  [1584.2683, 0.005],
    "P10":  [1585.4698, 0.005],
    "P11":  [1586.6993, 0.005],
    "P12":  [1587.9567, 0.006],
    "P13":  [1589.2422, 0.006],
    "P14":  [1590.5559, 0.006],
    "P15":  [1591.8978, 0.005],
    "P16":  [1593.2681, 0.006],
    "P17":  [1594.6669, 0.006],
    "P18":  [1596.0942, 0.007],
    "P19":  [1597.5502, 0.006],}
    return pd.DataFrame.from_dict(dic_CO12,orient='index',columns=[ "wavelength/nm", "expanded uncertainties/nm"])

def CO13_gascell_peaks():
    dic_CO13 = {"R21": [1595.37720, 6.0E-03],
    "R20": [1595.75540, 7.0E-03],
    "R19": [1596.15950, 7.0E-03],
    "R18": [1596.58950, 6.0E-03],
    "R17": [1597.04540, 7.0E-03],
    "R16": [1597.52710, 7.0E-03],
    "R15": [1598.03490, 5.0E-03],
    "R14": [1598.56860, 5.0E-03],
    "R13": [1599.12840, 5.0E-03],
    "R12": [1599.71410, 5.0E-03],
    "R11": [1600.32580, 5.0E-03],
    "R10": [1600.96360, 4.0E-03],
    "R9": [1601.62740, 5.0E-03],
    "R8": [1602.31740, 5.0E-03],
    "R7": [1603.03340, 4.0E-03],
    "R6": [1603.77560, 5.0E-03],
    "R5": [1604.54390, 5.0E-03],
    "R4": [1605.33850, 4.0E-03],
    "R3": [1606.15930, 6.0E-03],
    "R2": [1607.00640, 6.0E-03],
    "R1": [1607.87990, 4.0E-03],
    "R0": [1608.77990, 4.0E-03],
    "P1": [1610.65960, 4.0E-03],
    "P2": [1611.63930, 4.0E-03],
    "P3": [1612.64570, 6.0E-03],
    "P4": [1613.67880, 4.0E-03],
    "P5": [1614.73880, 6.0E-03],
    "P6": [1615.82550, 6.0E-03],
    "P7": [1616.93920, 5.0E-03],
    "P8": [1618.07970, 5.0E-03],
    "P9": [1619.24730, 5.0E-03],
    "P10": [1620.44200, 6.0E-03],
    "P11": [1621.66410, 5.0E-03],
    "P12": [1622.91320, 6.0E-03],
    "P13": [1624.18980, 5.0E-03],
    "P14": [1625.49400, 6.0E-03],
    "P15": [1626.82570, 7.0E-03],
    "P16": [1628.18510, 6.0E-03],
    "P17": [1629.57230, 7.0E-03],
    "P18": [1630.98730, 7.0E-03],
    "P19": [1632.43030, 6.0E-03],}
    return pd.DataFrame.from_dict(dic_CO13,orient='index',columns=[ "wavelength/nm", "expanded uncertainties/nm"])
