# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:15:01 2020

@author: eiserm01
"""
import numpy as np
from scipy.optimize import curve_fit

def linear(x,a,b,cen):
    return a + b*(x-cen)

def poly(x,a,b,c,cen):
    return a + b*(x-cen) + c*(x-cen)**2

def gauss(x, cen, amp, sigma):
    return amp*np.exp(-(x-cen)**2/(2*(sigma**2)))

def gaussconst(x, cen, amp, sigma, a):
    return amp*np.exp(-(x-cen)**2/(2*(sigma**2))) + a

def gausslin(x, cen, amp, sigma, a, b):
    return amp*np.exp(-(x-cen)**2/(2*(sigma**2))) + a + b*(x-cen)

def gausspoly(x, cen, amp, sigma, a, b, c):
    return amp*np.exp(-(x-cen)**2/(2*(sigma**2))) + a + b*(x-cen) + c*(x-cen)**2

def lorentzian(x, cen,amp, wid):
    return (amp*wid**2/((x-cen)**2+wid**2))

def lorentzianconst(x, cen,amp, wid, a):
    return (amp*wid**2/((x-cen)**2+wid**2)) + a

def lorentzianlin(x, cen, amp, wid, a, b):
    return (amp*wid**2/((x-cen)**2+wid**2)) + a + b*(x-cen)

def lorentzianpoly(x, cen, amp, wid, a, b, c):
    return (amp*wid**2/((x-cen)**2+wid**2)) + a + b*(x-cen) + c*(x-cen)**2

def voigt(x, cenG1, ampG1, sigmaG1, cenL1, ampL1, widL1):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) )
    
def voigtconst(x, cenG1, ampG1, sigmaG1, cenL1, ampL1, widL1, a):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) ) + a
    
def voigtlin(x, cenG1, ampG1, sigmaG1, cenL1, ampL1, widL1, a, b):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) ) + a + b*(x-cenL1)
    
def voigtlin_same_cen(x, ampG1, sigmaG1, cenL1, ampL1, widL1, a, b):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenL1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) ) + a + b*(x-cenL1)

def voigtlin_same_cen_const(x, ampG1, sigmaG1, cenL1, ampL1, widL1):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenL1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)))
    
def voigtpoly(x, cenG1, ampG1, sigmaG1, cenL1, ampL1, widL1, a, b, c):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/((2*sigmaG1)**2)))) +\
    ((ampL1*widL1**2/((x-cenL1)**2+widL1**2)) ) + a + b*(x-cenL1) + c*(x-cenL1)**2

def _3lorentzianlin(x, cen1, amp1, wid1, cen2, amp2, wid2, cen3, amp3, wid3, a, b, bgcen):
    return (-(amp1*wid1**2/((x-cen1)**2+wid1**2))+1) * (-(amp2*wid2**2/((x-cen2)**2+wid2**2))+1) * (-(amp3*wid3**2/((x-cen3)**2+wid3**2))+1) + a + b*(x-bgcen)

def _3lorentzianpoly(x, cen1, amp1, wid1, cen2, amp2, wid2, cen3, amp3, wid3, a, b, c, bgcen):
    return (-(amp1*wid1**2/((x-cen1)**2+wid1**2))+1) * (-(amp2*wid2**2/((x-cen2)**2+wid2**2))+1) * (-(amp3*wid3**2/((x-cen3)**2+wid3**2))+1) + a + b*(x-bgcen) + c*(x-bgcen)**2

def WL_fit_2order(x, x0, B1,B2):
    return (x0 +  B1 * x + B2 * x**2)  
      
def WL_fit_3order(x, x0, B1,B2,B3):
    return (x0 +  B1 * x + B2 * x**2 + B3 * x**3  )

def WL_fit_4order(x, x0, B1,B2,B3,B4):
    return (x0 +  B1 * x + B2 * x**2 + B3 * x**3 + B4 * x**4 )

def WL_fit_5order(x, x0, B1,B2,B3,B4,B5):
    return (x0 +  B1 * x + B2 * x**2 + B3 * x**3 + B4 * x**4 + B5 * x**5  )

def WL_fit_6order(x, x0, B1,B2,B3,B4,B5,B6):
    return (x0 +  B1 * x + B2 * x**2 + B3 * x**3 + B4 * x**4 + B5 * x**5 +B6 * x**6 )

def WL_fit_7order(x, x0, B1,B2,B3,B4,B5,B6,B7):
    return (x0 +  B1 * x + B2 * x**2 + B3 * x**3 + B4 * x**4 + B5 * x**5 +B6 * x**6 + B7 * x**7 )

def SFBG_fit(x,x_0,y_0,a,w,w2):
    return (y_0+a*np.exp(-np.exp((x-x_0)/w)+((x-x_0)/w2)))

def get_SFBG_fit(xdata,ydata,fit='SFBG_fit'):
    x_max       = xdata[np.argmax(ydata)]
    y_min       = np.amin(ydata)
    
    p0=[x_max,0.3,1,1,1]
    popt_SFBG, pcov_SFBG            = curve_fit(SFBG_fit, xdata=xdata,ydata=ydata,p0=p0)
    y_SFBG                          = SFBG_fit(xdata, *popt_SFBG)
    perr                            = np.sqrt(np.diag(pcov_SFBG))
    return popt_SFBG, pcov_SFBG,y_SFBG,perr 

def SFBG_fit_function(x,x_0,y_0,a,w,w2):    #no background
    return (y_0+a*np.exp(-np.exp((x-x_0)/w)+((x-x_0)/w2)))

def SFBG_fit_function_lin(x,x_0,y_0,a,w,w2,m):  #linear background
    return (y_0+a*np.exp(-np.exp((x-x_0)/w)+((x-x_0)/w2)))+m*(x-x_0)

def SFBG_fit_function_poly(x,x_0,y_0,a,w,w2,m,n):    #polynomial background
    return (y_0+a*np.exp(-np.exp((x-x_0)/w)+((x-x_0)/w2)))+m*(x-x_0)+n*(x-x_0)**2

def SFBG_fit(WL,AMP, bounds=([1500,0,0,0,0],[1600,1,2,10,10]) ):
    x_max=WL[np.argmax(AMP)]
    y_min=np.amin(AMP)
    param,param_cov=curve_fit(SFBG_fit_function,WL,AMP,p0=[x_max,y_min,1,1,1],bounds=bounds)
    y_fit                           = SFBG_fit_function(xdata, *param)
    perr                            = np.sqrt(np.diag(param_cov))
    return param,param_cov,y_fit,perr

def SFBG_fit_lin(WL,AMP, bounds=([1500,0,0,0,0,-1],[1600,1,2,10,10,1])):
    x_max=WL[np.argmax(AMP)]
    y_min=np.amin(AMP)
    param,param_cov=curve_fit(SFBG_fit_function_lin,WL,AMP,p0=[x_max,y_min,1,1,1,0.1],bounds=bounds)
    y_fit                           = SFBG_fit_function_lin(xdata, *param)
    perr                            = np.sqrt(np.diag(param_cov))
    return param,param_cov,y_fit,perr

def SFBG_fit_poly(WL,AMP, bounds=([1500,0,0,0,0,-1,-1],[1600,1,2,10,10,1,1])):
    x_max=WL[np.argmax(AMP)]
    y_min=np.amin(AMP)
    param,param_cov=curve_fit(SFBG_fit_function_poly,WL,AMP,p0=[x_max,y_min,1,1,1,0.1,0],bounds=bounds)
    y_fit                           = SFBG_fit_function_poly(xdata, *param)
    perr                            = np.sqrt(np.diag(param_cov))
    return param,param_cov,y_fit,perr