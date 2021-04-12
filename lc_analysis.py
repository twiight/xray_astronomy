#!/bin/bash
# -*- coding: utf-8 -*-
# written by Tong
# plot the phase-folded light curve from txt file (version for xmm)
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
#import correct as correct
from datetime import datetime
from scipy.interpolate import lagrange
from scipy import interpolate
from scipy.optimize import curve_fit
import pandas as pd
from astropy.timeseries import LombScargle
import stingray
from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum

font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 16, }

def filter_energy(time,energy,band):
    T=time
    E=energy
    i=0
    if len(time)!=len(energy):
        print('error')
        return None
    else:
        while i <len(E):
            if E[i]<band[0] or E[i]>band[1]:
                E=np.delete(E,i)
                T=np.delete(T,i)
            else:
                i+=1
    return T
def get_LS(time, flux,freq):
    x = time
    y = flux

    # LS = LombScargle(x, y, dy = 1, normalization = 'standard', fit_mean = True,
    #                  center_data = True).power(freq, method = 'cython')
    LS = LombScargle(x, y,normalization = 'standard')
    power = LS.power(freq)
    FP=LS.false_alarm_probability(power.max(),minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_99 = LS.false_alarm_level(0.0027, minimum_frequency = freq[0], maximum_frequency = freq[-1],method='baluev')
    FP_90 = LS.false_alarm_level(0.05,  minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    FP_68 = LS.false_alarm_level(0.32, minimum_frequency=freq[0],
                                 maximum_frequency=freq[-1], method='baluev')
    plt.plot([freq[0], freq[-1]], [FP_99, FP_99], '--')
    plt.plot([freq[0], freq[-1]], [FP_90, FP_90], '--')
    plt.plot([freq[0], freq[-1]], [FP_68, FP_68], '--')

    plt.title('FP={0}'.format(FP))
    plt.semilogx()
    # plt.xlim(1000.,1500.)
    plt.plot(freq, power)
    print(1./freq[np.where(power==np.max(power))])

    plt.show()
    res=1e5*power
    res=np.round(res,2)
    return [FP, 1. / freq[np.where(power == np.max(power))],np.max(power),res]

def plot_pds(time,flux):
    lc = Lightcurve(time, flux)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(lc.time, lc.counts, lw=2, color='blue')
    ax.set_xlabel("Time (s)", fontproperties=font1)
    ax.set_ylabel("Counts (cts)", fontproperties=font1)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.tick_params(which='major', width=1.5, length=7)
    ax.tick_params(which='minor', width=1.5, length=4)
    plt.show()

    ps = Powerspectrum(lc,norm='leahy')
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6), sharex=True)
    ax1.loglog()
    ax1.step(ps.freq, ps.power, lw=2, color='blue')
    ax1.set_ylabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()

    avg_ps = AveragedPowerspectrum(lc, 500,dt=lc.time[1]-lc.time[0],norm='leahy')
    print("Number of segments: %d" % avg_ps.m)
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
    ax1.loglog()
    ax1.step(avg_ps.freq, avg_ps.power, lw=2, color='blue')
    ax1.set_xlabel("Frequency (Hz)", fontproperties=font1)
    ax1.set_ylabel("Power (raw)", fontproperties=font1)
    ax1.set_yscale('log')
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.tick_params(which='major', width=1.5, length=7)
    ax1.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    plt.show()


def read_SAS_lc():
    # 1,2,3分别代表mos1,mos2,pn的light curve，也可以加起来用，记为_all;
    # 根据实际情况来决定lomb-scargle的输入
    dt = 1
    path='/Users/baotong/xmm/0302900101/cal/'
    os.chdir(path)
    mode=['mos1','mos2','pn']
    filename1=mode[0]+'_lccorr_bin{0}.lc'.format(int(dt));filename2=mode[1]+'_lccorr_bin{0}.lc'.format(int(dt));filename3=mode[2]+'_lccorr_bin{0}.lc'.format(int(dt))
    # lc1=fits.open(filename1);
    # lc2=fits.open(filename2);
    lc3=fits.open(filename3)
    # time1=lc1[1].data['TIME'][0:-100];rate1=lc1[1].data['RATE'][0:-100]
    # time2=lc2[1].data['TIME'][0:-100];rate2=lc2[1].data['RATE'][0:-100]
    time3=lc3[1].data['TIME'][20000:-50000];rate3=lc3[1].data['RATE'][20000:-50000]
    rate3=np.nan_to_num(rate3)
    rate3[np.where(rate3<0)]=0
    freq=np.arange(1./10000,0.5/dt,1./(10*100000))

    time_all=time3;rate_all=rate3
    get_LS(time_all,rate_all,freq)
    # plot_pds(time_all,rate_all)
read_SAS_lc()

