#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:58:58 2020

@author: jimmy
"""
import numpy as np
import matplotlib.pyplot as plt
from MathTools import Gaussian
from scipy.optimize import curve_fit

class SpectralFeatures():
    
    # Initialize input parameters
    def __init__(self,wavelength,offset,P,mean,sigma,sampleInterval,SNR):
        # Inputs:
        #   P: peak depth of spectrum
        #   sigma: width of spectrum (in Angstroms)
        #   mean: center of spectrum (in Angstroms)
        #   sampleInterval: spacing of pixels (in Angstroms)
        #   SNR: signal-to-noise ratio you want to generate
        
        self.SNR = SNR
        self.wavelength = wavelength
        self.offset = offset
        self.P = P
        self.mean = mean
        self.sigma = sigma
        self.sample = sampleInterval
    
    def GaussianNoise(self,plot=True):
        # Returns:
        #   Plot of Gaussian spectrum with user-defined properties
        
        # Generate list of wavelengths to draw spectrum over
        x = self.wavelength+np.arange(self.mean-5.0,self.mean+5.0+self.sample,self.sample)

        # Calculate value of Gaussian function at each x
        function = Gaussian(x,self.offset,self.P,self.mean,self.sigma)

        # Calculate photometric precision from SNR (SNR=1/sigma_prec)
        sigma_prec = 1.0/self.SNR
        
        # Generate Gaussian noise depending on SNR
        noise = np.random.normal(0,sigma_prec,len(x))
        
        # Add Gaussian noise to Gaussian function values
        noiseData = np.add(function,noise)

        if plot == True:
            # Plot noisy data
            plt.scatter(x,noiseData)
            plt.xlabel(r'Wavelength ($\AA$)',fontsize=14)
            plt.title('Gaussian Spectrum with SNR={0}:1'.format(self.SNR),fontsize=18)
            plt.ylim(0,1.4)
            plt.tight_layout()
    
        return(x,noiseData)
    
    # Function to fit a Gaussian model to Gaussian with noise
    def GaussianModel(self,init_val,plot=False):
        # Returns: best-fit parameters to model noisy Gaussian
        
        # Extract wavelength and Gaussian values from GaussianNoise
        x, y = self.GaussianNoise(plot=False)
        
        # Use curve_fit to find best parameters to fit model
        guesses = init_val  # for [P, mean, sigma]
        best_vals, covar = curve_fit(Gaussian, x, y, p0=guesses)
        #print('best_vals: {}'.format(best_vals))
        
        # Calculate y values of model
        y_model = Gaussian(x,best_vals[0],best_vals[1],best_vals[2],best_vals[3])
        
        if plot == True:
            # Plot data vs. model
            plt.scatter(x,y,label='Data',color='red')
            plt.plot(x,y_model,label='Model',color='blue')
            plt.legend()
            plt.xlabel(r'Wavelength ($\AA$)',fontsize=14)
            plt.title('Gaussian Spectrum with SNR={0}:1'.format(self.SNR),fontsize=18)
            plt.ylim(0,1.4)
            plt.tight_layout()
        
        mean_model = best_vals[2]
        perc_diff = np.abs(mean_model)/self.wavelength*100
        
        return(self.wavelength+mean_model)
    