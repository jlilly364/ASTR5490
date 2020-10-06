#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:56:37 2020

@author: jimmy
"""
# Import relevant modules/packages
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.timeseries import LombScargle
from ReadFile import Read

class Periodicity:
    
    def __init__(self,filename,objectname,numPoints,period=None):
        # Inputs:
        #     filename: file path or file name (if file in same folder as notebook)
        #     objectname: name of object you're plotting curve of
        #     numPoints: number of data points you want to use from the file
        #     period: period of plot feature (only used if xaxis='Phase' to fold the data)
        self.file = filename
        self.object = objectname
        self.numPoints = numPoints
        self.period = period
        
        # Extract time, flux, and error data from text file
        self.times,self.fluxes,self.errors = Read(self.file)


    # Function for plotting light curves from text file
    def LightCurve(self,plot=True,xaxis='Time',curve='Flux'):
        # Inputs:
        #     plot: boolean to decide whether to plot the data (True) or not (False)
        #     xaxis: decide which x parameter to calculate/plot ('Time or Phase')
        #     curve: string that indicates y parameter being plotted (used in axis label)
        # Returns:
        #     xdata: array with data from x-axis
        #     fluxes: array with associated y-axis data
        #     errors: measurement errors read from text file
        
        # Decide what times array to make
        if self.numPoints == None:
            self.times = [time - self.times[0] for time in self.times]
        else:
            self.times = [time - self.times[0] for time in self.times[:self.numPoints]]
            self.fluxes = self.fluxes[:self.numPoints]       
        
        # Define list of data to plot on x-axis (time or phase)
        xdata = []
        
        if plot == True:
            # Initialize axis figure and axis
            fig = plt.figure()
            ax = fig.add_subplot(111)
    
            # Decide what x-axis should be
            if xaxis == 'Time':
                xlabel = 'Time (days)'
    
                # Plot flux vs time
                ax.scatter(self.times,self.fluxes)
    
                xdata = self.times
    
            elif xaxis == 'Phase':
                xlabel = 'Phase'
    
                # Calculate phase from time data
                phases = [(time%self.period)/self.period for time in self.times]
    
                # Make a scatter plot of flux vs. phase
                ax.scatter(phases,self.fluxes,label='Period = {0:.3f} days'.format(self.period))
    
                xdata = phases
    
            # Add plot features
            ax.set_xlabel(xlabel,fontsize=14)
            ax.set_ylabel('{0}'.format(curve),fontsize=14)
            ax.set_title('{0} vs. {1} for {2}'.format(curve,xaxis,self.object),fontsize=18)
            ax.legend()
        else:
            # Decide what x-axis should be
            if xaxis == 'Time':
                xdata = self.times
            elif xaxis == 'Phase':
                # Calculate phase from time data
                phases = [(time%self.period)/self.period for time in self.times]
                xdata = phases
        
        return(xdata,self.fluxes,self.errors)


    # Function to generate power spectrum from a flux vs. time dataset
    def LS(self,minP,maxP,numIntervals,i,flux,error,plot=False,trueP=None,days=None):
    
        if len(flux) == 0:
            flux = self.fluxes
            error = self.errors
        
        if days != None:
            increment = days*50 # Approximately 50 data points per day
            self.times,self.fluxes,self.errors = self.times[:increment],self.fluxes[:increment],self.errors[:increment]
        
        # Define range of frequencies to search over
        minfreq = 1./maxP
        maxfreq = 1./minP
    
        # Make list of frequencies within the range
        frequency = np.linspace(minfreq,maxfreq,numIntervals)
    
        # Use LombScargle method to calculate power as a function of those frequencies
        power = LombScargle(self.times,flux,error,nterms=i).power(frequency)
        
        # Find maximum power and frequency/period of maximum power
        maxp = np.max(power) 
        maxind = np.argmax(power)
        
        maxfreq = frequency[maxind]
        best_period = 1./maxfreq
        
        if plot == True:
            # Plot power spectrum using lists from above
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(frequency,power)
            
            # Set axes limits
            ax.set(xlim=(frequency[0],frequency[-1]), ylim=(0,np.max(power)))
            ax.set_xlabel('Frequency (1/days)',fontsize=14)
            ax.set_ylabel('Power',fontsize=14)
            ax.set_title('Power vs. Freq. for {0}'.format(self.object),fontsize=18)
            
            # Plot line indicating period of system from SIMBAD
            if trueP != None:
                ax.vlines(1./trueP,0,1,linestyle='dashed',label='Published Period ({0:.3f} days)'.format(trueP),alpha=0.75)
            
            # Plot vertical line of best period
            ax.vlines(1./best_period,0,1,linestyle='dashed',label='Dominant Period ({0:.3f} days)'.format(best_period),color='red',alpha=0.5)
            ax.legend(loc='center left',bbox_to_anchor=(1, 0.5))
        
        return(maxp)

#Virb = Periodicity('Vogt2009_61Virb_vels.dat','61 Vir b',None,None)
#Virb.LS = Virb.LS(0.1,100,1000,1)