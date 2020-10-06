#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:22:55 2020

@author: jimmy
"""
# Import relevant modules/packages
import numpy as np
import matplotlib.pyplot as plt

# Class to calculate Fourier components
class Fourier:
    
    # Initialize the instance of this Class with the following properties
    def __init__(self,c,start,end,period=1.0,t0=0.0,parity='both'):
        # Inputs
        #    c: list of coefficients in front of Fourier components
        #    start: multiple of period to set as beginning of time list
        #    end: multiple of period to set as end of time list
        #    period: period of signal (default is 1)
        #    t0: reference time for signal (default is 0)
        #    parity: parameter to choose to plot even, odd, or even&odd ('both') terms
        
        # Define list of coefficients
        self.c = c
        
        # Define parity
        self.parity = parity
        
        # Calculate renormalization factor
        self.R = np.sum(c)
        
        # Define period and reference time
        self.period = period
        self.t0 = t0
        
        # Make a list of times to evaluate signal at (2 periods)
        self.t = np.linspace(start*period,end*period,100)
        
        # Calculate phase of Fourier components
        self.phase = [(t-self.t0)/self.period for t in self.t]

    # Function to calcuate signal vs. time
    def FourierSignal(self):
        
        # Define empty list of signal amplitude values
        S = []
        
        # Define multiple within cosine argument depending on which terms you're interested in
        if self.parity == 'even':
            
            # Makes evenly spaced list of even numbers
            multiples = np.arange(2,2*len(self.c),2)      
                
        elif self.parity == 'odd':
            
            # Makes evenly spaced list of odd numbers
            multiples = np.arange(1,2*len(self.c)+1,2)
            
        elif self.parity == 'both':
            
            # Make list of numbers from 0 to length of coefficient list
            multiples = np.arange(1,len(self.c),1)
            
        # Calculate signal value for each phase value
        for phi in self.phase:
            
            # Define C_0 as first term (assumed to be first term of user-entered coefficient array)
            terms = self.c[0]
            
            # Add all Fourier terms you're interested in
            for i in range(1,len(self.c)):
                terms += self.c[i]*np.cos(multiples[i-1]*2*np.pi*phi)
                
            # Append normalized signal value to signal array
            S.append(1/self.R*terms)

        return(self.c,self.phase,self.t,S)
    
    # Function to plot signal vs. time
    def FourierPlotter(self,xaxis,linestyle='solid',legend=True):
        
        c,phi,t,S = self.FourierSignal()
        
        # Plot signal vs. time
        ax = plt.subplot(111)
        
        # Plot S vs. t
        label = r'First {0} {1} coeff.'.format(len(c),self.parity)
        
        if xaxis == 'phase':
            ax.plot(phi,S,label=label,linestyle=linestyle)
            ax.set_xlabel('Phase',fontsize=14)
        elif xaxis == 'time':
            ax.plot(t,S,label=label,linestyle=linestyle)
            ax.set_xlabel('Time (period)',fontsize=14)

        # Add plot labels
        ax.axhline(0,color='black',linestyle='dashed')
        
        ax.set_ylabel('S(t)',fontsize=14)
        
        if legend == True:
            ax.legend(loc='center left',bbox_to_anchor=(1, 0.5),fontsize=12)
            
        ax.set_title('Signal vs. Time',fontsize=18)