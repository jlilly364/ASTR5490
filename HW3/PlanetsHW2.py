#!/usr/bin/env python
# coding: utf-8

# # <u> ASTR 5490: Homework 2 (Time-Domain Astronomy / Fourier Transforms) </u>

# # 1) Experimenting with Fourier Components
# ## $S(t) = C_0 + C_1\cos\left(1\frac{2\pi(t-t_0)}{P}\right) + C_2\cos\left(2\frac{2\pi(t-t_0)}{P}\right) + C_3\cos\left(3\frac{2\pi(t-t_0)}{P}\right) + ...$

# Import relevant modules/packages
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.timeseries import LombScargle

# Function for plotting light curves from text file
def LightCurve(filename,objectname,numPoints,period=None,plot=True,xaxis='Time',curve='Flux'):
    # Inputs:
    #     filename: file path or file name (if file in same folder as notebook)
    #     objectname: name of object you're plotting curve of
    #     numPoints: number of data points you want to use from the file
    #     period: period of plot feature (only used if xaxis='Phase' to fold the data)
    #     plot: boolean to decide whether to plot the data (True) or not (False)
    #     xaxis: decide which x parameter to calculate/plot ('Time or Phase')
    #     curve: string that indicates y parameter being plotted (used in axis label)
    # Returns:
    #     xdata: array with data from x-axis
    #     fluxes: array with associated y-axis data
    #     errors: measurement errors read from text file
    
    # Extract time, flux, and error data from text file
    times,fluxes,errors = np.loadtxt(filename,skiprows=1,unpack=True,usecols=[0,1,2])
    
    # Decide what times array to make
    if numPoints == None:
        times = [time - times[0] for time in times]
    else:
        times = [time - times[0] for time in times[:numPoints]]
        fluxes = fluxes[:numPoints]       
    
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
            ax.plot(times,fluxes)

            xdata = times

        elif xaxis == 'Phase':
            xlabel = 'Phase'

            # Calculate phase from time data
            phases = [(time%period)/period for time in times]

            # Make a scatter plot of flux vs. phase
            ax.scatter(phases,fluxes,label='Period = {0:.3f} days'.format(period))

            xdata = phases

        # Add plot features
        ax.set_xlabel(xlabel,fontsize=14)
        ax.set_ylabel('{0}'.format(curve),fontsize=14)
        ax.set_title('{0} vs. {1} for {2}'.format(curve,xaxis,objectname),fontsize=18)
        ax.legend()
    else:
        # Decide what x-axis should be
        if xaxis == 'Time':
            xdata = times
        elif xaxis == 'Phase':
            # Calculate phase from time data
            phases = [(time%period)/period for time in times]
            xdata = phases
    
    return(xdata,fluxes,errors)

P=258

# Extract time, phase, and velocity lists from LightCurve function at set period
time_list, vel_obs, errors_obs = LightCurve('HD89744_vels.dat','HD89744_vels',numPoints=None,period=P,plot=False,xaxis='Time',curve='Radial Velocity')
phase_obs, vel_obs, errors_obs = LightCurve('HD89744_vels.dat','HD89744_vels',numPoints=None,period=P,plot=False,xaxis='Phase',curve='Radial Velocity')

# Function to compare observed to modeled phased velocity curve
def Model(time_list,phase_obs,vel_obs,errors_obs,t0,P,e,K,w,gamma,plot=True):
    # Inputs:
    #     time_list = list of times (when observations were taken)
    #     phase_obs = phases corresponding to times in time_list (found with xaxis = 'Phase' from LightCurve func.)
    #     vel_obs = list of observed velocities 
    #     t0 = reference time for orbit
    #     P = guessed period of orbit
    #     e = eccentricity of orbit (0-1)
    #     K = velocity semi-amplitude (folds in a,i,e, and P)
    #     w = argument of periapsis
    #     gamma = constant in velocity equation (shifts graph up or down)
    #     plot = boolean telling the function to produce a plot (True) or not (False)
    # Returns:
    #     vel_model = velocities calculated for model
    
    # Redefine w in radians
    w *= np.pi/180
    
    # Calculate first term of E and theta relationship
    term1 = ((1+e**2)/(1-e**2))**(1/2)
    
    # Create empty lists of velocity values and residuals
    vel_model = []
    residuals = []
    
    # Calculate E at different times
    for i in range(len(time_list)):
        
        # Select particular time from list
        time = time_list[i]
        
        # Calculate mean anomaly (M) to use it as first guess for finding E
        M = 2*np.pi/P*(time-t0)
        
        # Define Kepler equation and its derivative
        f = lambda E: E - e*np.sin(E) - M
        df = lambda E: 1. - e*np.cos(E)
        
        # Use NewtonRaphson function to find E at each time
        E = NewtonRaphson(f,df,M,1e-6,1000)
        
        # Convert eccentric anomaly (E) to true anomaly (theta)
        phase = 2*np.arctan(term1*np.tan(E/2.))

        # Calculate velocity with e, theta, K, and w
        velocity = K*(np.cos(phase+w)+e*np.cos(w))+gamma
        vel_model.append(velocity)
        
    if plot == True:
        # Plot observed vs. modeled phase velocity curve
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.scatter(phase_obs,vel_obs,color='blue',label='Observation')
        ax.scatter(phase_obs,vel_model,color='r',marker='s',label='Model')
        ax.set_xlabel('Phase',fontsize=14)
        ax.set_ylabel('Radial Velocity',fontsize=14)
        ax.set_title('Fitting Model to HD89744',fontsize=18)
        ax.legend()
    
    return(vel_model)


# In[250]:


# Set fixed parameters in Kepler's equation
t0 = 107. # in days
P = 258. # in days
e = 0.7

# Set parameters in velocity equation
K = 250
w = 190
gamma = -270

# Plot observations vs. model for best set of parameters
vel_model = Model(time_list,phase_obs,vel_obs,errors_obs,107.0,P,e,K,w,gamma,plot=True)


# ### Best model results: $P=258$ days, $e=0.70$, $\omega=190^{\circ}$, $K=250$ m/s, $t_0=107$ days, $\gamma=-270$ m/s
# ### Published model used: $P=256.78$ days, $e=0.689$, $\omega=194^{\circ}$, $K=263$ m/s, $t_0=107.1$ days
# Wittenmyer et al., 2009 (https://arxiv.org/pdf/0706.1962.pdf)
# ### One difficulty I had, which was of my own creation, is that I was initially plotting the model velocity vs. the true anomaly ($\theta$) instead of vs. the phased time values from the data itself. This caused my model to not line up at all in the beginning.
# ### For estimating K, I estimated what the semi-amplitude was of the raw light curve which was a helpful first guess. For e, I started with 0.5 since it ranges from 0 to 1 and was pretty quickly able to see what a good value would be. I took a similar approach for $t_0$ because I knew it had to be somewhere between 0 and 500 days so I made the model vs. observation graph for 10 $t_0$ values over this interval and was able to visually converge on a good value. For the period, my first guess was that derived from the periodogram so I only had to minorly tweak it to improve the model fit. My $\omega$ value is clearly very different from the published results, but after adjusting it, 10 degrees seemed to fit best. The paper did not list a $\gamma$ value so I chose one that best fit my data. 

# In[61]:


# Using best parameters from paper
vel_published = Model(time_list,phase_obs,vel_obs,errors_obs,107.1,256.78,.689,264,194,gamma=-265,plot=True)


# ## 4d) Compute $\chi^2$ and $\tilde{\chi}^2$ (reduced chi-squared) for the best-fitting model from 4c. Interpolate the model at the phases of your data to calculate the residuals from the best-fitting model.

# In[211]:


# Calculate (reduced) chi squared for HD89744
ChiSq, RedChiSq, degrees = ChiSquared(vel_model,vel_obs, errors_obs,6)
print(ChiSq,RedChiSq,degrees)


# ### According to this table from NIST (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm), if the number of degrees of freedom is 88, then the probability less than the critical value is between 0.10 and 0.05. There may be other factors in the star-planet system itself that lead to such a high $\tilde{\chi}^2$. The provided errors may be too small if the surface of the star isn't uniform (star spots) or if the star is pulsating. Also, if another body is present, then this would affect what the errors should be too. The same team detected a second planet around HD89744 with a longer period (https://arxiv.org/abs/1901.08471).

# ## 4e) Use scipy.optimize.curve_fit or lmfit to obtain best parameters and report $\tilde{\chi}^2$. Use table of $\tilde{\chi}^2$ probabilities to give a probability of obtaining such a $\tilde{\chi}^2$ by chance and discuss if the model is a good fit to the data.
# ### lmfit: https://github.com/lmfit/lmfit-py/blob/master/README.rst

# In[83]:


from lmfit import minimize, Parameters

# Function to compare observed to modeled phased velocity curve
def ModelNew(params,time_list,phase_obs,vel_obs,errors_obs,plot=False):
    # Inputs:
    #     time_list = list of times (when observations were taken)
    #     phase_obs = phases corresponding to times in time_list (found with xaxis = 'Phase' from LightCurve func.)
    #     vel_obs = list of observed velocities 
    #     t0 = reference time for orbit
    #     P = guessed period of orbit
    #     e = eccentricity of orbit (0-1)
    #     K = velocity semi-amplitude (folds in a,i,e, and P)
    #     w = argument of periapsis
    #     gamma = constant in velocity equation (shifts graph up or down)
    #     plot = boolean telling the function to produce a plot (True) or not (False)
    # Returns:
    #     vel_model = velocities calculated for model
    
    # Extract parameters from params object
    t0 = params['t0']
    P = params['P']
    e = params['e']
    K = params['K']
    w = params['w']
    gamma = params['gamma']
    
    # Redefine w in radians
    w *= np.pi/180
    
    # Calculate first term of E and theta relationship
    term1 = ((1+e**2)/(1-e**2))**(1/2)
    
    # Create empty lists of velocity values and residuals
    vel_model = []
    residuals = []
    
    # Calculate E at different times
    for i in range(len(time_list)):
        
        # Select particular time from list
        time = time_list[i]
        
        # Calculate mean anomaly (M) to use it as first guess for finding E
        M = 2*np.pi/P*(time-t0)
        
        # Define Kepler equation and its derivative
        f = lambda E: E - e*np.sin(E) - M
        df = lambda E: 1. - e*np.cos(E)
        
        # Use NewtonRaphson function to find E at each time
        E = NewtonRaphson(f,df,M,1e-6,1000)
        
        # Convert eccentric anomaly (E) to true anomaly (theta)
        phase = 2*np.arctan(term1*np.tan(E/2.))

        # Calculate velocity with e, theta, K, and w
        velocity = K*(np.cos(phase+w)+e*np.cos(w))+gamma
        vel_model.append(velocity)
        
        residuals.append(vel_obs[i]-velocity/errors_obs[i])
        
    if plot == True:
        # Plot observed vs. modeled phase velocity curve
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.scatter(phase_obs,vel_obs,color='blue',label='Observation')
        ax.scatter(phase_obs,vel_model,color='r',marker='s',label='Model')
        ax.set_xlabel('Phase',fontsize=14)
        ax.set_ylabel('Radial Velocity',fontsize=14)
        ax.set_title('Fitting Model to HD89744',fontsize=18)
        ax.legend()
    
    return(residuals)

# Set parameters for ModelNew function
params = Parameters()
params.add('t0',value=107.0,min=0.01,max=500.)
params.add('P',value=258.0,min=0.01,max=500.)
params.add('e',value=0.70,min=0.01,max=0.999)
params.add('K',value=250.0)
params.add('w',value=190.0,min=0.01,max=360.0)
params.add('gamma',value=-265)

# Print results
results = minimize(ModelNew, params, args=(time_list,phase_obs,vel_obs,errors_obs))
print('Best t0 value = {0:.2f} days'.format(results.params['t0'].value))
print('Best P value = {0:.2f} days'.format(results.params['P'].value))
print('Best e value = {0:.2f}'.format(results.params['e'].value))
print('Best K value = {0:.5f} m/s'.format(results.params['K'].value))
print('Best w value = {0:.2f} deg'.format(results.params['w'].value))
print('Best gamma value = {0:.2f}'.format(results.params['gamma'].value))
print('Chi Squared = {0:.2f}'.format(results.chisqr))
print('Reduced Chi Squared = {0:.2f}'.format(results.redchi))


# In[79]:


# Plot parameters yielded from lmfit (not a good fit and parameters are much different than the published ones)
lmfit_result = Model(time_list,phase_obs,vel_obs,errors_obs,111.99,259.07,.77,875.94,200.30,gamma=-652.6,plot=True)


# In[85]:


# Trying scipy.optimize.curve_fit instead
from scipy.optimize import curve_fit

# Function to calculate model velocities
def ModelScipy(time_list,t0,P,e,K,w,gamma):
    
    # Redefine w in radians
    w *= np.pi/180
    
    # Calculate first term of E and theta relationship
    term1 = ((1+e**2)/(1-e**2))**(1/2)
    
    # Create empty lists of velocity values and residuals
    vel_model = []
    
    # Calculate E at different times
    for i in range(len(time_list)):

        # Select particular time from list
        time = time_list[i]

        # Calculate mean anomaly (M) to use it as first guess for finding E
        M = 2*np.pi/P*(time-t0)

        # Define Kepler equation and its derivative
        f = lambda E: E - e*np.sin(E) - M
        df = lambda E: 1. - e*np.cos(E)

        # Use NewtonRaphson function to find E at each time
        E = NewtonRaphson(f,df,M,1e-6,1000)

        # Convert eccentric anomaly (E) to true anomaly (theta)
        phase = 2*np.arctan(term1*np.tan(E/2.))

        # Calculate velocity with e, theta, K, and w
        velocity = K*(np.cos(phase+w)+e*np.cos(w))+gamma
        vel_model.append(velocity)
        
    return(vel_model)

# Set initial parameter guesses in Kepler's equation
t0 = 107. # in days
P = 258. # in days
e = 0.7

# Set initial parameter guesses in velocity equation
K = 250.0
w = 190
gamma = -270.0

# Calculate model velocitites with initial parameters
y = ModelScipy(time_list,t0,P,e,K,w,gamma)

# Choose bounds for parameters
ranges = ((5.0,1.0,0.01,1.0,0.01,-np.inf),(200.0,400.0,0.999,600,360.0,np.inf))

# Optimize fit
popt, pcov = curve_fit(ModelScipy,time_list,y,sigma=errors_obs,bounds=ranges)
print('Best t0 value = {0:.2f} days'.format(popt[0]))
print('Best P value = {0:.2f} days'.format(popt[1]))
print('Best e value = {0:.2f}'.format(popt[2]))
print('Best K value = {0:.5f} m/s'.format(popt[3]))
print('Best w value = {0:.2f} deg'.format(popt[4]))


# ### These parameters from scipy.optimize.curve_fit are the ones that I gave it so I wasn't able to get it working properly. I've spent quite a lot of time trying to get these optimizers to give me better results, but to no avail so I'm electing to finish the homework here given how many total hours I've already spent on it.

# In[ ]:




