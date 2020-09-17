# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt


class Planets:
    # Intialize important user-input information
    def __init__(self,filename,x,y):
       
        self.file = filename
        self.x = x
        self.y = y

    # Function to read text files
    def Read(self,filename):
    
        # Read and return data from NASA Exoplanet Archive of confirmed exoplanets
        data = np.genfromtxt(filename,dtype=None,delimiter=',',skip_header=96,names=True,invalid_raise=False,encoding=None)
        return data

    # Function to plot different planet properties vs. each other
    def Plot(self):

        # Save data using 'Read' function
        data = self.Read(self.file)
    
        if self.x == 'mass' and self.y == 'sma':
            xdata = data['pl_bmasse']
            ydata = data['pl_orbsmax']
        plt.scatter(xdata,ydata)

test = Planets('ConfirmedExoplanets.csv','mass','sma')
test.Plot()
