#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""
    
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Function to plot the solution
def plot(x_list, y_list, funcLambdify, barTitle, nContourLines = 29, colorBarFormat = "%.f", plotTitle=None):
    
    """ 
        Plots x,y data in a countor plot.
        Plots the results of T (x,y) in a color bar
        Note: Considers results to be in steady state
    """
    
    # Data for the plot
    x_value, y_value= np.meshgrid(x_list, y_list)
    z =  funcLambdify(x_value, y_value)
    
    # If z is a uniform field
    # Need to input some variation in the data for the contour plot to work
    if (isinstance(z, (int,float))):
        SMALL = 1e-12
        z = z*np.ones(x_value.shape)
        z[0] = z[0]-z[0]*SMALL 
        z[-1] = z[-1] + z[-1]*SMALL
    
    
    ## Plot
    # Plot properties
    fontSize = 15
    lineWidth = 2.0
    
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')
    plt.rcParams["figure.figsize"] = (8,6)
    plt.rcParams['axes.linewidth'] = lineWidth        
    
    
    fig,ax=plt.subplots(1,1)

    contourLevels = np.linspace(np.min(z), np.max(z), nContourLines, endpoint=True)
        
    colorOfMap = plt.cm.jet
    
       
    # Plots contour lines
    ax.contour(x_value, y_value, z, contourLevels, linewidths=0.5, colors='k')  
    
    # Plot contour map
    cp = ax.contourf(x_value, y_value, z, contourLevels, cmap=colorOfMap )
    
    # Add a colorbar with 8 entries to the plot
    tickMarks = []
    for i in range(0, len(contourLevels), 4):
        tickMarks.append(contourLevels[i])
    
           
    cbar = fig.colorbar(cp, ticks=tickMarks, format=colorBarFormat)
    
    # Size of ticks in color bar
    cbar.ax.tick_params(labelsize=fontSize)
    
    # Size of label in color bar
    cbar.ax.set_ylabel(barTitle, fontsize=fontSize +4)
    
    if (plotTitle is not None):
        ax.set_title(plotTitle, fontsize=fontSize+8)
        
    ax.tick_params(which = 'both', direction='out', length=5, width=lineWidth)
    plt.xlabel("x [m]", fontsize=fontSize+4)
    plt.ylabel("y [m]", fontsize=fontSize+4)
    plt.yticks(fontsize = fontSize+2)
    plt.xticks(fontsize = fontSize+2)
    plt.grid(visible=True, which='major', color='black', linestyle='-')
    plt.minorticks_on()
    plt.grid(visible=True, which='minor', color='darkgrey', linestyle='--', alpha=0.5)
    plt.tight_layout()
    

    # Add logo
    logoPath = 'pyMMSFoam/Logo/logo.png'
    logo = Image.open(logoPath)
    newax = fig.add_axes([0.84, 0.01, 0.15, 0.15], anchor='NE', zorder=-1)
    newax.imshow(logo, alpha=0.2)
    newax.axis('off')
    plt.subplots_adjust(left=0.075, bottom=0.15, right=0.95, top=0.95)

    plt.show()
    

