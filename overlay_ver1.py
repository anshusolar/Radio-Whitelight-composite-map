#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:19:54 2019

@author: anshu
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import tkinter as tk
from tkinter import Tk,Frame,messagebox,Menu
import os
from tkinter import filedialog
import cv2
import os, glob
import random
import math
from os import path
import sys
import subprocess
import pdb
import sunpy
import sunpy.data.sample
import sunpy.map
from sunpy.sun import constants as con
from sunpy.map import header_helper
from astropy.io import fits
import astropy.units as u
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord
from scipy import ndimage, misc
import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import math
from tkinter.constants import TOP, BOTH, BOTTOM, LEFT
from tkinter.ttk import *
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#import pylab


root= tk.Tk()
root.style = Style()
root.style.theme_use("clam")
root.minsize(1200,600)
root.title('Overlay Window')
mainFrame = tk.Frame(root)
mainFrame.grid()
entryFrame = tk.Frame(mainFrame, width=454, height=20)
entryFrame.grid(row=0, column=1)
entryFrame.columnconfigure(0, weight=10)  
entryFrame.grid_propagate(False)





def OpenLASCO():
    global LASCO_file

    LASCO_file=filedialog.askopenfile(initialdir='pwd',title='Select LASCO png file')
    if '.png' in LASCO_file.name:
        print('You have selected LASCO png file correctly')
        
    else:
        messagebox.showerror('Extension Error','You chose a wrong file!! Please choose a png file')
    return None
  
def OpenGRAPH():
    global GRAPH_file

    GRAPH_file=filedialog.askopenfile(initialdir='pwd',title='Select GRAPH fits file')
    if '.fits' in GRAPH_file.name:
        print('You have selected GRAPH fits file correctly')
        
    else:
        messagebox.showerror('Extension Error','You chose a wrong file!! Please choose a png file')
    return None

    
    
def plotLASCO():
    
    global LASCO_data 
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    
    figure = Figure(figsize=(12, 6), dpi=100)
    plot = figure.add_subplot(1, 1, 1)
    canvas = FigureCanvasTkAgg(figure, root)
    canvas.get_tk_widget().grid(row=0, column=0)
    
    plot.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    
    circle1 = plt.Circle((255, 244), 41, color='b', fill=False)
    plot.set_title('Whitelight CME LASCO-C2')
    plot.set_xlabel('X (pixels)')
    plot.set_ylabel('Y (pixels)')   

    return LASCO_data
                
            
def plotGRAPH():
    
    global grh
    global hdu
    global k
    global grh_rot
    global lev
    
    
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    print(a.deg)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    
    
    RA_c=hdu.header['CRVAL1']
    DEC_c=hdu.header['CRVAL2']
    
    figure = Figure(figsize=(12, 6), dpi=100)
    plot = figure.add_subplot(1, 1, 1)
    canvas = FigureCanvasTkAgg(figure, root)
    canvas.get_tk_widget().grid(row=0, column=0)
    
    plot.contour(grh_rot, lev, extend='both')
    circle1 = plt.Circle((255, 255), 41, color='b', fill=False)
    plot.add_artist(circle1)
    plot.set_aspect('equal', 'box')
    plot.set_title('Radio contours for Type II Burst')
    plot.set_xlabel('X (pixels)')
    plot.set_ylabel('Y (pixels)')   
    
    return grh,k,a,grh_rot
    
def plotAll():
    
    global LASCO_data
    global grh
    global hdu
    global k
    global grh_rot
    global lev
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    
    figure = Figure(figsize=(12, 6), dpi=100)
    plot = figure.add_subplot(1, 1, 1)
    canvas = FigureCanvasTkAgg(figure, root)
    canvas.get_tk_widget().grid(row=0, column=0)
    plot.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    plot.contour(grh_rot, lev, extend='both')
    circle1 = plt.Circle((255, 244), 41, color='b', fill=False)
    plot.add_artist(circle1)
    plot.set_title('Composite Image of Whitelight CME and Radio contours for Type II Burst')
    plot.set_xlabel('X (pixels)')
    plot.set_ylabel('Y (pixels)')

    
plt.ioff()   
    
def LASCO_png():
    
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    plt.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    plt.title('Whitelight CME LASCO-C2')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')
    plt.savefig(LASCO_file.name+'new'+'.png')
    plt.close()

def LASCO_pdf():
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    plt.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    plt.title('Whitelight CME LASCO-C2')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')
    plt.savefig(LASCO_file.name+'.pdf')  
    plt.close()
    
def GRAPH_png():
    
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    print(a.deg)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    

    plt.contour(grh_rot, lev, extend='both') 
    #plt.aspect('equal', 'box')
    plt.title('Radio contours for Type II Burst')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')   
    circle1 = plt.Circle((255, 255), 41, color='b', fill=False)
    #plt.plot(circle1)
    plt.savefig(GRAPH_file.name+'.png') 
    plt.close()
        

def GRAPH_pdf():
    
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    print(a.deg)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    plt.contour(grh_rot, lev, extend='both') 
    plt.title('Radio contours for Type II Burst')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')
    plt.savefig(GRAPH_file.name+'.pdf')  
    plt.close()
    
    
def Composite_png():
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    plt.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    plt.contour(grh_rot, lev, extend='both') 
    plt.title('Composite Image of Whitelight CME and Radio contours for Type II Burst')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')
    plt.savefig(k+'.png')
    plt.close()

def Composite_pdf():
    
    LASCO_data = cv2.imread(LASCO_file.name,0)
    hdulist = fits.open(GRAPH_file.name)
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape
    grh1=hdu.data
    grh=grh1[0,0,:,:]

    lev=np.linspace(10,100,12)*np.max(grh)/100
    k=hdu.header['DATE-OBS']
    a=sunpy.coordinates.sun.P(k)
    grh_rot = ndimage.rotate(grh, a.deg, reshape=False)
    plt.imshow(LASCO_data, cmap = 'gray', interpolation = 'bicubic')
    plt.contour(grh_rot, lev, extend='both') 
    plt.title('Composite Image of Whitelight CME and Radio contours for Type II Burst')
    plt.xlabel('X (pixels)')
    plt.ylabel('Y (pixels)')
    plt.savefig(k+'.pdf')  
    plt.close()

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate


def close():
        if messagebox.askyesno("Overlay plotting tool","Do you really want to exit?"):
            root.destroy()
        else:
            return None





menu = tk.Menu(root)

mfile = tk.Menu(menu)
mfile.add_command(label="OpenLASCO", command=OpenLASCO)
mfile.add_command(label="OpenGRAPH", command=OpenGRAPH)


mplot = tk.Menu(menu)
mplot.add_command(label="Plot LASCO", command=plotLASCO)
mplot.add_command(label="Plot GRAPH", command=plotGRAPH)
mplot.add_command(label="Overlay", command=plotAll)

manalysis = tk.Menu(menu)
manalysis.add_command(label="HT_plot", command=plotLASCO)
manalysis.add_command(label="Width", command=plotGRAPH)
manalysis.add_command(label="...")

msave = tk.Menu(menu)
msave.add_command(label="LASCO_png", command=LASCO_png)
msave.add_command(label="LASCO_pdf", command=LASCO_pdf)
msave.add_command(label="GRAPH_png", command=GRAPH_png)
msave.add_command(label="GRAPH_pdf", command=GRAPH_pdf)
msave.add_command(label="Composite_png", command=Composite_png)
msave.add_command(label="Composite_pdf", command=Composite_pdf)
msave.add_command(label="...")


menu.add_cascade(label='File', menu=mfile)
menu.add_cascade(label='Plot', menu=mplot)
menu.add_cascade(label='Analysis',menu=manalysis)
menu.add_cascade(label='Save',menu=msave)
menu.add_command(label="Close", command=close)


root.config(menu=menu)

#root.mainloop()


