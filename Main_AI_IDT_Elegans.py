import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import time
import scipy.io as sio
import numpy.lib.scimath as npsci

# This version uses Step0 and Step2, omits Step1
# Plotting and Step2_IDT_Poss require edits; see below.

## Parameters

wavelength = 0.515 # Wavelength
k=2*math.pi/wavelength # Wave number

NA = 0.65
Mag = 40
Pixelsize = 6.5/Mag
n_Medium = 1.33

Length_MN = 8 # the number of LEDs
Bright_Radius = 29.7 # the Radius of LED ring
fDL = Bright_Radius / np.tan(np.arcsin(NA)) # the distance of LED and object 

Calib = 0 # not calibrating LED postion
gpu = 0 # if using gpu

## Step Load measured intensity data and initial spectrum postion

I_Raw = sio.loadmat('IRaw_Elegans.mat')
I_Raw = I_Raw['I_Raw']

Cablib_Nx=400  
Cablib_Ny=400

Cablib_pointX=300
Cablib_pointY=300

Step0_IDT_Init()

## Step Implete the calibation of LED postion

# Calibration step skipped

## Implete the IDT
dz=1
Depth_Set = np.arange(-10,11,dz)

Alpha=1e2
Beta=1e2

Step2_IDT_Poss(Ini_NAx,Ini_NAy) # In need of some edits, see file.

## show RI Slice

# Plotting below is barely translated from Matlab. 
# Likely needs adjustment to properly function in Python.
plt.gray()
for ii in range(1,len(Depth_Set)+1):
    plt.subplot(121)
    plt.imshow((np.squeeze(RI[:,:,ii-1])).real)
    plt.subplot(122)
    plt.imshow((np.squeeze(RI[:,:,ii-1])).imag)
    time.sleep(0.5)
