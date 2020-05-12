## Load measured intensity data and define region for calibration

import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy.lib.scimath as npsci  

def Step0_IDT_Init():
    global Nx
    global Ny

    Nx = I_Raw.shape[0]
    Ny = I_Raw.shape[1]
    Nz=Nx

    I_Calib = I_Raw[Cablib_pointX-Cablib_Nx//2-1:Cablib_pointX+Cablib_Nx//2-1, Cablib_pointY-Cablib_Ny//2-1:Cablib_pointY+Cablib_Ny//2-1, : ]
      
    ## Calc the spectrum postion 

    #Each LED geometric postion
    Sorted_Pos = sio.loadmat('Sorted_Pos.mat')
    Sorted_Pos = Sorted_Pos['Sorted_Pos']

    Pos_X = Sorted_Pos[0,0::24//Length_MN]
    Pos_Y = Sorted_Pos[1,0::24//Length_MN]

    Sorted_Pos = np.row_stack((Pos_X,Pos_Y))

    # Frequency coordinate
    Max_frequency = NA/wavelength
    delta_x = 1/(Pixelsize*Nx)      # frequency sampling X.
    delta_y = 1/(Pixelsize*Ny)      # frequency sampling Y.
    delta_z = 1/(Pixelsize*Nz)      # frequency sampling Z.

    fx = np.arange(-np.floor(Nx/2),np.floor((Nx-1)/2)+1)*delta_x # frequency coordinate X.
    fy = np.arange(-np.floor(Ny/2),np.floor((Ny-1)/2)+1)*delta_y # frequency coordinate Y.
    global fx2D
    global fy2D
    (fx2D, fy2D) = np.meshgrid(fx,fy)

    # Generating coordinates on the surface of Ewald Sphere
    fz2D = (npsci.sqrt((n_Medium/wavelength)**2-fy2D**2-fx2D**2)).real

    fxy2D=fx2D+1j*fy2D
    R = np.abs(fxy2D)
    Theta = np.angle(fxy2D)
    
    Aperture = np.logical_not(R>Max_frequency)
    global Aperture_fun
    Aperture_fun = Aperture
 
    # figure
    # imshow(Aperture_fun)

    ##

    #initial guess of led postion and frequecy coord
    
    global Ini_NAx
    global Ini_NAy
    
    Ini_NAx = np.zeros((1,Length_MN))
    Ini_NAy = np.zeros((1,Length_MN))
    Ini_NAz = np.zeros((1,Length_MN))

    Ini_PixelShiftx = np.zeros((1,Length_MN))
    Ini_PixelShifty = np.zeros((1,Length_MN))

    for i in range(1,Length_MN+1):
        aii = Sorted_Pos[0,i-1]
        ajj = Sorted_Pos[1,i-1]
        pic_pos = i
            
        Ini_NAx[0,pic_pos-1] = (aii*Bright_Radius)/(npsci.sqrt((aii*Bright_Radius)**2+(ajj*Bright_Radius)**2+fDL**2))/wavelength
        Ini_NAy[0,pic_pos-1] = (ajj*Bright_Radius)/(npsci.sqrt((aii*Bright_Radius)**2+(ajj*Bright_Radius)**2+fDL**2))/wavelength
        Ini_NAz[0,pic_pos-1] = (npsci.sqrt((n_Medium/wavelength)**2-(Ini_NAx[0,pic_pos-1])**2-(Ini_NAy[0,pic_pos-1])**2)).real
