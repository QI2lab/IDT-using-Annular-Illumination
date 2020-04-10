def model_TF(source_pattern,deltaz,Pixelsize,wavelength,NA):
    #UNTITLED Summary of this function goes here
    #   Detailed explanation goes here

    import numpy as np
    import math

    (Nx,Ny) = np.shape(source_pattern)

    AA_Pha = np.ones((Nx,Ny))
    PP_Pha = np.zeros((Nx,Ny))
    PP_Pha[int(Nx/2),int(Ny/2)] = 0.1
    PurePha_Obj = AA_Pha*(np.exp(1j*PP_Pha))

    AA_Abs = np.ones((Nx,Ny))
    AA_Abs[int(Nx/2),int(Ny/2)] = 0.9
    PP_Abs = np.zeros((Nx,Ny))
    PureAbs_Obj = AA_Abs*(np.exp(1j*PP_Abs))

    S = np.fft.ifft2(np.fft.fftshift(source_pattern))
    SS = np.angle(S)
    U0_Pha = PurePha_Obj*(np.exp(1j*SS))
    U0_Abs = PureAbs_Obj*(np.exp(1j*SS))

    Uz1 = Numerical_Propagation(U0_Pha,deltaz,Pixelsize,wavelength,NA,'Angular Spectrum')
    Uz2 = Numerical_Propagation(U0_Abs,deltaz,Pixelsize,wavelength,NA,'Angular Spectrum')

    Iz1 = np.abs(Uz1)**2
    Iz2 = np.abs(Uz2)**2

    PTF_Output = np.fft.fft2(Iz1-1)/(np.fft.fft2(PP_Pha))/2
    ATF_Output = np.fft.fft2(Iz2+1)/(np.fft.fft2(AA_Abs))/2
    # Note: MatLab program converts floating point to 32-bit single precision in the above two lines
    
    return PTF_Output, ATF_Output
