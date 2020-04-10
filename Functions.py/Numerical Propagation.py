###############################################################################
def Numerical_Propagation(U0,z,pixelsize,wavelength,NA,method):
# function  Uz  = Numerical_Propagation(U0,z,pixelsize,lambda,method)
# Purpose: numerical propagation the complex field to another plane at a
# given distance using 'Angular Specturm' or 'Fresnel' method
#'Inputs': 'U0',Orignal complex field;
#          'z',Propagation distance 
#          'pixelsize'£¬pixelsize('mm'); 
#          'wavelength',Wavelength('mm'); 
#          'method', type of transfer function used ('Angular Spectrum' or
#          'Fresnel')
#'Outputs': 'Uz',Complex field after propagation; 
###############################################################################
# Version 2.0 - 
# Coded by Chao Zuo - 2012-11-9   
# Lastest edited by Chao Zuo - 2014-7-15 
###############################################################################

    import numpy as np
    import math
    
    (N,M) = np.shape(U0)
    
    x = np.arange(1,M+1)
    y = np.arange(1,N+1)
    
    L0X = pixelsize*M
    L0Y = pixelsize*N
    
    k = 2*math.pi/wavelength
    
    u = wavelength*(-M/(L0X*2)+(1/L0X)*(x-1))
    v = wavelength*(-N/(L0Y*2)+(1/L0Y)*(y-1))
    
    uu = u
    vv = v
    
    for i in range(len(v)-1):
        uu = np.row_stack((uu,u))

    for i in range(len(u)-1):
        vv = np.column_stack((vv,v))

    SHy = np.arange(-M/2,M/2)
    SHy = SHy/(M*pixelsize)
    a = SHy
    for i in range(N-1):
        SHy = np.row_stack((SHy,a))
    SHy = SHy*SHy

    SHx = np.arange(-N/2,N/2)
    SHx=SHx/(N*pixelsize)
    b = SHx
    for i in range(M-1):
        SHx = np.column_stack((SHx,b))
    SHx=SHx*SHx

    SH=SHx+SHy

    kk = k/(2*math.pi)*NA
    HHH = np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            if SH[i,j] <= kk**2:
                HHH[i,j] = 1

    FU0 = np.fft.fftshift(np.fft.fft2(U0))

    if method == 'Angular Spectrum':
        H = np.exp(1j*k*z*((1-uu**2-vv**2)**0.5)) # Angular Specturm method 
    elif method == 'Fresnel':
        H = np.exp(1j*k*z*(1-(uu**2+vv**2)/2)) # Fresnel method
    else:
        print('Type of transfer function must be <Angular Spectrum> or <Fresnel>','Error')

    Uz = np.fft.ifft2(np.fft.fftshift(FU0*H*HHH))


##
