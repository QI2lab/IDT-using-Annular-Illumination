##### IDT poss after LED Postion calibration #####
# Requires edits; see comments below.

import numpy as np
import math
import matplotlib.pyplot as plt
import time
import numpy.lib.scimath as npsci  

def Step2_IDT_Poss(Ini_NAx,Ini_NAy):

    # apply calibrated
    if Calib == 1:
        Ini_NAx = (np.transpose(freqXY3[:,1]) - (Cablib_Nx/2+1))/Cablib_Nx/Pixelsize
        Ini_NAy = (np.transpose(freqXY3[:,0]) - (Cablib_Ny/2+1))/Cablib_Ny/Pixelsize
    else:
        Ini_NAx = np.squeeze(Ini_NAx)
        Ini_NAy = np.squeeze(Ini_NAy)
        TempNA = Ini_NAx
        Ini_NAx = -Ini_NAy
        Ini_NAy = TempNA

    # correct frequency coord of LED 
    Ini_PixelShiftx = np.zeros((1,Length_MN))
    Ini_PixelShifty = np.zeros((1,Length_MN))

    for i in range (1,Length_MN+1):
        pic_pos = i
       
        Ini_PixelShiftx[0,pic_pos-1] = np.round(Ini_NAx[pic_pos-1]*Pixelsize*Nx)
        Ini_PixelShifty[0,pic_pos-1] = np.round(Ini_NAy[pic_pos-1]*Pixelsize*Ny)

    ## Generate Complex or pure phase 3D object, and correspoding intensity bright field iamges

    print('Calculating slice-wise transfer functions...')
    PTF_4D = np.zeros((Nx,Ny,len(Depth_Set),Length_MN), dtype=complex)
    ATF_4D = np.zeros((Nx,Ny,len(Depth_Set),Length_MN), dtype=complex)

    PTF_3D = np.zeros((Nx,Ny,len(Depth_Set)), dtype=complex)
    ATF_3D = np.zeros((Nx,Ny,len(Depth_Set)), dtype=complex)

    # Note: MatLab program converts floating point to 32-bit single precision in the above four lines

    ##

    uu = fx2D
    vv = fy2D
    k_Wavenum = k
    k_Medium = k_Wavenum*n_Medium

    #figure
    for i in range (1,Length_MN+1):
        pic_pos = i

        v_s = Ini_NAx[pic_pos-1]
        u_s = Ini_NAy[pic_pos-1]
       
        G = (1/(k_Medium*npsci.sqrt(1-(wavelength**2)*((uu-u_s)**2+(vv-v_s)**2)))).real
        Gf = (1/(k_Medium*npsci.sqrt(1-(wavelength**2)*((uu+u_s)**2+(vv+v_s)**2)))).real
       
        Pupil = np.roll(Aperture_fun, (int(Ini_PixelShiftx[0,pic_pos-1]),int(Ini_PixelShifty[0,pic_pos-1])), (0,1))
        Pupilf = np.roll(Aperture_fun, (int(-Ini_PixelShiftx[0,pic_pos-1]),int(-Ini_PixelShifty[0,pic_pos-1])), (0,1))
       
        uv_vector1 = (npsci.sqrt(1-(wavelength**2)*((uu-u_s)**2+(vv-v_s)**2))).real
        uv_vector2 = (npsci.sqrt(1-(wavelength**2)*((uu+u_s)**2+(vv+v_s)**2))).real
        uv_s = npsci.sqrt(1-(wavelength**2)*(u_s**2+v_s**2))
       
        for j in range (1,len(Depth_Set)+1):
            PTF_3D[:,:,j-1] = \
            (Pupil * np.sin(k_Medium*Depth_Set[j-1]*(uv_vector1 - uv_s))*G + \
                Pupilf * np.sin(k_Medium*Depth_Set[j-1]*(uv_vector2 - uv_s))*Gf) + \
            1j * (Pupil * np.cos(k_Medium*Depth_Set[j-1]*(uv_vector1 - uv_s))*G - \
                Pupilf * np.cos(k_Medium*Depth_Set[j-1]*(uv_vector2 - uv_s))*Gf)
       
            ATF_3D[:,:,j-1] = \
            -(Pupil * np.cos(k_Medium*Depth_Set[j-1]*(uv_vector1 - uv_s))*G + \
                Pupilf * np.cos(k_Medium*Depth_Set[j-1]*(uv_vector2 - uv_s))*Gf) + \
            1j * (Pupil * np.sin(k_Medium*Depth_Set[j-1]*(uv_vector1 - uv_s))*G - \
                Pupilf * np.sin(k_Medium*Depth_Set[j-1]*(uv_vector2 - uv_s))*Gf)      
           
            PTF_3D[:,:,j-1] = np.fft.fftshift(PTF_3D[:,:,j-1])
            ATF_3D[:,:,j-1] = np.fft.fftshift(ATF_3D[:,:,j-1])

        PTF_4D[:,:,:,pic_pos-1] = 0.5*dz*k_Wavenum**2 * PTF_3D
        ATF_4D[:,:,:,pic_pos-1] = 0.5*dz*k_Wavenum**2 * ATF_3D


    ## Calculate Eq.(7) and Eq.(8) in paper
    
    sum_PTF = 0
    sum_ATF = 0

    conj_PTF_Iten = 0
    conj_ATF_Iten = 0

    conj_term1 = 0
    conj_term2 = 0

    tic = time.perf_counter()
    for i in range(1, Length_MN+1):
        pic_pos = i

        Itmp = I_Raw[:,:,pic_pos-1]
        Ihat_tmp = np.fft.fft2((Itmp))

        sum_PTF = sum_PTF + np.abs(PTF_4D[:,:,:,pic_pos-1])**2
        sum_ATF = sum_ATF + np.abs(ATF_4D[:,:,:,pic_pos-1])**2

        repmatIhat = Ihat_tmp
        for j in range(len(Depth_Set)-1):
            repmatIhat = np.dstack((repmatIhat, Ihat_tmp))

        conj_PTF_Iten = conj_PTF_Iten + np.conj(PTF_4D[:,:,:,pic_pos-1]) * repmatIhat
        conj_ATF_Iten = conj_ATF_Iten + np.conj(ATF_4D[:,:,:,pic_pos-1]) * repmatIhat

        conj_term1 = conj_term1 + np.conj(PTF_4D[:,:,:,pic_pos-1]) * ATF_4D[:,:,:,pic_pos-1]
        conj_term2 = conj_term2 + np.conj(ATF_4D[:,:,:,pic_pos-1]) * PTF_4D[:,:,:,pic_pos-1]

    ## Repeat IDT algorithm

    global V_im
    
    Normalized_term = (sum_PTF+Alpha) * (sum_ATF+Beta) - (conj_term1 * conj_term2)
    V_re = ((sum_ATF+Beta) * conj_PTF_Iten - conj_term1 * conj_ATF_Iten) / Normalized_term
    V_im = ((sum_PTF+Alpha) * conj_ATF_Iten - conj_term2 * conj_PTF_Iten) / Normalized_term

    toc = time.perf_counter()
    print('Performing IDT spends...'+str(toc-tic)+' seconds')
    
    # NOTE: Up to this point, variables appear to be accurate.
    # Namely, V_im and V_re and consistent with their corresponding Matlab values.
    # However, np.fft.ifft2(V_im) and np.fft.ifft2(V_re) gives arrays different
    # from those of ifft2(V_im) and ifft2(V_re) (as performed in Matlab).
    # Consequently, all subsequent Python variables (v_im, v_re, n_re, n_im,
    # and most importantly, RI) do not match those of Matlab.
    # This discrepancy requires correction in order to produce accurate and coherent images.
    
    v_im = (np.fft.ifft2(V_im)).real
    v_re = (np.fft.ifft2(V_re)).real
    
    n_re = npsci.sqrt(((n_Medium**2 + v_re) + npsci.sqrt((n_Medium**2 + v_re)**2 + v_im**2)) / 2)
    n_im = (v_im / n_re) / 2

    global RI
    RI = n_re+(1j*n_im)
