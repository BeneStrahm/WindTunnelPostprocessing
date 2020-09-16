# ------------------------------------------------------------------------------
# Description:  Calculating wind speeds at different return periods
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-16
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  
import numpy as np
import scipy as sp
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------
from helpers.pyExtras import getKeyList
from helpers.txtEditor import writeToTxt
# ------------------------------------------------------------------------------
# Abbreviations
# ------------------------------------------------------------------------------

# p...  load
# r...  response
# ms... model scale
# fs... full scale 
# L...  lift
# D...  drag
# M...  moment
# F...  force
# H...  (at) height of building
# sp..  sample
# f...  frequency

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class response(object):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        """Scale time/freq. from model to full scale
        :param modelForce: obj w/ wind tunnel measurements
        :param RPeriod: str w/ return period of wind speed
        :param uH_fs, H_fs: flt w/ full scale building properties
        """
        self.modelForces = modelForces
        self.RPeriod= RPeriod
        self.H_fs   = H_fs
        self.uH_fs  = uH_fs
    
    def scaleTime(self):
        """Scale time/freq. from model to full scale
        """
        # Model properties
        self.modelForces.dT_ms    = 1 / self.modelForces.fsp_ms

        # Scaling factors
        self.lambda_u = self.uH_fs / self.modelForces.uH_ms
        self.lambda_g = self.H_fs / self.modelForces.H_ms
        self.lambda_f = self.lambda_u / self.lambda_g
        self.lambda_t = 1 / self.lambda_f

        # Scale quantities
        self.dT_fs    = self.lambda_t * self.modelForces.dT_ms
        self.fsp_fs   = self.lambda_f * self.modelForces.fsp_ms

    def scaleForces(self):
        """Scale base forces from model to full scale
        """
        # Scaling factors
        self.lambda_F = self.lambda_u ** 2 * self.lambda_g ** 2
        self.lambda_M = self.lambda_u ** 2 * self.lambda_g ** 3

        # Scale floor forces
        self.modelForces.F_p_fs_D= self.modelForces.F_p_ms_D * self.lambda_F
        self.modelForces.F_p_fs_L= self.modelForces.F_p_ms_L * self.lambda_F
        
        # Scale base forces
        self.modelForces.BF_p_fs_D= self.modelForces.BF_p_ms_D * self.lambda_F
        self.modelForces.BF_p_fs_L= self.modelForces.BF_p_ms_L * self.lambda_F
        self.modelForces.BM_p_fs_D= self.modelForces.BM_p_ms_D * self.lambda_M
        self.modelForces.BM_p_fs_L= self.modelForces.BM_p_ms_L * self.lambda_M

    def transToSpectralDomain(self, F_p, dT):
        """Transform time series of forces into the spectral domain
        """
        # Length of time series
        n    = np.shape(F_p)[0]
        N    = n//2 

        # Get the Spectral Density, only positive half
        S_p = abs(np.fft.fft(F_p)[0:N])

        # Compute the power spectra density
        S_p = S_p ** 2 
        
        # According to "Boggs - Wind Loading ...[1991], (p.237)": S(f)+ = 2 * S(f)+/-
        S_p = 2 * S_p

        # Scaling Factor
        S_p = S_p / n

        # Scaling Factor
        S_p = S_p * dT 

        # Compute the frequencies, only positive half
        f = np.fft.fftfreq(n, dT)[0:N]

        return S_p, f

    def calcSpectralResponse(self, f, S_p, omega_e, D):
        """Calculate the response spectrum
        """
        # Length of spectrum
        N = np.shape(S_p)[0]

        # Frequency - Circular frequency
        omega_p = 2 * np.pi * f

        # Apply Dynamic amplification factor
        S_r = np.zeros(N)

        for i in range(0, N):
            eta_i   = omega_p[i]/omega_e
            H_i     = response.mechanicalAdmittance(self, eta_i, D)
            S_r[i]  = abs(H_i)**2 * S_p[i]

        return S_r

    def mechanicalAdmittance(self, eta, D):
        """Mechanical admittance of the 1-DOF System
        """
        # Dynamic amplification factor
        H_f = 1 / np.sqrt((1-eta**2)**2 + (2*D*eta)**2)
        return H_f

    def numericalIntSpectrum(self, dT, S_r):
        """Integrate the response spectrum
        """
        # Length of spectrum
        N = np.shape(S_r)[0]

        # Nyquist frequency
        f_nyq = 1 / (2 * dT)

        # Sample spacing df
        df = f_nyq / N 

        # Perform the numerical integration with Simpson rule
        F_ms = sp.integrate.simps(S_r, dx=df) 

        # Return directly RMS
        F_rms = np.sqrt(F_ms)
        return F_rms

    def peakFactor(self, fe, T):
        """Compute the peak factor
        """
        g_peak = np.sqrt(2 * np.log(fe * T)) + 0.5772 / np.sqrt((2 * np.log(fe * T)))
        return g_peak


class BaseResponseForces(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)

    def calcLoadStats(self):
        # Investigated wind speed
        writeToTxt("loadStats.txt", "u_H_mean:      " + '{:02.3f}'.format(self.uH_fs))
        writeToTxt("loadStats.txt", "Return Period: " + self.RPeriod)
        
        # Base moment stats, drag direction
        writeToTxt("loadStats.txt", "Base moment in drag direction [kNm]")
        self.BM_p_fs_D_mean = np.mean(self.modelForces.BM_p_fs_D)
        writeToTxt("loadStats.txt", "F_p_mean:      " + '{:02.3f}'.format(self.BM_p_fs_D_mean))

        self.BM_p_fs_D_max  = np.max(self.modelForces.BM_p_fs_D)
        writeToTxt("loadStats.txt", "F_p_max:       " + '{:02.3f}'.format(self.BM_p_fs_D_max))

        self.BM_p_fs_D_min  = np.min(self.modelForces.BM_p_fs_D)
        writeToTxt("loadStats.txt", "F_p_min:       " + '{:02.3f}'.format(self.BM_p_fs_D_min))

        self.BM_p_fs_D_std  = np.std(self.modelForces.BM_p_fs_D)
        writeToTxt("loadStats.txt", "F_p_std:       " + '{:02.3f}'.format(self.BM_p_fs_D_std))

        # Base moment stats, lift direction
        writeToTxt("loadStats.txt", "Base moment in lift direction [kNm]")
        self.BM_p_fs_L_mean = np.mean(self.modelForces.BM_p_fs_L)
        writeToTxt("loadStats.txt", "F_p_mean:      " + '{:02.3f}'.format(self.BM_p_fs_L_mean))

        self.BM_p_fs_L_max  = np.max(self.modelForces.BM_p_fs_L)
        writeToTxt("loadStats.txt", "F_p_max:       " + '{:02.3f}'.format(self.BM_p_fs_L_max))

        self.BM_p_fs_L_min  = np.min(self.modelForces.BM_p_fs_L)
        writeToTxt("loadStats.txt", "F_p_min:       " + '{:02.3f}'.format(self.BM_p_fs_L_min ))

        self.BM_p_fs_L_std  = np.std(self.modelForces.BM_p_fs_L)
        writeToTxt("loadStats.txt", "F_p_std:       " + '{:02.3f}'.format(self.BM_p_fs_L_std))

class TipResponseDeflections(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)

class TipResponseAccelerations(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
