# ------------------------------------------------------------------------------
# Description:  Calculating base force coefficients from pressure coefficients
#               Used for the TPU Aerodynamic Database experiments
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-15
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  
import numpy as np
import hdf5storage as h5s    
import pandas as pd                    
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

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

class modelForces(object):
    def __init__(self, fname, rho_air_ms=1.25):
        """Loading wind tunnel data 
        :param fname: str w/ full path to file
        :param rho_air_ms: flt w/ air density of wind tunnel
        """        
        self.fname      = fname
        self.rho_air_ms = rho_air_ms

    def loadWindtunnelData(self):
        """Loading wind tunnel data and calculating forces from surface 
        pressure coeffs in model scale
        """        
        # Load matlab file
        mat         = h5s.loadmat(self.fname)

        # Pressure coeff.
        self.Cp     = mat['Wind_pressure_coefficients']
        self.Loc    = mat['Location_of_measured_points']

        # Wind speeds
        self.uH_ms  = float(mat['Uh_AverageWindSpeed'][0])

        # Time scales
        self.fsp_ms = float(mat['Sample_frequency'][0])
        self.T_ms   = float(mat['Sample_period'][0])
        self.nT     = np.shape(self.Cp)[0]

        # Geometrie [in m]
        self.H_ms   = float(mat['Building_height'][0])

        # Coordinates [in m], number and face of measurement points
        self.x      = self.Loc[0]
        self.z      = self.Loc[1]
        self.n      = self.Loc[2]
        self.face   = self.Loc[3]

        # Get levels
        self.z_lev  = list(set(self.z))
        self.nz     = len(self.z_lev)
        
        # Measurement area [in m]
        self.A_i    = 0.02 * 0.02

        # Calculate forces [in kN]
        self.F_p_ms = self.A_i * self.Cp * 0.5 * self.rho_air_ms  * self.uH_ms  ** 2 * 10 ** -3
        
    def calcModelBaseForces(self):
        """Calculating base forces from surface pressure coeffs in model scale
        """
        # Sort in drag / lift base forces
        self.BF_p_ms_D, self.BF_p_ms_L= np.zeros(self.nT), np.zeros(self.nT)
        self.BM_p_ms_D, self.BM_p_ms_L= np.zeros(self.nT), np.zeros(self.nT)
        
        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p_ms.T, 0):
            # Face 1 -> +DragF, +DragM
            if self.face[i] == 1:
                self.BF_p_ms_D = self.BF_p_ms_D + F_i
                self.BM_p_ms_D = self.BM_p_ms_D + F_i * self.z[i] 
            # Face 2 -> +LiftF, -LiftM
            elif self.face[i] == 2:
                self.BF_p_ms_L = self.BF_p_ms_L + F_i
                self.BM_p_ms_L = self.BM_p_ms_L - F_i * self.z[i] 
            # Face 3 -> -DragF, -DragM
            elif self.face[i] == 3:
                self.BF_p_ms_D = self.BF_p_ms_D - F_i
                self.BM_p_ms_D = self.BM_p_ms_D - F_i * self.z[i] 
            # Face 4 -> -LiftF, +LiftM
            elif self.face[i] == 4:
                self.BF_p_ms_L = self.BF_p_ms_L - F_i
                self.BM_p_ms_L = self.BM_p_ms_L + F_i * self.z[i] 

    def calcModelFloorForces(self):
        """Calculating floor forces from surface pressure coeffs in model scale
        """
        # Sort in drag / lift forces for each floor
        self.F_p_ms_D, self.F_p_ms_L= np.zeros((self.nz, self.nT)), np.zeros((self.nz, self.nT))

        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p_ms.T, 0):
            # get index where to sort to
            j = self.z_lev.index(self.z[i])

            # Face 1 -> +DragF, +DragM
            if self.face[i] == 1:
                self.F_p_ms_D[j,:] = self.F_p_ms_D[j,:] + F_i
            # Face 2 -> +LiftF, -LiftM
            elif self.face[i] == 2:
                self.F_p_ms_L[j,:] = self.F_p_ms_L[j,:] + F_i
            # Face 3 -> -DragF, -DragM
            elif self.face[i] == 3:
                self.F_p_ms_D[j,:] = self.F_p_ms_D[j,:] - F_i
           # Face 4 -> -LiftF, +LiftM
            elif self.face[i] == 4:
                self.F_p_ms_L[j,:] = self.F_p_ms_L[j,:] - F_i

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
