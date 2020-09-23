# ------------------------------------------------------------------------------
# Description:  
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-22
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  
import hdf5storage as h5s  
import numpy as np
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
# fq...  frequency
# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  
class wtModelAeroForces:
    """Class containing the model scale forces of the wind tunnel
    :cvar Cp: time series of wind pressure coefficients at measured points
    :vartype Cp: list[float, float]
    :cvar F_p: time series of wind forces at measured points
    :vartype F_p: list[float, float]
    :cvar BF_p_D: time series of base shear in drag direction
    :vartype BF_p_D: np.arr[float]
    :cvar BF_p_L: time series of base shear in lift direction
    :vartype BF_p_L: np.arr[float]
    :cvar BM_p_D: time series of base moment in drag direction
    :vartype BM_p_D: np.arr[float]
    :cvar BM_p_L: time series of base moment in lift direction
    :vartype BM_p_L: np.arr[float]
    :cvar F_p_D: time series of forces at each level in drag direction
    :vartype F_p_D: np.arr[float, float]
    :cvar F_p_L: time series of forces at each level in lift direction
    :vartype F_p_L: np.arr[float, float]
    """
    def __init__(self):
        """Inits the class
        """
        pass

    def loadTPUModelForces(self, wtModelProp):
        """Loading measured wind tunnel model forces
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        """
        # Load matlab file
        mat         = h5s.loadmat(wtModelProp.fname)
        
        # Pressure coeff.
        self.Cp     = mat['Wind_pressure_coefficients']

        # Calculate forces [in kN]
        self.F_p = wtModelProp.A_i * self.Cp * 0.5 * wtModelProp.rho_air  * wtModelProp.uH  ** 2 * 10 ** -3

        # Calculating base forces
        self.calcModelBaseForces(wtModelProp)

        # Calculating floor forces
        self.calcModelFloorForces(wtModelProp)

    def calcModelBaseForces(self, wtModelProp):
        """Calculating base forces from wind tunnel model forces in model scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        """
        # Sort in drag / lift base forces
        self.BF_p_D, self.BF_p_L= np.zeros(wtModelProp.nT), np.zeros(wtModelProp.nT)
        self.BM_p_D, self.BM_p_L= np.zeros(wtModelProp.nT), np.zeros(wtModelProp.nT)
        
        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p.T, 0):
            # Face 1 -> +DragF, +DragM
            if wtModelProp.face[i] == 1:
                self.BF_p_D = self.BF_p_D + F_i
                self.BM_p_D = self.BM_p_D + F_i * wtModelProp.z[i] 
            # Face 2 -> +LiftF, -LiftM
            elif wtModelProp.face[i] == 2:
                self.BF_p_L = self.BF_p_L + F_i
                self.BM_p_L = self.BM_p_L - F_i * wtModelProp.z[i] 
            # Face 3 -> -DragF, -DragM
            elif wtModelProp.face[i] == 3:
                self.BF_p_D = self.BF_p_D - F_i
                self.BM_p_D = self.BM_p_D - F_i * wtModelProp.z[i] 
            # Face 4 -> -LiftF, +LiftM
            elif wtModelProp.face[i] == 4:
                self.BF_p_L = self.BF_p_L - F_i
                self.BM_p_L = self.BM_p_L + F_i * wtModelProp.z[i] 

    def calcModelFloorForces(self, wtModelProp):
        """Calculating floor forces from wind tunnel model forces in model scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        """
        # Sort in drag / lift forces for each floor
        self.F_p_D, self.F_p_L= np.zeros((wtModelProp.nz, wtModelProp.nT)), np.zeros((wtModelProp.nz, wtModelProp.nT))

        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p.T, 0):
            # get index where to sort to
            j = np.where(wtModelProp.z_lev == wtModelProp.z[i])
            # Face 1 -> +DragF, +DragM
            if wtModelProp.face[i] == 1:
                self.F_p_D[j,:] = self.F_p_D[j,:] + F_i
            # Face 2 -> +LiftF, -LiftM
            elif wtModelProp.face[i] == 2:
                self.F_p_L[j,:] = self.F_p_L[j,:] + F_i
            # Face 3 -> -DragF, -DragM
            elif wtModelProp.face[i] == 3:
                self.F_p_D[j,:] = self.F_p_D[j,:] - F_i
           # Face 4 -> -LiftF, +LiftM
            elif wtModelProp.face[i] == 4:
                self.F_p_L[j,:] = self.F_p_L[j,:] - F_i

class buildAeroForces:
    """Class containing the building scale forces of the wind tunnel
    :cvar BF_p_D: time series of base shear in drag direction
    :vartype BF_p_D: np.arr[float]
    :cvar BF_p_L: time series of base shear in lift direction
    :vartype BF_p_L: np.arr[float]
    :cvar BM_p_D: time series of base moment in drag direction
    :vartype BM_p_D: np.arr[float]
    :cvar BM_p_L: time series of base moment in lift direction
    :vartype BM_p_L: np.arr[float]
    :cvar F_p_D: time series of forces at each level in drag direction
    :vartype F_p_D: np.arr[float, float]
    :cvar F_p_L: time series of forces at each level in lift direction
    :vartype F_p_L: np.arr[float, float]
    """
    def __init__(self, scalingFactors, wtModelAeroForces):
        """Inits the class
        :param scalingFactors: scaling factors from wind tunnel to building scale
        :type scalingFactors: :class:`~scaling.scalingFactors`
        :param scalingFactors: full scale properties of the building
        :type scalingFactors: :class:`~aeroForces.wtModelAeroForces`
        """
        # Scale forces
        self.scaleForces(scalingFactors, wtModelAeroForces)

    def scaleForces(self, scalingFactors, wtModelAeroForces):
        """Scale base forces from model to full scale
        :param scalingFactors: scaling factors from wind tunnel to building scale
        :type scalingFactors: :class:`~scaling.scalingFactors`
        :param scalingFactors: full scale properties of the building
        :type scalingFactors: :class:`~aeroForces.wtModelAeroForces`
        """
        # Scale floor forces
        self.F_p_D = wtModelAeroForces.F_p_D * scalingFactors.lambda_F
        self.F_p_L = wtModelAeroForces.F_p_L * scalingFactors.lambda_F
        
        # Scale base forces
        self.BF_p_D = wtModelAeroForces.BF_p_D * scalingFactors.lambda_F 
        self.BF_p_L = wtModelAeroForces.BF_p_L * scalingFactors.lambda_F
        self.BM_p_D = wtModelAeroForces.BM_p_D * scalingFactors.lambda_M
        self.BM_p_L = wtModelAeroForces.BM_p_L * scalingFactors.lambda_M

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   