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
# fq... frequency
# dn... direction
# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  
class wtModelAeroForces:
    """Class containing the model scale forces of the wind tunnel
    :cvar Cp: time series of wind pressure coefficients at measured points
    :vartype Cp: list[float, float]
    :cvar F_p: time series of wind forces at measured points
    :vartype F_p: list[float, float]
    :cvar BF_p: time series of base shear 
    :vartype BF_p: np.arr[float]
    :cvar BM_p: time series of base moment 
    :vartype BM_p: np.arr[float]
    :cvar LF_p: time series of forces at each level
    :vartype LF_p: np.arr[float, float]
    """
    def __init__(self):
        """Inits the class
        """
        pass

    def loadTPUModelForces(self, wtModelProp, buildProp):
        """Loading measured wind tunnel model forces
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Load matlab file
        mat         = h5s.loadmat(wtModelProp.fname)
        
        # Pressure coeff.
        self.Cp     = mat['Wind_pressure_coefficients']

        # Calculate forces [in kN]
        self.F_p = wtModelProp.A_i * self.Cp * 0.5 * wtModelProp.rho_air  * wtModelProp.uH  ** 2 * 10 ** -3

        # Calculating base forces
        self.calcModelBaseForces(wtModelProp, buildProp)

        # Calculating floor forces
        self.calcModelFloorForces(wtModelProp, buildProp)

    def calcModelBaseForces(self, wtModelProp, buildProp):
        """Calculating base forces from wind tunnel model forces in model scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Sort in forces
        self.BF_p = np.zeros(wtModelProp.nT)
        self.BM_p = np.zeros(wtModelProp.nT)
        
        # shape(F) = (nT, i) -> Loop over columns of array
      
        for i, F_i in enumerate(self.F_p.T, 0):
            # Face 1 -> +DragF, +DragM
            if wtModelProp.face[i] == 1 and buildProp.dn == 'D':
                self.BF_p = self.BF_p + F_i
                self.BM_p = self.BM_p + F_i * wtModelProp.z[i] 
            # Face 2 -> +LiftF, -LiftM
            elif wtModelProp.face[i] == 2 and buildProp.dn == 'L':
                self.BF_p = self.BF_p + F_i
                self.BM_p = self.BM_p - F_i * wtModelProp.z[i] 
            # Face 3 -> -DragF, -DragM
            elif wtModelProp.face[i] == 3 and buildProp.dn == 'D':
                self.BF_p = self.BF_p - F_i
                self.BM_p = self.BM_p - F_i * wtModelProp.z[i] 
            # Face 4 -> -LiftF, +LiftM
            elif wtModelProp.face[i] == 4 and buildProp.dn == 'L':
                self.BF_p = self.BF_p - F_i
                self.BM_p = self.BM_p + F_i * wtModelProp.z[i] 

    def calcModelFloorForces(self, wtModelProp, buildProp):
        """Calculating floor forces from wind tunnel model forces in model scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Sort in drag / lift forces for each floor
        self.LF_p = np.zeros((wtModelProp.nz, wtModelProp.nT))

        # shape(F) = (nT, i) -> Loop over columns of array
        for i, F_i in enumerate(self.F_p.T, 0):
            # get index where to sort to
            j = np.where(wtModelProp.z_lev == wtModelProp.z[i])
            # Face 1 -> +DragF, +DragM
            if wtModelProp.face[i] == 1 and buildProp.dn == 'D':
                self.LF_p[j,:] = self.LF_p[j,:] + F_i
            # Face 2 -> +LiftF, -LiftM
            elif wtModelProp.face[i] == 2 and buildProp.dn == 'L':
                self.LF_p[j,:] = self.LF_p[j,:] + F_i
            # Face 3 -> -DragF, -DragM
            elif wtModelProp.face[i] == 3 and buildProp.dn == 'D':
                self.LF_p[j,:] = self.LF_p[j,:] - F_i
           # Face 4 -> -LiftF, +LiftM
            elif wtModelProp.face[i] == 4 and buildProp.dn == 'L':
                self.LF_p[j,:] = self.LF_p[j,:] - F_i

class buildAeroForces:
    """Class containing the building scale forces of the wind tunnel
    :cvar BF_p: time series of base shear
    :vartype BF_p: np.arr[float]
    :cvar BM_p: time series of base moment
    :vartype BM_p: np.arr[float]
    :cvar LF_p: time series of forces at each level
    :vartype LF_p: np.arr[float, float]
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
        self.LF_p = wtModelAeroForces.LF_p * scalingFactors.lambda_F
        
        # Scale base forces
        self.BF_p = wtModelAeroForces.BF_p * scalingFactors.lambda_F 
        self.BM_p = wtModelAeroForces.BM_p * scalingFactors.lambda_M

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   