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
class wtModelProp:
    """Class containing the model scale properties of the building
    :cvar fname: full path to file w/ wind tunnel data
    :vartype fname: string
    :cvar H: model-scale building height [m]
    :vartype H: float
    :cvar Loc: location of measured points
    :vartype Loc: list [float, float, float, float]
    :cvar x: x-coordinate of measured points [m]
    :vartype x: list[float]
    :cvar z: z-coordinate of measured points [m]
    :vartype z: list[float]
    :cvar n: number of measured points
    :vartype n: list[float]
    :cvar face: face of measured points
    :vartype face: list[float]
    :cvar A_i: influence area of measured points [m2]
    :vartype A_i: float
    :cvar z_lev: z-coordinate of measurement levels [m]
    :vartype z_lev: list[float]
    :cvar nz: number of measurement levels
    :vartype nz: int
    :cvar uH: model scale wind speed at top of the building [m/s]
    :vartype uH: float
    :cvar fq_sp: measurement frequency [Hz]
    :vartype fq_sp: float
    :cvar dT: measurement time stepping [s]
    :vartype dT: float
    :cvar T: sample period [s]
    :vartype T: string
    :cvar nT: number of time steps
    :vartype nT: string
    """
    def __init__(self, fname):
        """Inits the class
        :param fname: full path to file w/ wind tunnel data
        :type fname: string
        """
        self.fname  = fname

    def loadTPUModelProp(self):
        """Loading wind tunnel model characteristics
        """        
        # Load matlab file
        mat     = h5s.loadmat(self.fname)

        # Geometrie
        self.H  = float(mat['Building_height'][0])

        # Measurement locations
        self.Loc    = mat['Location_of_measured_points']

        # Coordinates [in m], number and face of measurement points
        self.x      = self.Loc[0]
        self.z      = self.Loc[1]
        self.n      = self.Loc[2]
        self.face   = self.Loc[3]

        # Measurement area [in m2]
        self.A_i    = 0.02 * 0.02

        # Air density
        self.rho_air= 1.25

        # Get levels & sort them
        self.z_lev  = np.unique(self.z)
        self.z_lev  = self.z_lev[::-1]
        self.nz     = len(self.z_lev)

        # Wind speeds
        self.uH     = float(mat['Uh_AverageWindSpeed'][0])

        # Time scales
        self.fq_sp  = float(mat['Sample_frequency'][0])
        self.dT     = 1 / self.fq_sp
        self.T      = float(mat['Sample_period'][0])
        self.nT     = int(self.T * self.fq_sp) 

class buildProp():
    """Class containing the full scale properties of the building
    :cvar H: full-scale building height [m]
    :vartype H: float
    :cvar E: E-Modulus [kN/m²]
    :vartype E: float
    :cvar I: moment of inertia in drag direction [m4]
    :vartype I: float
    :cvar mue: mass distribution [t/m]
    :vartype mue: float
    :cvar D: damping [%]
    :vartype D: float
    :cvar uH: full scale wind speed at top of the building [m/s]
    :vartype uH: float
    :cvar x: x-coordinate of measured points [m]
    :vartype x: list[float]
    :cvar z: z-coordinate of measured points [m]
    :vartype z: list[float]
    :cvar n: number of measured points
    :vartype n: list[float]
    :cvar face: face of measured points
    :vartype face: list[float]
    :cvar A_i: influence area of measured points [m2]
    :vartype A_i: float
    :cvar z_lev: z-coordinate of measurement levels [m]
    :vartype z_lev: list[float]
    :cvar nz: number of measurement levels
    :vartype nz: int
    :cvar fq_sp: measurement frequency [Hz]
    :vartype fq_sp: float
    :cvar dT: measurement time stepping [s]
    :vartype dT: float
    :cvar T: sample period [s]
    :vartype T: string
    :cvar nT: number of time steps
    :vartype nT: string
    """
    def __init__(self, H, E, I, mue, D, uH):
        """Inits the class.
        :param H: full-scale building height [m]
        :type H: float
        :param E: E-Modulus [kN/m²]
        :type E: float
        :param I: moment of inertia in drag direction [m4]
        :type I: float
        :param mue: mass distribution [t/m]
        :type mue: float
        :param D: Damping [%]
        :type D: float
        :param uH: full scale wind speed at top of the building [m/s]
        :type uH: float
        """
        # Geometric properties
        self.H      = H
        # Static properties
        self.I    = I
        self.E      = E

        # Dynamic properties
        self.mue    = mue
        self.D      = D

        # Atmospheric properties
        self.uH     = uH

    def scaleBuildProp(self, wtModelProp, scalingFactors):
        """scale wind tunnel model data to full scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param scalingFactors: scaling factors from wind tunnel to building scale
        :type scalingFactors: :class:`~scaling.scalingFactors`
        """
        # Coordinates [in m], number and face of measurement points
        self.x      = wtModelProp.x * scalingFactors.lambda_g
        self.z      = wtModelProp.z * scalingFactors.lambda_g

        # Measurement area [in m2]
        self.A_i    = wtModelProp.A_i * scalingFactors.lambda_g ** 2

        # Air density
        self.rho_air= 1.25

        # Get levels & sort them
        self.z_lev  = wtModelProp.z_lev * scalingFactors.lambda_g
        self.nz     = len(self.z_lev)

        # Time scales
        self.fq_sp  = wtModelProp.fq_sp * scalingFactors.lambda_fq
        self.dT     = 1 / self.fq_sp
        self.T      = wtModelProp.T * scalingFactors.lambda_t
        self.nT     = int(self.T * self.fq_sp) 

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   

def test():
    a = wtModelProp("T115_4_000_1.mat")
    a.loadTPUModelProp()