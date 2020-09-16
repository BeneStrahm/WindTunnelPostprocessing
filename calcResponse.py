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
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------
from helpers.pyExtras import getKeyList
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

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class response(object):
    def __init__(self, modelForces, uH_fs, H_fs):
        """Scale time/freq. from model to full scale
        :param modelForce: obj w/ wind tunnel measurements
        :param uH_fs, H_fs: flt w/ full scale building properties
        """
        self.modelForces = modelForces
        self.H_fs   = H_fs
        self.uH_fs  = uH_fs
    
    def scaleTime(self):
        """Scale time/freq. from model to full scale
        :param uH_f, H_f: float w/ full scale building properties
        """
        # Model properties
        self.modelForces.dT_ms    = 1 / self.modelForces.fs

        # Scaling factors
        lambda_u = self.uH_fs / uH_m
        lambda_g = self.H_fs / H_m 
        lambda_f = lambda_u / lambda_g
        lambda_t = 1 / lambda_f

        # Scale quantities
        dT_f    = lambda_t * dT_m
        f_f     = lambda_f * f_m

        return dT_f, f_f, nT

class baseResponse(response):
    def __init__(self):
        super.__init__()


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
