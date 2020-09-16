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
# sp..  sample
# f...  frequency

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

# class baseResponse(response):
#     def __init__(self):
#         super.__init__()


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
