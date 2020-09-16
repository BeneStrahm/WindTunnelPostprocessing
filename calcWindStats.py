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

class windStats(object):
    def __init__(self, uH_fs_50):
        """Setting initial wind speed for return period R = 50y 
        :param u_fs_h_50: flt w/ wind speed at z = H
        """        
        # Define desired return periods
        self.uH_fs  = {#"uH_fs_001":0, \
                        "uH_fs_002":0, \
                        "uH_fs_005":0, \
                        "uH_fs_010":0, \
                        "uH_fs_020":0, \
                        "uH_fs_050":uH_fs_50, \
                        "uH_fs_100":0}

    def calcWindStats(self):
        """Calculating wind speeds at different return periods
        acc. to DIN EN 1991-1-4 / NA Abs. 4.2
        """        

        for key in getKeyList(self.uH_fs):
            # Get return period / probability of excedence
            R = float(key[-3:])
            p = 1 / R

            # According to DIN EN 1991-1-4 / NA Abs. 4.2, Anmerkung 5
            K = 0.1
            n = 1

            # Get probability factor
            cprob   = ((1 - K * np.log(-np.log(1-p))) / (1 - K * np.log(-np.log(0.98)))) ** n

            # Calculate wind speeds at different return periods
            self.uH_fs[key] = cprob * self.uH_fs["uH_fs_050"]  

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
