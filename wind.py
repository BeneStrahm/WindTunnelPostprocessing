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
# fq...  frequency
# dn... direction
# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class windStats(object):
    def __init__(self, uH_50):
        """Inits the class
        :param uH_50: flt w/ wind speed at z = H
        """        
        self.calcWindStatsWithNorm(uH_50)

        # self.calcWindStatsRange(uH_50)

    def calcWindStatsWithNorm(self, uH_50):
        """Calculating wind speeds at different return periods, 
        setting initial wind speed for return period R = 50y 
        acc. to DIN EN 1991-1-4 / NA Abs. 4.2
        :param uH_50: flt w/ wind speed at z = H
        """   
        # Define desired return periods
        self.uH  = {#"uH_001":0, \
                        "uH_002":0, \
                        "uH_005":0, \
                        "uH_010":0, \
                        "uH_020":0, \
                        "uH_050":uH_50, \
                        "uH_100":0, \
                        "uH_200":0, \
                        "uH_500":0}     

        for key in getKeyList(self.uH):
            # Get return period / probability of excedence
            R = float(key[-3:])
            p = 1 / R

            # According to DIN EN 1991-1-4 / NA Abs. 4.2, Anmerkung 5
            K = 0.1
            n = 1

            # Get probability factor
            cprob   = ((1 - K * np.log(-np.log(1-p))) / (1 - K * np.log(-np.log(0.98)))) ** n

            # Calculate wind speeds at different return periods
            self.uH[key] = cprob * self.uH["uH_050"]  

    def calcWindStatsRange(self, uH_50):
        """Calculating a range of wind speeds 
        """   
        u_min   = uH_50 * 0.5
        u_max   = uH_50 * 2.0
        nu      = 40

        # Define desired return periods
        self.uH  = {}

        u_list = np.linspace(u_min, u_max, nu)
        # Define desired return periods
        for i, u in enumerate(u_list):
            key = 'uH_'+'{:03d}'.format(i)
            item = u
            self.uH[key] = item

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
