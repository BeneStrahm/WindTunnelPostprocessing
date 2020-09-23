# ------------------------------------------------------------------------------
# Description:  
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-22
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  

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

class scalingFactors:
    """Class containing scaling factors from the wind tunnel to the building scale
    :cvar lambda_u: wind velocity scaling
    :vartype lambda_u: float
    :cvar lambda_g: geometric scaling
    :vartype lambda_g: float
    :cvar lambda_fq: frequency scaling
    :vartype lambda_fq: float
    :cvar lambda_t: time scaling
    :vartype lambda_t: float
    :cvar lambda_F: force scaling
    :vartype lambda_F: float
    :cvar lambda_M: moment scaling
    :vartype lambda_M: float
    """
    def __init__(self, wtModelProp, buildProp):
        """Inits the class
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # call time scaling function
        self.calcScaleFactorsTime(wtModelProp, buildProp)
        # call force scaling function
        self.calcScaleFactorsForce(wtModelProp, buildProp)

    def calcScaleFactorsTime(self, wtModelProp, buildProp):
        """Scale time/freq. from model to full scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Scaling factors
        self.lambda_u = buildProp.uH / wtModelProp.uH
        self.lambda_g = buildProp.H / wtModelProp.H
        self.lambda_fq= self.lambda_u / self.lambda_g
        self.lambda_t = 1 / self.lambda_fq

    def calcScaleFactorsForce(self, wtModelProp, buildProp):
        """Scale base forces from model to full scale
        :param wtModelProp: model scale properties of the building
        :type wtModelProp: :class:`~modelProp.wtModelProp`
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Scaling factors
        self.lambda_F = self.lambda_u ** 2 * self.lambda_g ** 2
        self.lambda_M = self.lambda_u ** 2 * self.lambda_g ** 3

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   