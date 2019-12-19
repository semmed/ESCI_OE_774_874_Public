import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from numpy import pi, cos, sin, log, exp
import numpy as np

from mycode.baseline import Baseline
from mycode.position import Position
class Network:
    """A Class for handling Baseline Data"""

    def __init__(self,dx,dy,dz):

        # The data attributes
        self.bl = np.matrix([[dx],[dy],[dz]])

        
        
        