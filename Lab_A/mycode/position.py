import os
from datetime import datetime, timezone, timedelta
import matplotlib.pyplot as plt
import numpy as np


## Created by Semme J. Dijkstra 8/28/2019
    #  Added method: read_hypack_raw_file, Semme J. Dijkstra 10/4/2019

class Position:
    """A Class for handling Position Data"""

    def __init__(self):
       
        # The data attributes
        self.times = list()
        
        # The geodetic coordinates - these are curvilinear so do not put them
        # in vectors, as that is a linear concept
        self.latitudes = list()
        self.longitudes = list()
        self.heights = list()
        self.qualities = list()
        self.n_sats = list()
        self.hdops = list()
        self.separations = list()
        self.corr_ages = list()
        self.corr_stations = list()
        
        # The projected coordinates are (usually) in Cartesian coordinates (which are linear), store them in numpy arrays
        self.proj_pos=np.array([])
        self.data_path=str()
        self.metadata = dict()
        self.metadata["geodetic_units"] = "rad"
        self.metadata["height_units"] = "m"
        self.metadata["proj_units"] = "m"
        self.metadata["geoid_name"] = None
        self.metadata["ellipsoid_name"] = None
        self.metadata["height_relative_to"] = None
        self.metadata["time_basis"] = "UTC"
        self.metadata["proj_str"] = None


    def __str__(self): 
        txt  = "Geoid Used             : %s\n" % (self.metadata["location_name"])
        txt += "Ellipsoid Used         : %s\n" % (self.metadata["datum_type"])
        txt += "Reference Surface Name : %s\n" % (self.metadata["datum_name"])
        txt += "Observation Time Basis : %s\n" % (self.metadata["time_basis"])
        txt += "Observations Units     : %s\n" % (self.metadata["units"])

        if len(self.times):
            txt += "Start Time             : %s\n" % (min(self.times))
        else:
            txt += "No time data present\n"
        
        # If there are no latitudes, logically there are no longitudes either
        
        if len(self.latitudes):
            txt += "Minimum latitude       : %.2f\n" % (min(self.latitudes))           
            txt += "Maximum latitude       : %.2f\n" % (max(self.latitudes))       
            txt += "Minimum longitude      : %.2f\n" % (min(self.longitude))           
            txt += "Maximum longitude      : %.2f\n" % (max(self.longitude))       
        else:
            txt += "No geodetic coordinates are present\n"
            
        # Test heights separately, as we may just have horizontal positioning
            
#        if len(self.heights):
#             txt += "Minimum height       : %.2f%s\n" % (min(self.heights), self.metadata["height_units"]         
#             txt += "Minimum height       : %.2f%s\n" % (min(self.heights), self.metadata["height_units"]) 

        return txt

    def read_jhc_file(self, fullpath):

        pass # replace this with code that causes the object to read a data file