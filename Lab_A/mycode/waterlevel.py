import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt

class WaterLevel:
    """A Class for handling Water Level Data"""

    def __init__(self):

        # The data attributes
        self.times = list()
        self.water_levels = list()
        self.metadata = dict()
        self.data_path=str()
        self.metadata["units"] = "m"
        self.metadata["datum_type"] = None
        self.metadata["datum_name"] = None
        self.metadata["time_basis"] = "UTC"        
        self.metadata["location_name"] = "Unknown"


    def __str__(self): 
        txt  = "Location name          : %s\n" % (self.metadata["location_name"])
        txt += "Reference Surface Type : %s\n" % (self.metadata["datum_type"])
        txt += "Reference Surafce Name : %s\n" % (self.metadata["datum_name"])
        txt += "Observation Time Basis : %s\n" % (self.metadata["time_basis"])
        txt += "Observations Units     : %s\n" % (self.metadata["units"])

        if len(self.times):
            txt += "Start Time             : %s\n" % (min(self.times))
        else:
            txt += "No time data present\n"
            
        if len(self.water_levels):
            txt += "Minimum Water Level    : %.2f%s\n" % (min(self.water_levels), self.metadata["units"])           
            txt += "Maximum Water Level    : %.2f%s\n" % (max(self.water_levels), self.metadata["units"])       
        else:
            txt += "No water level data present\n"
        return txt
    
    def read_jhc_file(self, fullpath):

        pass # replace this with code that causes the object to read a file
         
    def draw(self):
        pass # replace this with code that causes the object to draw its contents