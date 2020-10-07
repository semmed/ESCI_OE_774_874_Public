import os
from datetime import datetime, timezone
from numpy import pi, cos, sin, log, exp, arccos, tan, arctan, tanh, arctanh
import numpy as np
from mycode.position import *

class SSP:
    """A Class for handling Sound Speed Profile data"""

    def __init__(self):

        # The data attributes
        self.obs_time = None
        self.log_time = None
        self.obs_latitude = None
        self.obs_longitude = None
        self.vessel_latitude = None
        self.vessel_longitude = None
        self.obs_sample = list()
        self.obs_depth = list()
        self.obs_ss = list()
        self.proc_depth = np.array([])
        self.proc_ss = np.array([])
        self.g = np.array([])
        self.proc_ss = np.array([])
        self.twtt_layer=np.array([])
        self.vessel_speed = None
        self.bot_depth = None

        self.metadata = dict()
        self.metadata["units"] = "rad"
        self.metadata["count"] = None
        self.metadata["geoid"] = None
        self.metadata["ellipsoid"] = None
        self.metadata["chart_datum"] = None
        self.metadata["time_basis"] = "UTC"
        self.metadata["name"] = None


    def read_jhc_file(self, fullpath):
        # Check to see whether data already exists in the object
        
        if self.obs_depth:
            raise RuntimeError('SSP object already contains a profile')
            
        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        svp_file = open(fullpath)
        svp_content = svp_file.read()
        svp_file.close

        # Tokenize the contents
        svp_lines = svp_content.splitlines()
        self.obs_time = datetime.fromtimestamp(
            float(svp_lines[1].split()[0]), timezone.utc)
        self.log_time = datetime.fromtimestamp(
            float(svp_lines[2].split()[0]), timezone.utc)
        self.obs_latitude = float(svp_lines[3].split()[0])
        self.obs_longitude = float(svp_lines[3].split()[1])
        self.vessel_latitude = float(svp_lines[4].split()[0])
        self.vessel_longitude = float(svp_lines[4].split()[1])
        self.metadata["count"] = int(svp_lines[5].split()[0])

        count = 0  # initialize the counter for the number of rows read

        for svp_line in svp_lines[16:]:
            observations = svp_line.split()  # Tokenize the stringS
            self.obs_sample.append(float(observations[0]))
            self.obs_depth.append(float(observations[1]))
            self.obs_ss.append(float(observations[2]))
            count += 1

        if self.metadata["count"] != count:
            raise RuntimeError('Nr of Samples read ('+str(count) +
                               ') does not match metadata count (' +
                               str(self.metadata["count"])+')')




    def draw(self):
        pass # replace this with code that causes the object to draw its contents

