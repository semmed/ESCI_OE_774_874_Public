import numpy as np
from scipy.interpolate import interp1d
from numpy import pi, cos, sin, log, exp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mycode.EchosounderData import EchosounderData
from mycode.Motion import Motion
from mycode.SSP import SSP
from mycode.waterlevel import WaterLevel
from mycode.vessel import Vessel
from mycode.position import Position

class Integration:
    """A Class for Integrating Data to Create Soundings"""

    def __init__(self, twtt, pos, motions, sound_speed_profile, water_levels, vessel):
        
        # 3.2.1
        # For now we can only integrate if the Positions have been projected - we will use 
        # either UTM or UPS depending on the latitude
        
        self.vessel = vessel
        self.twtt = twtt
        self.ssp = sound_speed_profile

        # 3.2.2
        # The number of twtts
        

        #3.2.3
        # the variables determined in the integration
        R_tx = list()
        R_rx = list()
        self.lever_arm_pos_tx=np.zeros([3,n_twtt_times])
        
              
        # 3.2.4
        # Determine the transmit times as posix times
        

        # Obtain the times of the various observed data as posix times
        

        # 3.2.5.1 
        # Interpolate the time series for the TWTT transmit times
        
        
        
        # Determine the interpolation function for the positions
        

        # 3.2.5.2
        # Determine the reception times as posix times
        

        # 3.2.5.3
        # Interpolate for the time of reception
        

        # 3.2.6
        # For each ping (twtt) determine a sounding solution
        

        # 3.2.7
        ping = 0
        for t in t_twtt:

            # 3.2.8
            # Calculate the right-handed Euclidean Euler andgle Rotation matrices at the
            # time of transmit around the x,y and z axes
            Rx_tx = np.array([[1, 0,                     0                   ],
                              [0, cos(self.r_tx[ping]), -sin(self.r_tx[ping])],
                              [0, sin(self.r_tx[ping]),  cos(self.r_tx[ping])]])

            

            # 3.2.8.2
            # Calculate the total rotation matrix at transmit in the order x, y, z
            

            # 3.2.8.3
            # Rotation matrices at time of reception
            

            # Calculate the total rotation matrix at receive in the order x, y, z
            

            # 3.2.9
            # Calculate the georeferenced lever_arms at the transmit times
            # i.e., rotate the vessel lever arms using the rotation matrix R_tx

            # 3.2.10
            # Calculate the geo referenced lever_arms at the reception times using R_rx

            # 3.2.11
            # Calculate depth as height relative to the ships transducer
            # Note that this requires both the waterline and the lever_arm to get the right starting point 
            # in the water column - not heave - as we assume that the profile moves up and down with the surface.
            
            
            # 3.2.12
            # Now determine the rp and transducer position at the transmits
            # Be careful about axis alignment!

            # Now determine the rp and transducer position at reception
            
            # 3.2.13
            # Calculate the location of the virtual transucer located at the mean position
            # of the transmitter and receiver
            
            # 3.2.14 Calculate the sounding, that is the
            
            self.sounding[ping] = self.virtual_txrx[2,[ping]]-self.depth[ping]


            # update the ping count
            ping += 1
        

    def draw(self):
        pass