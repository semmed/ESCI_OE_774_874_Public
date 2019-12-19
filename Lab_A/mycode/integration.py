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
        n_twtt_times = len(twtt.times)

        #3.2.3
        # the variables determined in the integration
        R_tx = list()
        R_rx = list()
        self.lever_arm_pos_tx=np.zeros([3,n_twtt_times])
        self.lever_arm_pos_rx=np.zeros([3,n_twtt_times])
        self.lever_arm_trans_tx=np.zeros([3,n_twtt_times])
        self.lever_arm_rec_rx=np.zeros([3,n_twtt_times])
        self.pos_rp_tx=np.zeros([3,n_twtt_times])
        self.pos_rp_rx=np.zeros([3,n_twtt_times])
        self.pos_trans_tx=np.zeros([3,n_twtt_times])
        self.pos_rec_rx=np.zeros([3,n_twtt_times])
              
        # 3.2.4
        # Determine the transmit times as posix times
        t_twtt = np.array([e.timestamp() for e in twtt.times])

        # Obtaine the times of the various observed data as posix times
        t_pos = np.array([e.timestamp() for e in pos.times])
        t_mru = np.array([e.timestamp() for e in motions.times])
        t_wl = np.array([e.timestamp() for e in water_levels.times])

        # 3.2.5.1 
        # Interpolate the time series for the TWTT transmit times
        
        self.p_tx = np.interp(t_twtt, t_mru, motions.pitch)
        self.r_tx = np.interp(t_twtt, t_mru, motions.roll)
        self.y_tx = np.interp(t_twtt, t_mru, motions.yaw)
        self.h_tx = np.interp(t_twtt, t_mru, motions.heave)
        self.wl_tx = np.interp(t_twtt, t_wl, water_levels.water_levels)
        
        # Determine the interpolation function for the positions
        f=interp1d(t_pos,pos.proj_pos,bounds_error=False)
        self.pos_proj_ant_tx=f(t_twtt)

        # 3.2.5.2
        # Determine the reception times as posix times
        t_twtt += twtt.twtts

        # 3.2.5.3
        # Interpolate for the time of reception
        self.pos_proj_ant_rx=f(t_twtt)
        self.p_rx = np.interp(t_twtt, t_mru, motions.pitch)
        self.r_rx = np.interp(t_twtt, t_mru, motions.roll)
        self.y_rx = np.interp(t_twtt, t_mru, motions.yaw)
        self.h_rx = np.interp(t_twtt, t_mru, motions.heave)
        self.wl_rx = np.interp(t_twtt, t_wl, water_levels.water_levels)

        # 3.2.6
        # For each ping (twtt) determine a sounding solution
        self.depth = np.zeros((n_twtt_times))
        self.sounding = np.zeros((n_twtt_times))
        self.virtual_txrx = np.zeros((3,n_twtt_times))

        # 3.2.7
        ping = 0
        for t in t_twtt:

            # 3.2.8
            # Calculate the right-handed Euclidean Euler andgle Rotation matrices at the
            # time of transmit around the x,y and z axes
            Rx_tx = np.array([[1, 0,                     0                   ],
                              [0, cos(self.r_tx[ping]), -sin(self.r_tx[ping])],
                              [0, sin(self.r_tx[ping]),  cos(self.r_tx[ping])]])

            Ry_tx = np.array([[cos(self.p_tx[ping]),  0, sin(self.p_tx[ping]) ],
                              [0,                     1,  0                   ],
                              [-sin(self.p_tx[ping]), 0,  cos(self.p_tx[ping])]])

            Rz_tx = np.array([[cos(self.y_tx[ping]), -sin(self.y_tx[ping]), 0],
                              [sin(self.y_tx[ping]),  cos(self.y_tx[ping]), 0],
                              [0,                     0,                    1]])

            # 3.2.8.2
            # Calculate the total rotation matrix at transmit in the order x, y, z
            R_tx.append( Rz_tx@Ry_tx@Rx_tx)

            # 3.2.8.3
            # Rotation matrices at time of reception
            Rx_rx = np.array([[1, 0,                     0                   ],
                              [0, cos(self.r_rx[ping]), -sin(self.r_rx[ping])],
                              [0, sin(self.r_rx[ping]),  cos(self.r_rx[ping])]])

            Ry_rx = np.array([[cos(self.p_rx[ping]),  0,  sin(self.p_rx[ping])],
                              [0,                     1,  0                   ],
                              [-sin(self.p_rx[ping]), 0,  cos(self.p_rx[ping])]])

            Rz_rx = np.array([[cos(self.y_rx[ping]), -sin(self.y_rx[ping]), 0],
                              [sin(self.y_rx[ping]),  cos(self.y_rx[ping]), 0],
                              [0,                     0,                    1]])

            # Calculate the total rotation matrix at receive in the order x, y, z
            R_rx.append( Rz_rx@Ry_rx@Rx_rx)

            # 3.2.9
            # Calculate the georeferenced lever_arms at the transmit times
            # i.e., rotate the vessel lever arms using the rotation matrix R_tx
            self.lever_arm_pos_tx[:,[ping]]=R_tx[ping]@vessel.lever_arm_pos
            self.lever_arm_trans_tx[:,[ping]]=R_tx[ping]@vessel.lever_arm_trans

            # 3.2.10
            # Calculate the geo referenced lever_arms at the reception times using R_rx
            self.lever_arm_pos_rx[:,[ping]]=R_rx[ping]@vessel.lever_arm_pos
            self.lever_arm_rec_rx[:,[ping]]=R_rx[ping]@vessel.lever_arm_rec

            # 3.2.11
            # Calculate depth as height relative to the ships transducer
            # Note that this requires both the waterline and the lever_arm to get the right starting point 
            # in the water column - not heave - as we assume that the profile moves up and down with the surface.
            self.depth[ping] =  np.mean(self.ssp.obs_ss) / 2 * twtt.twtts[ping]

            
            # 3.2.12
            # Now determine the rp and transducer position at the transmits
            # Be careful about axis alignment!
            self.pos_rp_tx[0,[ping]]=self.pos_proj_ant_tx[0,[ping]]-self.lever_arm_pos_tx[1,[ping]]
            self.pos_rp_tx[1,[ping]]=self.pos_proj_ant_tx[1,[ping]]-self.lever_arm_pos_tx[0,[ping]]
            self.pos_rp_tx[2,[ping]]=self.pos_proj_ant_tx[2,[ping]]-self.lever_arm_pos_tx[2,[ping]]
            self.pos_trans_tx[:,[ping]]=self.pos_rp_tx[:,[ping]]+self.lever_arm_trans_tx[:,[ping]]

            # Now determine the rp and transducer position at reception
            self.pos_rp_rx[0,[ping]]=self.pos_proj_ant_rx[0,[ping]]-self.lever_arm_pos_rx[1,[ping]]
            self.pos_rp_rx[1,[ping]]=self.pos_proj_ant_rx[1,[ping]]-self.lever_arm_pos_rx[0,[ping]]
            self.pos_rp_rx[2,[ping]]=self.pos_proj_ant_rx[2,[ping]]-self.lever_arm_pos_rx[2,[ping]]
            self.pos_rec_rx[:,[ping]]=self.pos_rp_rx[:,[ping]]+self.lever_arm_rec_rx[:,[ping]]
            
            # 3.2.13
            # Calculate the location of the virtual transucer located at the mean position
            # of the transmitter and receiver
            self.virtual_txrx[:,[ping]] = (self.pos_trans_tx[:,[ping]]+self.pos_rec_rx[:,[ping]])/2
            
            # 3.2.14 
            
            self.sounding[ping] = self.depth[ping] - self.virtual_txrx[2,[ping]]


            # update the ping count
            ping += 1
        

    def draw(self):
        fig=plt.figure(figsize=(12, 6))
        ax1 = plt.subplot(3,1,1)
        plt.plot(self.twtt.times, self.depth)
        plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2)
        plt.plot(self.twtt.times[0:-100], self.sounding[0:-100])

#     def draw(self):
#         fig=plt.figure(figsize=(12, 6))
#         ax1 = plt.subplot(3,1,1)
#         plt.plot(self.twtt.times, self.twtt.twtts*np.mean(self.ssp.obs_ss)/2)
# #         plt.plot(self.twtt.times, self.twtt.twtts *
# #              np.mean(self.ssp.obs_ss)/2+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:])
# #         plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:])

#         plt.title('Depths [m]')
#         plt.ylabel('Depths [m] →')
#         plt.xlabel('Time ('+self.twtt.metadata['time_basis']+') →')
#         ax1.invert_yaxis()
    
#         ax2 = plt.subplot(3,1,2)
# #         plt.plot(self.twtt.times, (self.twtt.twtts * \
# #          np.mean(self.ssp.obs_ss)/2+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:])-(self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:]))
    
#         plt.title('Depths [m]')
#         plt.ylabel('Depths [m] →')
#         plt.xlabel('Time ('+self.twtt.metadata['time_basis']+') →')
#         ax2.invert_yaxis()
        
#         ax3 = plt.subplot(3,1,3)
#         twtt_times=self.twtt.times[0:100]
#         twtts=self.twtt.twtts[0:100]
#         obs_ss=self.ssp.obs_ss[0:100]
#         h_tx=self.h_tx[0:100]
#         h_rx=self.h_rx[0:100]
# #         la=self.la_trans_rec_txrx[2,0:100]
#         depth=self.depth[0:100]
#         ant_tx_full=self.pos_proj_ant_rx[:,0:100]
#         rp_tx_full=self.pos_rp_tx[:,0:100]
#         trans_tx_full=self.pos_trans_tx[:,0:100]
        
# #         plt.plot(twtt_times, (twtts * \
# #                   np.mean(obs_ss)/2+(h_tx+h_rx)/2+la)- \
# #                    (depth+(h_tx+h_rx)/2+la))
    
#         plt.plot(twtt_times,(h_tx+h_rx)/300+.125)
#         plt.title('Depths [m]')
#         plt.ylabel('Depths [m] →')
#         plt.xlabel('Time ('+self.twtt.metadata['time_basis']+') →')
#         ax3.invert_yaxis()
        
#         # Space out the plots so that they do not overlap
#         plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=.5, wspace=0.4)
#         plt.show()        
        
#         # Plot the navigation
#         fig=plt.figure(figsize=(12, 6))
#         ax4 = plt.subplot(2,1,1)
#         plt.plot(self.pos_proj_ant_rx[0,:],self.pos_proj_ant_rx[1,:],'b.',label='Antenna')
#         plt.plot(self.pos_rp_tx[0,:],self.pos_rp_tx[1,:],'r.',label='RP')
#         plt.plot(self.pos_trans_tx[0,:],self.pos_trans_tx[1,:],'k.',label='Transducer')
#         plt.legend()
#         plt.axis('equal')
#         ax5 = plt.subplot(2,1,2)
#         plt.plot(ant_tx_full[0,0:100],ant_tx_full[1,0:100],'b.',label='Antenna')
#         plt.plot(rp_tx_full[0,0:100],rp_tx_full[1,0:100],'r.',label='RP')
#         plt.plot(trans_tx_full[0,0:100],trans_tx_full[1,0:100],'k.',label='Transducer')
#         plt.axis('equal')
#         plt.legend()
#         plt.show()
    
#     def draw_depths(self):
#         fig=plt.figure(figsize=(12, 6))
#         plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx[2,:])

#         plt.title('Depths [m]')
#         plt.ylabel('Depths [m] →')
#         plt.xlabel('Time ('+self.twtt.metadata['time_basis']+') →')
#         plt.gca().invert_yaxis()
    
            
    def heave_gnss():
        # Determine heave from the rp trajectorty
        # Use an 30 second filtering window
         return heave
        
    def qc(self):
        # The match of the IMU heave to the positioning z-component is one quality control indicator
        
        fig=plt.figure(figsize=(12, 6))
        ax1 = plt.subplot(3,1,1)
#         plt.plot(self.twtt.times, self.pos_ant
#         plt.plot(self.twtt.times, self.depth+(self.h_tx+self.h_rx)/2+self.la_trans_rec_txrx)
        