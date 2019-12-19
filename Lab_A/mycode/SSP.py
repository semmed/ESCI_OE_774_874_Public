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

    # The I/O methods:
    
    def read_mvp_file(self, fullpath):
        
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
        
        # Save the file name as meta data
        
        self.metadata["name"] = os.path.basename(fullpath)

        # Tokenize the contents
        svp_lines = svp_content.splitlines()
                
        # Create a Position object
        pos = Position()
        
        n_lines = 0
        for line in svp_lines:
            n_lines += 1

            # Parse the time
            if line[0:8] == "GPS Time":
                # Extract the ZDA record
                obs = line.split()

                # Extract the UTC time string
                self.obs_time = ParseNMEA0183_ZDA(obs[2])

            # Parse the position
            if line[0:12] == "GPS Position":
                # Extract the ZDA record
                obs = line.split()

                # Extract the UTC time string
                _, self.obs_latitude, self.obs_longitude, _, _, _, _, _, _, _ = ParseNMEA0183_GGA(obs[2])

            # Parse the Depth
            if line[0:13] == "Bottom Depth:":
                obs = line.split()
                self.bot_depth=float(obs[2])
        
            # Parse the Vessel Speed
            if line[0:11] == "Ship Speed:":
                obs = line.split()
                self.vessel_speed=float(obs[2])/1.852
            
            if line == '':
                break
    
        # Split the record types line
        rec_type = svp_lines[n_lines].split(',')
    
        # Find the index of the Dpth(m) records
        index_depth = rec_type.index('Dpth(m)')
        
        # Find the index of the Sound Speed records
        index_ss = rec_type.index('SV(m/s)')
        
        for line in svp_lines[ n_lines + 1:]:
            obs = line.split(',')
            self.obs_depth.append( float(obs[index_depth]))
            self.obs_ss.append( float(obs[index_ss]))

        # Make sure that there are no depth reversals due to heaving

        temp = sorted(zip(self.obs_depth, self.obs_ss), key=lambda x: x[0])
        self.obs_depth, self.obs_ss = map(list, zip(*temp))

        # Remove any duplicate depths with associated sound speeds
        d_p=self.obs_depth[0]
        index = 0
        unwanted = []
        for d in self.obs_depth[1:]:
            index += 1
            if d  == d_p:
                unwanted.append(index)
            d_p = d
            
        for e in sorted( unwanted, reverse = True):
            del self.obs_depth[ e]
            del self.obs_ss[ e]
        
        # Add Numpy arrays to work with the data
        self.d = np.array(self.obs_depth)
        self.c = np.array(self.obs_ss)
        
        # Extend the profiles to the surface
        
        if self.d[0] > 0:
            self.d = np.insert(self.d,0,0)
            self.c = np.insert(self.c,0,self.c[0])
   
        # Calculate the gradients
        
        self.g = (self.c[1:] - self.c[0:-1])/(self.d[1:] - self.d[0:-1])

        # Extend the profiles to full ocean depth
        
        if self.d[-1] < 12000:
            self.d = np.append(self.d,12000)
            self.c = np.append(self.c,self.c[-1]+0.017*(12000-self.d[-2]))
            self.g = np.append(self.g,0.017)
           
        # To avoid dividing by zero assign layers with gradient g = 0 a very small value
        
        self.g[self.g == 0] = 0.0001
        
        # Recalculate the sound speed profile
        
        for i in range(1,len(self.g)):
            self.c[i]=self.c[i-1]+(self.d[i] - self.d[i-1])*self.g[i-1]
        
    def determine_twtt(self, d_start, d_end):
        depth=0
        rad_dist=0
        layer_s=0
        layer_e=0 
        
     
        # return 
        
        return depth, rad_dist, layer_s, layer_e;
    
    # B1.3.0.0 Function Declaration
                       
    def determine_depth(self, d_start, th_start, ss_start, twtt):
        
        # B1.3.0.0 Initialization
        depth=0
        rad_dist=0
        layer_s=0
        layer_e=0
        
        swap=False
        if th_start > pi/2:
            swap = True
            th_start = pi - th_start
        
        # B.1.3.0.1 Determine the start layer
        layer_s = sum( d_start >= self.d) - 1
        
        # B.1.3.0.2 Determine the Ray Constant
        ray_c = cos(th_start)/ss_start
           
        # B.1.3.0.3 Calculate Ray Path Properties for Each Layer  
        delta_d = (self.d[1:] - self.d[0:-1])
        r_curve = -1/(self.g[0:]*ray_c)
        th = arccos(self.c[0:]*ray_c)

        dx = r_curve * (sin(th[1:]) - sin(th[0:-1]))
        dz = delta_d     
        hm = 2/self.g*log(self.c[1:]/self.c[0:-1])
        dt = hm+2/self.g*log((1+sin(th[0:-1]))/(1+sin(th[1:])))
        
        # B.1.3.0.4 Determine properties for start layer from the top to the start
        dx_init = r_curve[layer_s]*(sin(th_start)-sin(th[layer_s]))
        dz_init = d_start - self.d[layer_s]
        dt_init = 2/self.g[layer_s]*log(ss_start/self.c[layer_s])
        dt_init += 2/self.g[layer_s]*log((1+sin(th[layer_s]))/(1+sin(th_start)))
        
        # B.1.3.0.5 Accumulate From the Start Layer
        sum_dx = np.cumsum(dx[layer_s:])
        sum_dt = np.cumsum(dt[layer_s:])
        sum_dz = np.cumsum(dz[layer_s:])

        # B.1.3.0.6 Offset to the Start Depth
        sum_dx -= dx_init
        sum_dz -= dz_init
        sum_dt -= dt_init
        
        # B.1.3.0.7 Determine the Number of Boundaries Crossed and the End Layer Index        
        n_bounds =  sum( twtt >= sum_dt[:-1])
        layer_e = n_bounds + layer_s 
        
        # B.1.3.0.8 Determine Properties for the final layer from the top to the end
        t = twtt-sum_dt[n_bounds-1]
        
        if n_bounds >= 0:
            t = twtt-sum_dt[n_bounds-1]
        else:
            raise RuntimeError('SSP start depth in same layer as reflector - not yet implemented')
            
        th_end = 2*arctan(tanh(-t*self.g[layer_e]/4+arctanh(tan(th[layer_e]/2))))
        dx_end = r_curve[layer_e] * (sin(th_end)-sin(th[layer_e]))
        dz_end = -r_curve[layer_e] * (cos(th_end)-cos(th[layer_e]))

        # B.1.3.0.9 The Final Results
        if n_bounds >= 0:
            depth = sum_dz[n_bounds-1] + dz_end + d_start
            rad_dist = sum_dx[n_bounds-1] + dx_end
            
        if swap:
            rad_dist = - rad_dist
        
        #  B1.3.0.0 Return values as tuple

        return depth, rad_dist, layer_s, layer_e; 
    
    def determine_c( self, d):
        layer = sum( d >= self.d) - 1
        ss = self.c[layer]+(d-self.d[layer])*self.g[layer]
        return ss
    
    def determine_twtt(self, d_start, th_start, ss_start, depth):
        
        # B1.3.0.0 Initialization
        dist_z= depth - d_start
        rad_dist=0
        layer_s=0
        layer_e=0
        
#         swap=False
#         if th_start > pi/2:
#             Swap = True
#             th_start = pi - th_start
        
        # B.1.3.0.1 Determine the start layer

        layer_s = sum( d_start >= self.d) - 1
        
        # B.1.3.0.2 Determine the Ray Constant
        ray_c = cos(th_start)/ss_start
                
        # B.1.3.0.3 Calculate Ray Path Properties for Each Layer
        
        delta_d = (self.d[1:] - self.d[0:-1])
        r_curve = -1/(self.g[0:]*ray_c)
        th = arccos(self.c[0:]*ray_c)
        dx = r_curve * (sin(th[1:]) - sin(th[0:-1]))
        dz = delta_d     
        hm = 2/self.g*log(self.c[1:]/self.c[0:-1])
        dt = hm+2/self.g*log((1+sin(th[0:-1]))/(1+sin(th[1:])))
                             
        # B.1.3.0.4 Determine properties for start layer from the top to the start
        
        dx_init = r_curve[layer_s]*(sin(th_start)-sin(th[layer_s]))
        dz_init = d_start - self.d[layer_s]
        dt_init = 2/self.g[layer_s]*log(ss_start/self.c[layer_s])
        dt_init += 2/self.g[layer_s]*log((1+sin(th[layer_s]))/(1+sin(th_start)))
        
        # B.1.3.0.5 Accumulate From the Start Layer
        
        sum_dx = np.cumsum(dx[layer_s:])
        sum_dt = np.cumsum(dt[layer_s:])
        sum_dz = np.cumsum(dz[layer_s:])

        # B.1.3.0.6 Offset to the Start Depth
        sum_dx -= dx_init
        sum_dz -= dz_init
        sum_dt -= dt_init
        
        # B.1.3.1.0 Determine the Number of Boundaries Crossed and the End Layer Index

        layer_e = sum( depth >= np.cumsum(dz))
        n_bounds = layer_e - layer_s 
        
        # B.1.3.1.1 Determine the Vertical Distance Traversed in the Last Layer
        d_z = depth - self.d[layer_e]

        # B.1.3.1.2 Determine the Sound Speed at the End
        c_end = self.c[layer_e]+self.g[layer_e]*d_z

        # B.1.3.1.3 Determine the Depression Angle th_end At the End
        
        th_end = arccos(ray_c*c_end)


        # B.1.3.1.4 Determine the final TWTT and dx
        
        d_t = (d_z/sin(th[layer_e]))/self.c[layer_e]*2
        dx_end = r_curve[layer_e] * (sin(th_end)-sin(th[layer_e]))
        
        if n_bounds >= 0:
            twtt = sum_dt[n_bounds-1]+d_t
            rad_dist = sum_dx[n_bounds-1] + dx_end           

#         if swap:
#             rad_dist = - rad_dist
            
        #  B1.3.0.0 Return values as tuple
        
        return twtt, rad_dist, layer_s, layer_e;               



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
        svp_content = motion_file.read()
        svp_file.close

        # Tokenize the contents
        motion_lines = motion_content.splitlines()
        self.obs_time = datetime.fromtimestamp(
            float(motion_lines[1].split()[0]), timezone.utc)
        self.log_time = datetime.fromtimestamp(
            float(motion_lines[2].split()[0]), timezone.utc)
        self.obs_latitude = float(motion_lines[3].split()[0])
        self.obs_longitude = float(motion_lines[3].split()[1])
        self.vessel_latitude = float(motion_lines[4].split()[0])
        self.vessel_longitude = float(motion_lines[4].split()[1])
        self.metadata["count"] = int(motion_lines[5].split()[0])

        count = 0  # initialize the counter for the number of rows read

        for motion_line in motion_lines[16:]:
            observations = motion_line.split()  # Tokenize the stringS
            self.obs_sample.append(float(observations[0]))
            self.obs_depth.append(float(observations[1]))
            self.obs_ss.append(float(observations[2]))
            count += 1

        if self.metadata["count"] != count:
            raise RuntimeError('Nr of Samples read ('+str(count) +
                               ') does not match metadata count (' +
                               str(self.metadata["count"])+')')

        # Process the data - in the jhc data files this is already a one-way profile,
        # this just for illustration
        self.proc_ss = np.zeros((count, 3))

        # Sort the data samples by depth
        sorted_ss = sorted(zip(self.obs_depth, self.obs_ss))

        layer = 0
        for d, ss in sorted_ss:
            self.proc_ss[[layer], [0]] = d
            self.proc_ss[[layer], [1]] = ss
            layer += 1

        # Identify all the depths for which there are multiple observations
        mask = np.full((count, 1), True)
        mask[1:, [0]] = np.diff(self.proc_ss[:, [0]], axis=0) != 0

        # Remove the duplicates - You really should get statistical representations here
        # but to keep this short just remove the duplicates
        self.proc_ss = self.proc_ss[mask[:, 0], ...]

        # Determine the gradients - Note the indexing: the gradient of the first layer 
        # is contained at the same index as the data for the TOP of the layer.
        self.proc_ss[0:-1, [2]] = np.diff(self.proc_ss[:, [1]],
                                          axis=0)/np.diff(self.proc_ss[:, [0]], axis=0)

        # Estimate gradient for last layer assuming that the temperature and salinity remain the same
        # gradient solely a function of pressure (depth)
        self.proc_ss[-1, [2]] = 0.017

        # Extend to 12000 m if necesarry - this is to get around some manufcturers requirements
        if self.obs_depth[-1] < 12000:
            ss = self.proc_ss[-1:, [1]] + self.proc_ss[-1:, [2]] \
             * (12000-self.proc_ss[-1:, [0]])
            self.proc_ss = np.vstack((self.proc_ss, [12000, ss, 0.017]))

        # Extend to 0 m if necesarry - assume well mixed
        if self.obs_depth[0] > 0:
            self.proc_ss = np.vstack(
                ([0, self.proc_ss[0, [1]], 0.], self.proc_ss))
            
        # Step 5 Create a look-up array of twtts for each full layer
        # Allows for great gain in efficiency (do not have to calculate for each ping)
        self.twtt_layer = np.zeros((count, 1))
        
        for layer in range(0,self.metadata["count"]-1):
            if self.proc_ss[layer, [2]] == 0:
                self.twtt_layer[layer] = 2 * \
                    (self.proc_ss[layer+1, [0]] - self.proc_ss[layer, [0]])/ \
                     self.proc_ss[layer, [1]]
            else:
                self.twtt_layer[layer] = 2 / self.proc_ss[layer, [2]] * \
                 log(self.proc_ss[layer+1, [1]]/self.proc_ss[layer, [1]])

    def depth(self, start_depth, twtt):
        # Determine depth relative to the transducer for vertical incidence data 
        # using the harmonic mean sound speed

        # Find the index to the start layer using Boolean logic
        # Note that the profile is extended to full ocean depth
        # As layer is indexed by the top depth this is a guarantee that there is a next layer
        start_i = sum(start_depth >= self.proc_ss[:, [0]])-1
        start_i = int(start_i[0])  # Turn the index into an integer

        # The sound speed at the top of the layer
        c_start=self.proc_ss[start_i, [1]]
        
        # Calculate the time it takes to get through this layer - avoid division by zero
        if self.proc_ss[start_i, [2]] == 0:
            twtt_layer = 2*(self.proc_ss[start_i+1, [0]] - \
                start_depth)/self.proc_ss[start_i, [1]]
        else:
            # Update c_start to comepensate for the change in speed to the transducer depth
            c_start += self.proc_ss[start_i, [2]]*(start_depth-self.proc_ss[start_i, [0]])
            
            twtt_layer = 2/self.proc_ss[start_i, [2]] * \
             log(self.proc_ss[start_i+1, [1]]/c_start)
            
        # Create a cumulative sum of the contributions of each layer
        twtt_cum=np.zeros((self.metadata["count"]))
        twtt_cum[start_i+1]=twtt_layer
        twtt_cum[start_i+2:]=twtt_layer+np.cumsum(self.twtt_layer[start_i+2:])
        
        # Determine where twtt_cum starts to exceed twtt using Boolean logic
        end_i=np.max(np.where(twtt_cum<twtt),axis=1)

        # Deal with the twtt in the final layer. There are two cases 
        # 1) the signal originated in a higher layer 
        # 2) the signal originated in the same layer  
        
        if start_i!=end_i:
            twtt_last=(twtt-twtt_cum[end_i])
            d_start=self.proc_ss[end_i, [0]]
        else: 
            twtt_last=twtt
            d_start=start_depth
            
        if  self.proc_ss[[end_i], [2]]!=0:
            # gradient is non-zero
            # From twtt=2/g*log(C2/C1) => C2=C1*exp(g*twtt/2)
            c_end = c_start*exp( self.proc_ss[end_i, [2]]*twtt_last/2)
            
            # The associated depth
            depth = d_start + (c_end-c_start)/self.proc_ss[end_i, [2]]
            
        else:
            # gradient is zero
            depth=d_start + self.proc_ss[end_i, [1]]*twtt_last/2

        # the depth is relative to the instantaneous (heaving) water surface, 
        # bring the value back to the transducer.
        
        depth-=start_depth   
        return depth

    def draw(self, full_profile=False, ax1=False, depth_range=False, ss_range=False, label=True):

        if ax1 == False:
            fig, ax1 = plt.subplots()   
            
        if full_profile:
            if depth_range == False:
                depth_range = (min(self.d), max(self.d))
            if ss_range == False:
                ss_range = (min(self.c), max(self.c))
            plt.plot(self.c[0:], self.d[0:])
        else:
            if depth_range == False:
                depth_range = (min(self.d[0:-1]), max(self.d[0:-1]))
            if ss_range == False:
                ss_range = (min(self.c[0:-1]), max(self.c[0:-1]))
            plt.plot(self.c[1:-1], self.d[1:-1])
            
        plt.ylim(depth_range)
        plt.xlim(ss_range)
        
        if label:
            plt.ylabel('← Depth [m]')
        else:
            labels = [item.get_text() for item in ax1.get_yticklabels()]
            empty_string_labels = ['']*len(labels)
            ax1.set_yticklabels(empty_string_labels)
            
        plt.xlabel('Sound Speed [m/s] →')
        ax1.invert_yaxis()
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('bottom')
        
        # Set the title from the file name that contained the data
        ax1.title.set_text(os.path.splitext(self.metadata['name'])[0])

