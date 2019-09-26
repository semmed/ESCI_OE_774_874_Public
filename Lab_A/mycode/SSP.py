import os
from datetime import datetime, timezone


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
        self.proc_ss = np.array([])
        self.twtt_layer=np.array([])

        self.metadata = dict()
        self.metadata["units"] = "rad"
        self.metadata["count"] = None
        self.metadata["geoid"] = None
        self.metadata["ellipsoid"] = None
        self.metadata["chart_datum"] = None
        self.metadata["time_basis"] = "UTC"

    # The I/O methods:

    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening sound speed profile data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        motion_file = open(fullpath)
        motion_content = motion_file.read()
        motion_file.close

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

    

   