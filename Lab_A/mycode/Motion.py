import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from numpy import pi, cos, sin, log, exp
import numpy as np

class Motion:
    """A Class for handling motion Data"""

    def __init__(self):

        # The data attributes
        self.times = list()
        self.yaw = list()
        self.roll = list()
        self.pitch = list()
        self.heave = list()
        self.metadata = dict()
        self.metadata["units"] = "rad"
        self.metadata["start_time"] = None
        self.metadata["end_time"] = None
        self.metadata["count"] = None
        self.metadata["time_basis"] = "UTC"

    # The I/O methods:

    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening motion data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        motion_file = open(fullpath)
        motion_content = motion_file.read()
        motion_file.close

        # Tokenize the contents
        motion_lines = motion_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for motion_line in motion_lines:
            observations = motion_line.split()  # Tokenize the string
            time = datetime.fromtimestamp(
                float(observations[5]), timezone.utc)
            self.times.append(time)
            self.yaw.append(float(observations[6])*pi/180)
            self.roll.append(float(observations[7])*pi/180)
            self.pitch.append(float(observations[8])*pi/180)
            self.heave.append(float(observations[9]))
            count += 1

    def draw(self):
        print('Drawing Motion Data')
        plt.figure(figsize=(20, 10))
        ax1 = plt.subplot(4, 1, 1)
        plt.plot(self.times, np.degrees(self.yaw))
        plt.ylabel('Heading [deg] →')
        ax2 = plt.subplot(4, 1, 2, sharex=ax1)
        plt.plot(self.times, self.heave)
        plt.ylabel('Heave [m] →')
        ax3 = plt.subplot(4, 1, 3, sharex=ax1)
        plt.plot(self.times, np.degrees(self.roll))
        plt.ylabel('Roll [deg] →')
        ax4 = plt.subplot(4, 1, 4, sharex=ax1, sharey=ax3)
        plt.plot(self.times, np.degrees(self.pitch))
        plt.ylabel('Pitch [deg] →')
        plt.xlabel('Time ('+self.metadata['time_basis']+') →')
        plt.xticks(rotation='60')

        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)