import os
from datetime import datetime, timezone
from numpy import pi, cos, sin, log, exp
import numpy as np
import matplotlib.pyplot as plt

class EchosounderData:
    """A Class for handling Two Way Travel Time Data"""

    def __init__(self):

        # The data attributes
        self.times = np.array([])
        self.twtts = np.array([])
        self.metadata = dict()
        self.metadata["units"] = "s"
        self.metadata["start_time"] = None
        self.metadata["end_time"] = None
        self.metadata["count"] = None
        self.metadata["time_basis"] = "UTC"
        self.metadata["Source_File"]=str()
        
    # The I/O methods:

    def read_jhc_file(self, fullpath):

        # Check the File's existence
        if os.path.exists(fullpath):
            self.metadata["Source File"] = fullpath
            print('Opening Two Way Travel Time (TWTT) data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        twtt_file = open(fullpath)
        twtt_content = twtt_file.read()
        twtt_file.close
        
        times=list()
        twtts=list()
        
        # Tokenize the contents
        twtt_lines = twtt_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for twtt_line in twtt_lines:
            observations = twtt_line.split()  # Tokenize the string
            #time=datetime.fromtimestamp(float(observations[5]), timezone.utc)
            times.append(datetime.fromtimestamp(float(observations[5]), timezone.utc))
            twtts.append(float(observations[6]))
            count += 1
            
        self.times=np.asarray(times)
        self.twtts=np.array(twtts)
     
    def draw(self):
        plt.figure(figsize=(10, 10))
        print('Drawing TWTT Data')
     # plotting the points  
        plt.plot(self.times, self.twtts) 
        plt.title('Two Way Travel Times in [s]') 
        plt.ylabel('TWTT in [s] →') 
        plt.xlabel('Time time ('+self.metadata['time_basis']+') →') 
        plt.xticks(rotation='60')