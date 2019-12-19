import os
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from numpy import pi, cos, sin, log, exp, arctan2, sqrt, nan
import numpy as np

class Baseline:
    """A Class for handling Baseline Data"""

    def __init__(self,arg1, *argv):
        line = 'Baseline() must be constructed as Baseline(Baseline), ' + \
               'Baseline(numpy.array, std. deviation),  ' + \
               'numpy(el[1], el[2], el[...], el[n], (std. dev,), or ' + \
               'numpy(el[1], el[2], el[...], el[n], (a,b), or\n\n' + \
               'An example of a 1D tuple is (0.1,) - NOTE the comma!'
        n_args = len(argv)
        if n_args == 0:
            if not isinstance(arg1,Baseline):        # Copy constructor
                raise RuntimeError(line)
            self.bl = arg1.bl
            self.a = arg1.a
            self.b = arg1.b
            self.sd = arg1.sd          
            return
        
        # Constructor using scalars, at least 2D baseline
        
        elif n_args == 1:
            if not isinstance(arg1,np.ndarray):          # Constructor using numpy array and std. dev.
                raise RuntimeError(line)
            self.bl = arg1
            self.sd = argv[0]
            self.a = nan
            self.b = nan
            return                
        
        # Constructor using elements (there must be at least two) and a tuple holding either a and b or sd
        if not isinstance(argv[-1],tuple):          # Constructor using numpy array and std. dev.
            raise RuntimeError(line)
        
        self.bl = np.zeros(n_args).reshape((n_args,1))
        self.bl[0] = [arg1]
        for i in range(n_args-1):
            self.bl[i+1] = argv[i]
            
        if len(argv[-1]) == 1:
            self.sd = argv[-1][0]
            self.a = nan
            self.b = nan
        elif len(argv[-1]) == 2:  
            self.a = argv[-1][0]
            self.b = argv[-1][1]
            self.sd = self.a+self.b*np.asscalar(sqrt(self.bl.T@self.bl))
        return
                                       
    
    def __str__ (self):
        return 'Baseline vector: '+str(self.bl.T)+'.T\n'+'Std dev: '+str(self.sd)+'\n'+'a: '+str(self.a)+', b: '+str(self.b)
        
    # adding two baselines  
    def __add__(self, o): 
        return Baseline(self.bl + o.bl)
    
    # subtracting two baselines  
    def __sub__(self, o): 
        return Baseline(self.bl - o.bl)
    
    # Azimuth of baseline
    def azimuth(self):
        return arctan2(self.bl[0],self.bl[1])
    
    # Length of the baseline
    def norm(self):
        return np.asscalar(sqrt(self.bl.T@self.bl))
    
    # Have the baseline draw itself
    def draw(self,*argv):    
        # For now a 2D plot
        n_args = len(argv)
        origin = [0], [0] # By default origin point is 0,0
        for arg in argv:
            if isinstance(arg,tuple):
                if len(arg) == 2:
                    origin = [arg[0]], [arg[1]]
 
        # Plot the baseline as an arrow
        plt.quiver(*origin, self.bl[0,0],self.bl[1,0], color=['r'], angles='xy', scale_units='xy', scale=1)

    # Find the max baseline element
    def max(self):
        return np.max(self.bl)
    
    # Determine the misclosure between this baseline and the baseline bl passed in. Report the 
    # misclosure as a baseline object w, and also as a np.array ppm holding the element differences divided by
    # the norm of this baseline as parts per million
    def misclosure(self,bl):
        w = self - bl
        ppm = w.bl/self.norm()*10**6
        return w,ppm
        


        
        
        