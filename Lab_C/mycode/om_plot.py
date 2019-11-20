from numpy import nan, arctan2, pi, sin
import numpy as np
import matplotlib.pyplot as plt

#Function to display array as a wedge display as commonlu used for multibeam and other swath data
# To make this more efficient - say in the case we want to display the data for every ping, we will leave memory allocation
# to the caller - this may obviate the need to reallocate memory all the time

# Semme J. Dijkstra     11/20/2019

class om_map:
    """A Class for mapping rectangular arrays to visualizations used in Ocean Mapping"""
    def __init__(self):
        self.im = np.array([])     # Array to hold image data
        self.X = np.array([])      # Array with X indexes for all pixels in im
        self.Y = np.array([])      # Array with Y indexed for all pixels in im
        self.row = np.array([])    # Array with mapping of Y indexes from data matrix to image matrix
        self.col = np.array([])    # Array with mapping of X indexes from data matrix to image matrix
        self.scale = 1.
        
        self.type = None
        
    def map_wedge(self,th, r, y_max):
        
        self.type = 'wedge'
    
        # Allocate memory for the image and populate it with nan's
    
        self.im = np.empty((y_max,2*y_max+1))
        self.im[:] = nan
    
        # Create indexing arrays
    
        self.X = np.tile(np.asarray(range(2*y_max+1)).reshape(1,2*y_max+1),(y_max,1))
        self.Y = np.tile(np.asarray(range(1*y_max)).reshape(y_max,1),(1,2*y_max+1))
   
        # Determine the ranges associated to the data
    
        th_min = np.min(th)
        th_max = np.max(th)
        th_int = th[1]-th[0]
        r_min = np.min(r)
        r_max = np.max(r)
        r_int = r[1]-r[0]
    
        # Determine the scale factor
    
        self.scale = self.Y[-1,0]/r_max

        # Allocate memory for the angles and ranges
    
        self.col = np.copy(self.X)
        self.col[:] = 0
        th_im = np.empty(self.X.shape)
        th_im[:] = nan
        self.row = np.copy(self.X)
        self.row[:] = 0
        r_im = np.empty(self.X.shape)
        r_im[:] = nan
        X_2 = self.X[0,-1]/2

    
        for i in range(self.X[0,-1]+1):
            for j in range(self.Y[-1,0]+1):

                # Calculate the depression angle and range of the point in data
                # space represented by the point X[j,i],Y[j,i] in image space
                th_im[j,i]=arctan2(self.Y[j,i],self.X[j,i]-X_2)
                r_im[j,i]=(((self.X[j,i]-X_2)**2+self.Y[j,i]**2)**.5)/self.scale
            
                # Look up the data value associated to the range r and
                # depression angle th in the data matrix, but only if they
                # are in the range r<r_max
            
                if (r_min<=r_im[j,i] and r_im[j,i]<=r_max) and (th_min<=th_im[j,i] and th_im[j,i]<=th_max):              
                    self.col[j,i]=int((th_im[j,i]-th_min)/th_int)
                    self.row[j,i]=int((r_im[j,i]-r_min)/r_int)
        

    def plot(self,data,title,func):
        
        data[0,0] = nan
        for i in range(len(self.row[0,:])):
            for j in range(len(self.col[:,0])):
                self.im[self.Y[j,i],self.X[j,i]]=data[self.row[j,i],self.col[j,i]]
   
        fig, ax = plt.subplots(figsize=(10, 12))

        cax = ax.imshow(self.im,)
        
        # Set the labels for the axes
        
        if self.type == 'wedge':
            locs, labels = plt.xticks()
            locs = locs[1:-1]
            labels = list()
            for loc in locs:
                labels.append('%.0f' % ((loc-self.Y[-1,0] - 1)/self.scale))
            plt.xticks(locs, labels)
            ax.set_xlabel('Horizontal Distance from Tx [m]' )
            
            locs, labels = plt.yticks()
            locs = locs[1:]
            labels = list()
            for loc in locs:
                labels.append('%.0f' % (loc/self.scale))
            plt.yticks(locs, labels)
            ax.set_ylabel('Vertical Distance from Tx [m]' )
        ax.set_title(title)
    
        mi = np.nanmin(data)
        ma = np.nanmax(data)
    
        tcks = np.arange(mi,ma,(ma-mi)/10)
        
        if func == '1':
            tck_label = str(tcks)
        elif func == '10**':
            print(tcks)
            tck_label = str(10**tcks)
        else:
            raise RuntimeError('scaling function currently not available')
        
        cbar = fig.colorbar(cax, ticks=tcks, orientation='vertical',shrink=.3)
        cbar.ax.set_xticklabels(tck_label)  # horizontal colorbar
        

        
        
    
    def colorbar(self,data,func):
        mi = np.min(data)
        ma = np.max(data)
        
        self.tcks = np.arange(mi,ma,(ma-mi)/5)
        
        if func == '1':
            cbar = fig.colorbar(cax, ticks=tcks, orientation='vertical',shrink=.3)
            cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
            
        elif func == '10**':
            pass