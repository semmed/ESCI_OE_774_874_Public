from numpy import log10,cos,pi
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def modified_lambertian( th, Sb, verbose ):
    # modified Lambertian backscater angular response curve
    # th: depression angle
    # Sb: Bottom backscatter strength (Bottom gain)

    print(th)
    # Determine the number of bottom types
    n_bot=len(Sb[0])

    # Convert th to Incidence angles
    # I prefer everything as depression angles, but unfortunately many in the
    # literature mix coordinate systems - we will follow suit...
    th=np.tile(pi/2-th,(n_bot,1))    

    ## Lambertian component

    rough = Sb[1].reshape(n_bot,1)
    BSl = rough+10*log10(cos(th)**2)

    ## Specular component
    # Calculate specular contributions from all incidence angles
    spec = Sb[2].reshape(n_bot,1)
    crit = Sb[3].reshape(n_bot,1)
    BSs=spec*(1+cos(pi*(th/crit)))/2
    
    # Determine the range of incidence angles inside the critical angle
    # for all bottom types
    rCA=abs(th)<crit

    # Now Set all contributions from outside the critical angle to zero
    BSs[~rCA]=0
    
    ## Total Backscatter Angular Response
    BS=BSs+BSl

    if verbose:
        
        # Create a figure
        fig = plt.figure(figsize=(10, 12))
        
        # Add a figure Title
        fig.suptitle("Modified Lambertian Scattering")
        ax=list()
        for i in range(n_bot):

            if i == 0:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)-1))
            else:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)-1,sharex=ax[0],sharey=ax[0]))
                
            plt.plot(th[i,:]*180/pi,BSs[i,:],'b',linewidth=2)
        
            # Plot a grid
            plt.grid()
            
            # Add a title
            ax[-1].set_title('Specular Response: '+Sb[0][i])
            
            # Add the labels
            plt.xlabel('Angle of Incidence [deg]')
            plt.ylabel('Backscatter Strength BS [dB]')

            if i == 0:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1)))
            else:
                ax.append(fig.add_subplot(n_bot,2,2*(i+1),sharex=ax[1],sharey=ax[1]))
            
            # Plot the La,bertian Scatter
            plt.plot(th[i,:]*180/pi,BSl[i,:],'g',linewidth=2)
            
            # And the modified Lambertian Scatter
            plt.plot(th[i,:]*180/pi,BS[i,:],'k',linewidth=2)
            
            # Add a legend
            plt.legend(['Lambertian','Combined'])
        
            # Plot a grid
            plt.grid()
            
            # Add a title
            ax[-1].set_title('Modified Lambertian Response: '+Sb[0][i])
            
            # Add the labels
            plt.xlabel('Angle of Incidence [deg]')
            plt.ylabel('Backscatter Strength BS [dB]')

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=.5)
        plt.show()
        
    return BS
