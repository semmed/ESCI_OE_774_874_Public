import numpy as np
from numpy import log, log10


## NoiseLevel - Function estimating Noise Level based the Wenz curves
# Based on USNL Lecture notes found at:
# https://www.usna.edu/Users/physics/ejtuchol/documents/SP411/Chapter11.pdf
# On November 6, 2016 (Underwater Acoustics and Sonar SP411 Handouts and
# Notes Fall 2006)

# You'll find that many comments are copied verbatim from this source

# Semme J. Dijkstra  November 6,2016

def noise_level( SS,SD,fc,bw,verbose):
    # SS: Sea State
    # SD: Shipping density as ranging [1,4]

    ## 1 Set the frequency independent noise levels NL100 and NL1K

    if SD <= 1 or SD >= 7:
        error('NL: function not implemented to handle shipping densities outside the range SD=[1,7]')


    # Set the shipping density noise level
    NL100=60+(SD-1)*5;

    # Set the sea state noise level

    if SS==0:
        NL1K=null
    elif SS<=0.5:
        NL1K=null
    elif SS<=1:
        NL1K=null
    elif SS<=2:
        NL1K=null
    elif SS<=3:
        NL1K=null
    elif SS<=4:
        NL1K=null
    elif SS<=5:
        NL1K=null
    else:
        NL1K=null


    ## 2 Determine the upper and lower cutoff frequencies

    # To be able to determine the noise at either end of the signal's spectrum,
    # and thereby the mean noise level, we need to determine the cutoff 
    # frequencies. The cutoff frequency calculation depends on the bw relative 
    # to the central frequency fc. For a narrow band signal the fc is the
    # geometric mean of the cutoff frequencies, whereas for a wide band signal
    # the fc is the arithmetic mean.

    # if f2/f1 > 1.1 then the center frequency is the geometric mean
    # fc=sqrt(f1*f2), i.e. fc=sqrt(f1*(f1+bw) thus 
    # f1=(-bw+sqrt(bw^2+4*fc^2))/2 

    # If f2/f1 < 1.1 then the center frequency is the arithmetic mean 
    # fc=(f1+f2)/2, i.e. fc=(f1+f1+100)/2 thus
    # f1=(2*fc-bw)/2=fc-bw/2 

    # Obviously f2=f1+bw in both cases

    # Determine the cutoff frequencies
    
    f=np.zeros(2)
    if bw > fc/10.5:
        f[0]=null
    else:
        f[0]=null

    f[1]=null

    ## Estimate the noise for the lower and upper frequency

    NLship=np.zeros(2)
    NLss=np.zeros(2)
    NLss[:]=NL1K

    for i in range(2):
        if f[i] < 10:
            raise RuntimeError('NL: function not implemented for frequencies < 10 Hz')
        
        # 10-100 Hz ? Noise levels depend heavily on shipping density and industrial 
        # activities. Levels are typically in range of 60-90 dB with very little 
        # frequency dependence.

        if 10<f[i] and f[i]<100:
            NLship[i]=NL100
        
        # 100-1000 Hz ? Noise in this band is dominated by shipping (decreasing 
        # intensity with frequency increases). A significant contribution is also 
        # from sea surface agitation. Urick (1986) developed a model for predicting 
        # this shipping noise:

        if 100<f[i] and f[i]<1000:
            NLship[i]=NL100-null
        
    
        # 1-100 kHz ? Sea surface agitiation is now the dominant factor, unless 
        # marine mammals or rain is present. Knudsen (1948) presented a model 
        # to predict this contribution - A new model has been developed by 
        # APL (1994), it is more accurate but is more complex.
    
        if 1000<f[i] and f[i]<100000:
            NLss[i]=NL1K-null
        
    
        # Voluntary Challenge: Implement the 1994 APL model
    

    ## Noise Level by SS state as estimated by Wenz Curves 
    # (or rather the numerical models on which they are based)

    # The total NL=(NLship,ave+10*log10(bw)) 'level sum' (NLss,ave+10*log10(bw))

    # Combine the incoherent noise levels Ln using the 'level sum'
    # SPL=10*log10(10^(L1/10)+10^(L2/10)+...+10^(Ln/10)

    NL=...

    if verbose:
        # Let the user know what was found
        print('Lower cutoff frequency               : '+ str(f[0])+' Hz')
        print('Upper cutoff frequency               : '+ str(f[1])+' Hz')
        print('Average Shipping noise (Urick,1986)  : '+ str(sum(NLship)/2)+' dB re. 1uPa')
        print('Average Surface noise (Knudsen, 1984): '+ str(sum(NLss)/2)+' dB re. 1uPa')
        print('Total Ambient noise                  : '+ str(NL)+' dB re. 1uPa')

    return NL
