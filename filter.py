import scipy
from scipy import signal
from matplotlib import pyplot as plt
import numpy

def boxfilter(input, cutoff):
    fft=scipy.fft(signal) # (G) and (H)  
    bp=fft[:]  
    for i in range(len(bp)): # (H-red)  
        if i>=cutoff:bp[i]=0  
    ibp=scipy.ifft(bp) # (I), (J), (K) and (L)
    return ibp
         
def boxfilter2(input, cutoff):
    fil = signal.boxcar(cutoff)
    output = signal.convolve(input, fil)
    return output

input = numpy.array([0, 0.1, 0.2, 0.5, 0.45, 0.3, 0.8, 0.2, 0.3, 0.31])
ans = boxfilter2(input, 2)
plt.plot(range(len(input)), input, range(len(ans)), ans)
plt.show