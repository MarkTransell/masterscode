import numpy
import scipy
from scipy import signal
from matplotlib import pyplot as plt

class dataset():
	"""A class for dataset objects with filtering"""
	def __init__(self, filename):
		self.recs = numpy.recfromcsv(filename + '.csv', delimiter = ',')
		
		self.ts = []
		self.reclist = []
	
		for rec in self.recs:
			self.ts.append(rec[0])
			self.rectemp = []
			for i in range(1, len(rec)):
				self.rectemp.append(rec[i])
			self.reclist.append(self.rectemp)
		self.recarray = numpy.array(self.reclist)
		del self.rectemp
    
	def __eval__(self, t):
		return self.recarray
	
	def __call__(self, t):
		return self.recarray
	
	def boxfilter(self, cutoff):
		fft=scipy.fft(self.recarray[:,0]) # (G) and (H)  
		bp=fft[:]  
		for i in range(len(bp)): # (H-red)  
			if i>=cutoff:bp[i]=0  
		ibp=scipy.ifft(bp) # (I), (J), (K) and (L)
		return ibp
         
	def boxfilter2(self, cutoff):
		fil = signal.boxcar(cutoff)
		output = signal.convolve(self.recarray[:,0], fil)
		return output
		
x = dataset('Data')
##print x.recarray
print x.boxfilter2(1)
##print x.ts
