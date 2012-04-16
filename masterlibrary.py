import numpy
import scipy
from scipy import signal
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIStandalone
import subprocess
import time
import csv

def simpsons(x, fx, h):
	"""An algorithm to compute the weighted average value for fx over h*2 data points for the purposes of naive PAA"""
	n = (max(x) - min(x))/(h)
	nint = int(n)
	ranger = range(nint)
	incrementsize = h/(x[2] - x[0])
	Integ = 0
	for i in ranger:
		Integ = Integ + h*(fx[2*i*int(incrementsize)] + 4*fx[(2*i + 1)*int(incrementsize)] + fx[(2*i +2)*int(incrementsize)])/6
	return Integ

def dustmask(filename):
	"""Mask a FASTA file in preparation for being turned into a searchable BLAST database using dustmasker"""
	subprocess.Popen('dustmasker -in '+filename+'.fa -out '+filename+'mask.fa'+' -infmt fasta -outfmt fasta', shell=True)

def makedb(filename, dbname):
	"""Make a BLAST database out of a masked fasta file, named dbname"""
	subprocess.Popen('makeblastdb -in ' + filename + 'mask.fa -out ' + filename + 'db.fa -title ' + dbname, shell=True)

def psiblastme(queryname, dbname, resultname):
	"""Pass the query and database to a shell command which will execute PsiBLAST matching and output the result in resultname.csv"""
	subprocess.Popen('psiblast -db ' + dbname + '.fa -query ' + queryname + '.fa -out ' + resultname + '.csv -outfmt 10', shell=True)

class TPAdistance():
	"""A class to generate match scoring matrices for different alphabet sizes, from Keogh, Chiu and Lonardi"""
	def __init__(self, alphabetsize, dataset):
		self.n = alphabetsize
		p = numpy.sort(dataset.recarray.flatten())
		axis = numpy.linspace(0, 1, num=len(dataset.recarray), endpoint=True)
		self.breakpoints = numpy.interp(numpy.linspace(0, 1, self.n + 1, endpoint=True), axis, p)[1:-1]
	
		self.lookuptable = numpy.zeros((self.n, self.n))
	
		for q in range(self.n):
			for z in range(self.n):
				if numpy.abs(q - z) > 1:
					self.lookuptable[q, z] = self.breakpoints[numpy.max([q,z])-1] - self.breakpoints[numpy.min([q,z])]
				else:
					self.lookuptable[q, z] = 0
	
class Dataset():
	"""A class for dataset objects with filtering"""
	def __init__(self, filename):
		self.recs = numpy.recfromcsv(filename + '.csv', delimiter = ',')
		self.ts = []
		self.reclist = []
	
		for rec in self.recs:
#			if isinstance(rec[0], int)| isinstance(rec[0], float):
#				self.ts.append(rec[0])
			try:
				self.ts.append(time.mktime(time.strptime(rec[0], '"%Y/%m/%d %H:%M"')))
			except:
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
		"""Filter the data using a boxcar filter and store the values in filteredseries"""
		for i in range(len(self.recarray[0])):
			fil = signal.boxcar(cutoff)
			output = signal.convolve(self.recarray[:,i]/cutoff, fil, mode='same')
			self.filteredseries[:,i] = output
		return self.filteredseries
	
	def hamming(self, cutoff):
		"""Filter the data using a hamming filter and store the values in filteredseries"""
		for i in range(len(self.recarray[0])):
			fil = signal.hamming(cutoff)
			output = signal.convolve(self.recarray[:,i], fil, mode='same')
			self.filteredseries[:,i] = output
		return self.filteredseries
	
	def gaussian(self, cutoff):
		"""Filter the data using a hamming filter and store the values in filteredseries"""
		for i in range(len(self.recarray[0])):
			fil = signal.gaussian(cutoff, cutoff/6)
			output = signal.convolve(self.recarray[:,i], fil, mode='same')
			self.filteredseries[:,i] = output
		return self.filteredseries
		
	def plotraw(self):
		"""Plot the unfiltered dataset"""
#		plt.figure(figurenr)
		plt.plot(self.ts, self.recarray)
		plt.show()
		
	def plotfiltered(self):
		"""Plot the filtered dataset (After a filtering method has been applied to the dataset)"""
#		plt.figure(figurenr)
		plt.plot(self.ts, self.filteredseries)
		plt.legend('123456789')
		plt.show()
	
	def plot(self):
		"""Plot the original and filtered datasets together"""
#		plt.figure(figurenr)
		plt.plot(self.ts, self.recarray, self.ts, self.filteredseries)
		plt.legend('123456789')
		plt.show()
		
	def simplestringconvert(self, nsampleseg):
		"""Use a naive form of PAA to convert the dataset to a character string for BLAST"""
		fx = numpy.copy(self.filteredseries)
		nsampleseg = nsampleseg*2
		rem = (len(fx) - 1) % nsampleseg
		fx = fx[rem:len(fx)]
		print "Removed ", rem, " data points from beginning of series"
		stringreturn = ""
		tsegment = nsampleseg
		nsegments = len(fx)/nsampleseg
		h = tsegment/2
		numbertodo = 1
		while numbertodo < len(fx):
			avgval = simpsons(range(0, (nsampleseg+1)), fx[(numbertodo-1):(numbertodo+nsampleseg)], h)/nsampleseg
			if avgval <= 0.25:
				charac = "E"
			elif avgval <= 0.43:
				charac = "D"
			elif avgval <= 0.57:
				charac = "T"
			elif avgval <= 0.75:
				charac = "M"
			else:
				charac = "C"
			stringreturn = stringreturn + charac
			numbertodo = numbertodo + nsampleseg
			        	
		return stringreturn

	def mutateshift(self, start=0, end=-1, offset=0.5):
		self.recarray[start:end] = self.recarray[start:end] + offset
		
	def mutatenoise(self, stddev):
		self.recarray = self.recarray + numpy.random.normal(0, stddev, len(self.recarray))

class Sequence():
	"""A class for string objects as sequences, a wrapper for BioPython SeqRecord objects"""
	def __init__(self, filename='default', format='fa', seq=None, identity='unknown',descrip='unknownsequence'):
		self.seqs = []
        
		if seq:
			self.list = SeqRecord(MutableSeq(seq), id=identity, description=descrip)
		else:
			handle = open(filename + '.' + format)
			self.list = list(SeqIO.parse(handle, "fasta"))
			handle.close()
		for seqi in self.list:
			self.seqs.append(MutableSeq(seqi))
            
	def __eval__(self):
		return self.dict
    
	def __call__(self, num=0):
		return self.dict[num]
    
	def save(self, filename):
		"""Write the sequence objects to a fasta file for access by PsiBLAST"""
		handle = open(filename + '.fa', "w")
		SeqIO.write(self.list, handle, "fasta")
		handle.close()
		
	def simpledataconvert(self, resolution):
		"""Convert the sequence data back to approximated dataset by naive PAA, (process info is lost, only for demonstration purposes"""
		data = []
		for i in range(len(self.seqs)):
			if self.seqs[i] == "E":
				val = 0.125
			elif self.seqs[i] == "D":
				val = 0.34
			elif self.seqs[i] == "T":
				val = 0.5
			elif self.seqs[i] == "M":
				val = 0.64
			else:
				val = 0.875
			for i in range(resolution*2):
				data.append(val)
				
		return data
		
        
class Match():
	def __init__(self, filename='default', darray = None):
		if darray:
			self.matchrecs = darray
		else:
			handle = open(filename + '.csv')
			recs = csv.reader(handle, delimiter = ',')
			self.matchrecs = recs.next()
			handle.close()
		
		self.QuerySeqID = numpy.copy(self.matchrecs[0])
		self.SubjectSeqID = numpy.copy(self.matchrecs[1])
		self.Pident = float(self.matchrecs[2])
		self.length = int(self.matchrecs[3])
		self.mismatches = int(self.matchrecs[4])
		self.gaps = int(self.matchrecs[5])
		self.Qstart = int(self.matchrecs[6])
		self.Qend = int(self.matchrecs[7])
		self.Sstart = int(self.matchrecs[8])
		self.Send = int(self.matchrecs[9])
		self.evalue = float(self.matchrecs[10])
		self.bitscore = float(self.matchrecs[11])
		
		del self.matchrecs

#TEST CODE
#NOTE: This will be refactored to work with the unittest library ASAP


x = Dataset('QueryDat')
query = Dataset('QueryDat')
query.boxfilter(2)
q = Sequence(seq=query.simplestringconvert(2), descrip='QUERYTEST', identity='Q1')
q.save('querytest')
x.boxfilter(2)
p = x.simplestringconvert(2)
y = Sequence(seq=p, descrip='CONVERSIONTEST', identity='T1')
y.save('conversiontest')
dustmask('conversiontest')
makedb('conversiontest', 'DATABASETEST')
psiblastme('conversiontest', 'conversiontestdb', 'fulltest')
#handle = open('fulltest.csv')
#recs = csv.reader(handle, delimiter = ',')
#print [float(x) if '.' in x else int(x) if isempty(RWE(x)) for x in recs.next()]
#print [try float(x) except x for x in recs.next()]
blastmatch = Match(filename='fulltest')
#print blastmatch.evalue
#print blastmatch.Qstart
#print len(y.simpledataconvert(2))
#print len(x.ts)
#print len(p)
#print y.seqs[0]
#Q = TPAdistance(12, x)
#print Q.breakpoints
#print Q.lookuptable
z = Dataset('Data')
z.mutatenoise(0.5)
z.plotraw()
plt.show()



