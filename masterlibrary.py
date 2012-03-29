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
	subprocess.Popen('makeblastdb -in ' + filename + '.fa -input_type fasta -dbtype prot -parse_seqids -out ' + filename + 'db -title ' + dbname, shell=True)


def psiblastme(queryname, dbname, resultname):
	"""Pass the query and database to a shell command which will execute PsiBLAST matching and output the result in resultname.csv"""
	subprocess.Popen('psiblast -db ' + dbname + ' -query ' + queryname + '.fa -out ' + resultname + '.csv -outfmt 10', shell=True)
	
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
		
		#Normalising the table to approximate the BLOSUM62 matrix original format
		
		self.lookuptable = 1 + numpy.round(-self.lookuptable/float(self.lookuptable[0, 2]))
					
	def save(self, filename='BLOSUM62'):
		"""A method to rewrite the psiBLAST substitution matrix stored at /usr/share/ncbi/data/(filename)"""
		scorefile = open('/usr/share/ncbi/data/' + filename, "w")
		protstring = 'ARNDCQEGHILKMFPSTWYVBJZX*'
		savematrix = numpy.zeros((25,25))
		savematrix[0:len(self.lookuptable), 0:len(self.lookuptable)] = self.lookuptable[:,:]
		j = 0
		scorefile.write("   " + " ".join(protstring) + "\n")
		for char in protstring:
			scorefile.write(char + '  ')
			for i in range(25):
				scorefile.write(str(int(savematrix[j, i])) + '  ')
			scorefile.write('\n')
			j = j + 1
		scorefile.close()
		
	
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
		self.filteredseries = numpy.copy(self.recarray)
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
		
	def simplestringconvert(self, nsampleseg, TPAobject):
		"""Use a naive form of PAA to convert the dataset to a character string for BLAST"""
		protstring = 'ARNDCQEGHILKMFPSTWYVBJZX*'
		try:
			fx = numpy.copy(self.filteredseries)
		except:
			print "No filtered series, using raw data"
			fx = numpy.copy(self.recarray)
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
			if avgval < TPAobject.breakpoints[0]:
				charac = protstring[0]
			
			elif avgval >= TPAobject.breakpoints[-1]:
				charac = protstring[len(TPAobject.breakpoints) + 1]
			
			else:
				for i in range(1, len(TPAobject.breakpoints)):
					if avgval < TPAobject.breakpoints[i] and avgval >= TPAobject.breakpoints[i-1]:
						charac = protstring[i]
		
			stringreturn = stringreturn + charac
			numbertodo = numbertodo + nsampleseg
			        	
		return stringreturn

	def mutateshift(self, start=0, end=-1, offset=0.5):
		self.recarray[start:end] = self.recarray[start:end] + offset
		
	def mutatenoise(self, stddev):
		self.recarray[:,:] = self.recarray[:,:] + numpy.random.normal(0, stddev, numpy.shape(self.recarray))
		
	def mutateinsert(self, dataset, datastart=0, dataend=-1, start=0):
		end = start + len(dataset.recarray[datastart:dataend, 0])
		self.recarray[start:end,0] = dataset.recarray[datastart:dataend, 0]
		
		
class Sequence():
	"""A class for string objects as sequences, a wrapper for BioPython SeqRecord objects"""
	def __init__(self, filename='default', format='fa', seq=None, identity='unknown',descrip='unknownsequence'):
		self.seqs = []
        
		if seq is not None:
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
		
	def simpledataconvert(self, resolution, TPAobject):
		"""Convert the sequence data back to approximated dataset by naive PAA, (process info is lost, only for demonstration purposes"""
		
		protstring = 'ARNDCQEGHILKMFPSTWYVBJZX*'
		data = [];
		for i in range(len(self.seqs)):
			for j in range(len(protstring)):
				if cmp(self.seqs[i], protstring[0]) == 0:
						value = TPAobject.breakpoints[0] - 0.5*(TPAobject.breakpoints[1] - TPAobject.breakpoints[0])
					
				elif cmp(protstring[-1], self.seqs[i]):
						value = TPAobject.breakpoints[-1] + 0.5*(TPAobject.breakpoints[-1] - TPAobject.breakpoints[-2])
					 
				elif cmp(protstring[i], self.seqs[i]) == 0:
						value = 0.5*(TPAobjec.breakpoints[j] - TPAobject.breakpoints[j-1])
					
			for i in range(resolution*2):
				data.append(value)
				
		return data


        
class Match():
	def __init__(self, filename='default', darray=None):
		if darray is not None:
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
Q = TPAdistance(6, x)
Q.save('BLOSUM62')
print Q.lookuptable
q = Sequence(seq=query.simplestringconvert(6, Q), descrip='QUERYTEST', identity='Q1')
q.save('querytest')
p = x.simplestringconvert(6, Q)
y = Sequence(seq=p, descrip='CONVERSIONTEST', identity='T1')
print y.simpledataconvert(6, Q)
y.save('dbtest')
dustmask('dbtest')
makedb('dbtest', 'dbtest')
psiblastme('querytest', 'dbtestdb', 'fulltest')
blastmatch = Match(filename='fulltest')
print blastmatch.evalue
print blastmatch.Qstart
#print len(y.simpledataconvert(2))
#print len(x.ts)
#print len(p)
#print y.seqs[0]

#print Q.breakpoints
#print Q.lookuptable
z = Dataset('Data')
p = Dataset('QueryDat')


