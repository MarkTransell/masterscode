import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
import time
import os
from matplotlib import pyplot as plt


ts = numpy.recfromcsv('LPG Data Set_1_n.csv')
tstamps = []
recarray = []
for rec in ts:
	tstamps.append(rec[0])
	rectemp = []
	for i in range(1, len(rec)):
		rectemp.append(rec[i])
	recarray.append(rectemp)

#ts = numpy.loadtxt('synthetic.data')
#t = ts[:, 0]
#fx = ts[:, 1]
#fx = (fx - numpy.nanmin(fx))/(numpy.nanmax(fx) - numpy.nanmin(fx))
#ans = ''
#recs = range(10)

def convert(fx, nsampleseg):
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
            charac = "G"
        elif avgval <= 0.5:
            charac = "C"
        elif avgval <= 0.75:
            charac = "A"
        else:
            charac = "T"
        stringreturn = stringreturn + charac
        numbertodo = numbertodo + nsampleseg
    return stringreturn
    
def simpsons(x, fx, h):
    n = (max(x) - min(x))/(h)
    nint = int(n)
    ranger = range(nint)
    incrementsize = h/(x[2] - x[0])
    Integ = 0
    for i in ranger:
        Integ = Integ + h*(fx[2*i*int(incrementsize)] + 4*fx[(2*i + 1)*int(incrementsize)] + fx[(2*i +2)*int(incrementsize)])/6
    return Integ
    
def makerecord(sequence, idin = "gi|00000001|gb|AEE00001.001", descriptionin = "unknown"):
  recor = SeqRecord(Seq(sequence, generic_dna), id=idin, description=descriptionin)
  return recor
  
def readrecord(filename):
    handle = open(filename + '.fa')
    dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return dict
  
def writemulti(recordlist, filename):
  handle = open(filename + '.fa', "w")
  SeqIO.write(recordlist, handle, 'fasta')
  handle.close()
  
def makedbnucleo(filename, title = "Sequence"):
	os.system("dustmasker -in " + filename +".fa -parse_seqids -outfmt fasta -out " + filename +".mfa")
	os.system("dustmasker -in " + filename +".fa -parse_seqids -outfmt maskinfo_asn1_bin -out " + filename +"_dust.asnb")
	os.system('convert2blastmask -in ' + filename +'.mfa -parse_seqids -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out ' + filename + '_mfa.asnb')
	os.system('windowmasker.exe -in ' + filename + '.fa -infmt fasta -mk_counts true -parse_seqids -out ' + filename +'_mask.count')
	os.system('windowmasker.exe -in ' + filename + '.fa -infmt fasta -ustat ' + filename + '_mask.count -outfmt maskinfo_asn1_bin -parse_seqids -out '+ filename +'_mask.asnb')
	print("Masking done")
	os.system('makeblastdb -in ' + filename + '.mfa -dbtype nucl -parse_seqids -mask_data ' + filename + '_mfa.asnb -out ' + filename + '_mfa -title '+ title)
	print("Database made")
	masktitle = title + '_masked'
	os.system('makeblastdb -in ' + filename + '_mfa -dbtype nucl -parse_seqids -mask_data ' + filename + '_mask.asnb -out ' + filename + '_mfa -title ' + masktitle)
	os.system('blastdbcmd.exe -db ' + filename + '_mfa -info')
	

def queryblastn(query, db, outfmt = '7'):
    os.system("dustmasker -in " + query +".fa -parse_seqids -outfmt fasta -out " + query +".mfa")
    os.system('blastn.exe -query ' + query + '.mfa -task megablast -db ' + db + '_mfa -outfmt ' + str(outfmt) + '-out ' + query + 'out')
    os.system('blastn.exe -query ' + query + '.fa -task megablast -db ' + db + '.fa -outfmt ' + str(outfmt) + '-out ' + query + 'out2')
    
def querysubstr(querystring, filename):
    processdict = readrecord(filename)
    its = processdict.items()
    longsets = []
    for i in range(len(its)):
        data = LCSubstr_lenloc(querystring, its[i][1].seq.tostring())
        longsets.append((data[0], (str(its[i][0]) + ': ' + str(data[1]))))
    return longsets    
    
    

def LCSubstr_lenloc(S, T): #Wikibooks, edited to include location
     m = len(S); n = len(T)
     L = [[0] * (n+1) for i in xrange(m+1)]
     lcs = 0
     loc = -1
     for i in xrange(m):
         for j in xrange(n):
             if S[i] == T[j]:
                 L[i+1][j+1] = L[i][j] + 1
                 lcs = max(lcs, L[i+1][j+1])
                 loc = j
     return (lcs, loc)

def LCSubstr_set(S, T): #Wikibooks
     m = len(S); n = len(T)
     L = [[0] * (n+1) for i in xrange(m+1)]
     LCS = set()
     longest = 0
     for i in xrange(m):
         for j in xrange(n):
             if S[i] == T[j]:
                 v = L[i][j] + 1
                 L[i+1][j+1] = v
                 if v > longest:
                     longest = v
                     LCS = set()
                 if v == longest:
                     LCS.add(S[i-v+1:i+1])
     return LCS

def findmaxmatch(results):
    maxi = 0
    loc = []
    for i in range(len(results)):
        if results[i][0] == maxi:
            loc.append(results[i][1])
        elif results[i][0] > maxi:
            maxi = results[i][0]
            loc = [(results[i][1])]
    return (loc, maxi)

for i in range(10):
    fx = ts[:, i]
    fx = (fx - numpy.nanmin(fx))/(numpy.nanmax(fx) - numpy.nanmin(fx))
    ans = convert(fx, 1)
    recs[i] = makerecord(ans, idin=('gi|000000' + str(i) + '|gb|AEE00001.' + str(i)), descriptionin=("Synthrepeat" + str(i)))

writemulti(recs, "Repeatdata")
processdict = readrecord("Repeatdata")
test = querysubstr("AAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTGGGGG", "Repeatdata")

maxi = findmaxmatch(test)
print "Found maximum substring of length " + str(maxi[1]) + ' in ' + str(len(maxi[0][:])) + ' locations: '
for i in range(len(maxi[0][:])):
    print maxi[0][i]

    

#makedbnucleo("Repeatdata", title = "RepeatingSeqData")
#os.system('blastn.exe -query Querytest.mfa -task megablast -db Repeatdata_mfa -outfmt 7 -out Querytest.out')
#os.system('blastn.exe -query Querytest.fa -task megablast -db Repeatdata.fa -outfmt 7 -out Querytest.out2')
#queryblastn("Querytest", "Repeatdata")

