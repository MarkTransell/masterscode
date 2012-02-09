import numpy
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord


class Sequence():
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
        handle = open(filename + '.fa', "w")
        SeqIO.write(self.list, handle, "fasta")
        handle.close

j = MutableSeq('blahlalala')
print j

x = Sequence(filename='test')
print x.list
p = MutableSeq(x.list[0])
x.seqs[0]
handle = open('test2.fa', "w")
SeqIO.write(x.list[0], handle, "fasta")
handle.close
x.save('test3')
a = Sequence(seq="ABASGWTHYETHSGSHS")
a.save('test4')
