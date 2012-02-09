import numpy
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import time
from matplotlib import pyplot as plt

class querysequence(dbsequence):
    def __init__(self, seq):
        self.dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        return self.dict
    
    def __eval__(self):
        return self.dict
    
    def __call__(self, num=0):
        return self.dict[num]
    

        