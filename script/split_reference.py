#!/usr/bin/env python3
import glob
import os
import sys

ref_fasta_file=sys.argv[1]
output_dir=sys.argv[2]

def read_fasta(file):   
    fastas={}
    fasta_file=open(file,"r")
    for line in fasta_file:
        if (line[0]=='>'):
            header=line
            fastas[header]=''
        else:
            fastas[header]+=line
    fasta_file.close()
    return fastas

fasta_list=read_fasta(ref_fasta_file) 
 

for header,sequence in fasta_list.items():
    filename=output_dir+header[1:].strip()+".fasta"
    fasta_file=open(filename,'w')
    fasta_file.write(header.rstrip("\n")+"\n")
    fasta_file.write(sequence.rstrip("\n")+"\n")
    fasta_file.close() 



