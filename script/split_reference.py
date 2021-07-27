#!/usr/bin/env python3
import glob
import os
import sys
import string


ref_fasta_file=sys.argv[1]
output_dir=sys.argv[2]
if (output_dir[-1] != "/"):
	output_dir=output_dir+"/"

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
    sequence_ok=sequence.translate({ord(c): None for c in string.whitespace}).upper()
    fasta_file.write(header.rstrip("\n")+"\n")
    fasta_file.write(sequence_ok+"\n")
    fasta_file.close() 

R_tableBLAST=open(output_dir+"R_table_analysis.csv","w")
R_tableBLAST.write("REF_NUM;REF_NAME\n")
cpt=1
for header,sequence in fasta_list.items():
    R_tableBLAST.write(str(cpt)+";"+header[1:].strip()+"\n")
    cpt=cpt+1
R_tableBLAST.close()
