# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 21:56:56 2020

@author: haris
"""
"""
This script is to filter the CCS fastq from PacBio into individual files for each
barcoded sequence
"""

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import SingleLetterAlphabet
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

# random_seq = Seq('GGTACACGCTGTGACAA', generic_dna)
# print(random_seq)
# print(random_seq.reverse_complement())
# print(random_seq.complement())

#List of PacBio Prefix Linkers:
Prefix_A='TCAGACGATGCGTCAT'
Prefix_B='CTATACATGACTCTGC'
Prefix_C='TACTAGAGTAGCACTC'
Prefix_D='TGTGTATCAGTACATG'
Prefix_E='ACACGCATGACACACT'
Prefix_F='GATCTCTACTATATGC'
Prefix_G='ACAGTCTATACTGCTG'
Prefix_H='ATGATGTGCTACATCT'

#List of PacBio Suffix Linkers:
Suffix_1='CATAGCGACTATCGTG'
Suffix_2='CATCACTACGCTAGAT'
Suffix_3='CGCATCTGTGCATGCA'
Suffix_4='TATGTGATCGTCTCTC'
Suffix_5='GTACACGCTGTGACTA'
Suffix_6='CGTGTCGCGCATATCT'
Suffix_7='ATATCAGTCATGCATA'
Suffix_8='GAGATCGACAGTCTCG'
Suffix_9='CACGCACACACGCGCG'
Suffix_10='CGAGCACGCGCGTGTG'
Suffix_11='GTAGTCTCGCACAGAT'
Suffix_12='GAGACTCTGTGCGCGT'


#The Target DNA sequenced via PacBio
TargetDNA='GTCCACAATTTTCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGCTAAAAGATTTGACAGCTAGCTCAGTCCTAGGGATTGTGCTAGCGCGTCCGGCGTAGAGGATCGAGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGGTACTAGAGGGCTCGTTGAACACCGTCTCAGGTAAGTATCAGTTGTAAAAAGAGGAGAAATAGTCCATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGATCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAAGGCTCGATCGGTGTGAAAAGTCAGTATCCAGTCGTGTAGTTCTTATTACCTGTCCCCTAGCATAACCCCGCGGGGCCTCTTCGGGGGACTCGCGGGGTTTTTTGCTGAAAGAATTATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGCTGCATTACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGGGCTCG'

#Create barcoded reference sequences for downstream alignment by combining Prefix+TargetDNA+Suffix:
Ref_1= str(Prefix_A)+str(TargetDNA)+str(Suffix_1)
print(Ref_1)
#print(len(Ref_1))

Ref_2= str(Prefix_B)+str(TargetDNA)+str(Suffix_2)
#print(len(Ref_2))

Ref_3= str(Prefix_C)+str(TargetDNA)+str(Suffix_3)
#print(len(Ref_3))

Ref_4= str(Prefix_D)+str(TargetDNA)+str(Suffix_4)
#print((Ref_4))

Ref_5= str(Prefix_E)+str(TargetDNA)+str(Suffix_5)
print((Ref_5))

Ref_6= str(Prefix_F)+str(TargetDNA)+str(Suffix_6)
#print(len(Ref_6))

Ref_7= str(Prefix_G)+str(TargetDNA)+str(Suffix_7)
#print(len(Ref_7))

Ref_8= str(Prefix_H)+str(TargetDNA)+str(Suffix_8)
print((Ref_8))

#SequenceRecords=[]
#M187_1=[]

# for index, record in enumerate(SeqIO.parse("PacBio_testfile.fastq", "fastq")):
#     print(
#         "index %i, ID = %s, length %i, with %i features"
#         % (index, record.id, len(record.seq), len(record.features))
#     )
#     print(record.format('fastq'))
#     #dir(record)


M187_1 = []
M187_2 = []
M180 = []
M153 = []
ATFU_1= []
ATFU_2 = []
T7pol = []
T86 = []
Wrong_barcode = []
Incorrect_barcode = 0

for record in SeqIO.parse('m54053_191003_153045.fastq', "fastq"):
    
    sequence = str(record.seq)
    
    #if sequence.find(Prefix_A) != -1 and sequence.find(Suffix_1) != -1:
    if fuzz.partial_ratio(Prefix_A, sequence) > 70 and fuzz.partial_ratio(Suffix_1, sequence) > 70:
        M187_1.append(record)
        #print(sequence)
        continue
    
    #elif sequence.find(Prefix_B) != -1 and sequence.find(Suffix_2) != -1:
    elif fuzz.partial_ratio(Prefix_B, sequence) > 70 and fuzz.partial_ratio(Suffix_2, sequence) > 70:
        M187_2.append(record)
        
        continue
    
    #elif sequence.find(Prefix_C) != -1 and sequence.find(Suffix_3) != -1:
    elif fuzz.partial_ratio(Prefix_C, sequence) > 70 and fuzz.partial_ratio(Suffix_3, sequence) > 70:
        M180.append(record)
        
        continue
        
    #elif sequence.find(Prefix_D) != -1 and sequence.find(Suffix_4) != -1:
    elif fuzz.partial_ratio(Prefix_D, sequence) > 70 and fuzz.partial_ratio(Suffix_4, sequence) > 70:
        M153.append(record)
        
        continue
        
    #elif sequence.find(Prefix_E) != -1 and sequence.find(Suffix_5) != -1:
    elif fuzz.partial_ratio(Prefix_E, sequence) > 70 and fuzz.partial_ratio(Suffix_5, sequence) > 70:
        ATFU_1.append(record)
        
        continue
        
    #elif sequence.find(Prefix_F) != -1 and sequence.find(Suffix_6) != -1:
    elif fuzz.partial_ratio(Prefix_F, sequence) > 70 and fuzz.partial_ratio(Suffix_6, sequence) > 70:
        ATFU_2.append(record)
        
        continue
        
    #elif sequence.find(Prefix_G) != -1 and sequence.find(Suffix_7) != -1:
    elif fuzz.partial_ratio(Prefix_G, sequence) > 70 and fuzz.partial_ratio(Suffix_7, sequence) > 70:
        T7pol.append(record)
        
        continue
        
    #elif sequence.find(Prefix_H) != -1 and sequence.find(Suffix_8) != -1:
    elif fuzz.partial_ratio(Prefix_H, sequence) > 70 and fuzz.partial_ratio(Suffix_8, sequence) > 70:
        T86.append(record)
        
    else:
        Wrong_barcode.append(record)
        Incorrect_barcode += 1

    
print('Number of MM187_1 sequences: ' + str(len(M187_1)))
print('Number of MM187_2 sequences: ' + str(len(M187_2)))
print('Number of MM180 sequences: ' + str(len(M180)))
print('Number of MM153 sequences: ' + str(len(M153)))
print('Number of ATFU_1 sequences: ' + str(len(ATFU_1)))
print('Number of ATFU_2 sequences: ' + str(len(ATFU_2)))
print('Number of T7pol sequences: ' + str(len(T7pol)))
print('Number of TS86 sequences: ' + str(len(T86)))
print(len(Wrong_barcode))
print('Wrong Barcoded Sequence Count = '+ str(Incorrect_barcode))

with open('M187_1.fastq', 'w') as output_file1:
    SeqIO.write(M187_1, output_file1, 'fastq')
    
output_file1.close()

with open('M187_2.fastq', 'w') as output_file2:
    SeqIO.write(M187_2, output_file2, 'fastq')
    
output_file2.close()

with open('M180.fastq', 'w') as output_file3:
    SeqIO.write(M180, output_file3, 'fastq')
    
output_file3.close()

with open('M153.fastq', 'w') as output_file4:
    SeqIO.write(M153, output_file4, 'fastq')
    
output_file4.close()

with open('ATFU_1.fastq', 'w') as output_file5:
    SeqIO.write(ATFU_1, output_file5, 'fastq')
    
output_file5.close()

with open('ATFU_2.fastq', 'w') as output_file6:
    SeqIO.write(ATFU_2, output_file6, 'fastq')
    
output_file6.close()

with open('T7pol.fastq', 'w') as output_file7:
    SeqIO.write(T7pol, output_file7, 'fastq')
    
output_file7.close()

with open('T86.fastq', 'w') as output_file8:
    SeqIO.write(T86, output_file8, 'fastq')
    
output_file8.close()

with open('Wrong_Barconde.fastq', 'w') as output_file9:
    SeqIO.write(Wrong_barcode, output_file9, 'fastq')
    
output_file9.close()


# ID = []
# Sequence = []
# Name = []
    


# for records in SeqIO.parse('PacBio_testfile.fastq', "fastq"):
   
#     #print((records.id))
#     print(records.letter_annotations['phred_quality'])
#     ID.append(str(records.id))
#     Sequence.append(str(records.seq))
#     Name.append(str(records.name))
    
   
# #print(ID, Sequence, Name)

# NewRecord = SeqRecord(str(Sequence), id=str(ID), name= (str(Name)))


# for data in ID and Sequence and Name:
#     NewRecord = (SeqRecord(
        
#         Seq = Sequence,
#         id = ID,
#         name= Name
#         ))
    
# print(NewRecord)
    
    
# #learned that you can only append one feature of SeqRecord into one list.. how to append all?    
#     if len(records.seq) > 1202:
#         print(records)
#     else:
#         print('sequence is shorter')
    
#     #print(records)


# for record in SeqIO.parse('PacBio_testfile.fastq', "fastq"):
#     if record.seq.find(Prefix_A) != -1 and record.seq.find(Suffix_1) != -1:
#         print(record)
#         #with open('M187_1.fastq') as handle:
#             #SeqIO.write(record, handle, 'fastq')
#         SeqIO.write(record, 'M187_1.fastq', 'fastq')
#     else:
#         print('Barcodes not found')
        

    


# def datafilter(records):
#     SequenceRecords=[]
#     M187_1=[]
#     M187_2=[]
#     M180=[]
#     M153=[]
#     ATFU_1=[]
#     ATFU_2=[]
#     T7pol=[]
#     T86=[]
    
#     for records in SeqIO.parse('PacBio_testfile.fastq', "fastq"):
#         SequenceRecords.append(records)
#     print(SequenceRecords)
        
#     return(records)



    
