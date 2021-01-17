# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 20:01:11 2020

@author: haris
"""

from Bio import SeqIO
from Bio import pairwise2
import numpy as np
import os


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
print(len(TargetDNA))
#Create barcoded reference sequences for downstream alignment by combining Prefix+TargetDNA+Suffix:
Ref_1= str(Prefix_A)+str(TargetDNA)+str(Suffix_1)
#print(Ref_1)
#print(len(Ref_1))

Ref_2= str(Prefix_B)+str(TargetDNA)+str(Suffix_2)
#print(Ref_2)

Ref_3= str(Prefix_C)+str(TargetDNA)+str(Suffix_3)
#print(Ref_3)

Ref_4= str(Prefix_D)+str(TargetDNA)+str(Suffix_4)
#print((Ref_4))

Ref_5= str(Prefix_E)+str(TargetDNA)+str(Suffix_5)
#print((Ref_5))

Ref_6= str(Prefix_F)+str(TargetDNA)+str(Suffix_6)
#print(Ref_6)

Ref_7= str(Prefix_G)+str(TargetDNA)+str(Suffix_7)
#print(Ref_7)

Ref_8= str(Prefix_H)+str(TargetDNA)+str(Suffix_8)
#print((Ref_8))

#reference for PacBio alignment (len is 1202):
reference='ACACGCATGACACACTGTCCACAATTTTCGAAAAAACCCGCTTCGGCGGGTTTTTTTATAGCTAAAAGATTTGACAGCTAGCTCAGTCCTAGGGATTGTGCTAGCGCGTCCGGCGTAGAGGATCGAGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGGTACTAGAGGGCTCGTTGAACACCGTCTCAGGTAAGTATCAGTTGTAAAAAGAGGAGAAATAGTCCATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGATCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAAGGCTCGATCGGTGTGAAAAGTCAGTATCCAGTCGTGTAGTTCTTATTACCTGTCCCCTAGCATAACCCCGCGGGGCCTCTTCGGGGGACTCGCGGGGTTTTTTGCTGAAAGAATTATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGCTGCATTACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGGGCTCGGTACACGCTGTGACTA'
print(len(reference))
test='ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAA'
nonsense_sequence_1 = 'TTTTTTTTTTTTTTTTTTTTTT'
nonsense_sequence_2 = 'NNNNNNNNNN'
max_position_for_mutation_heatmap = 1400
alignment_condition = 1191


def id_mutation(reference_mut,sequence_mut):
    mutationlist = []
    positionlist = [0] * max_position_for_mutation_heatmap
    for base in range(0,len(reference_mut)):
        if reference_mut[base]!=sequence_mut[base]:
            #print(reference[base],'->',sequence_mut[base])
            mutationlist.append(reference_mut[base]+sequence_mut[base])
            #print(base)
            positionlist[base] +=1
    return mutationlist, positionlist

# creates countlists for each mutation type and iterate through sequences in selected fastqfile counting encountered mutations
#can the if statement be extended past the first mutation being detected
#nested for loop maybe?

def countmut(fastqfile):
    startlist=[]
    positionlist = [0] * max_position_for_mutation_heatmap
    countlist={}
    countlist['AT'] = 0
    countlist['AC'] = 0
    countlist['AG'] = 0
    countlist['TA'] = 0
    countlist['TC'] = 0
    countlist['TG'] = 0
    countlist['CT'] = 0
    countlist['CA'] = 0
    countlist['CG'] = 0
    countlist['GT'] = 0
    countlist['GC'] = 0
    countlist['GA'] = 0
    countlist['T-'] = 0
    countlist['A-'] = 0
    countlist['G-'] = 0
    countlist['C-'] = 0
    countlist['-T'] = 0
    countlist['-A'] = 0
    countlist['-G'] = 0
    countlist['-C'] = 0
    sequences_with_additions = 0
    sequences_with_deletions = 0
    sequences_with_substitutions = 0
    nonsense_sequence_count = 0
    unchanged_sequence_count = 0
    total_number_of_mutations_in_file = 0
    number_of_sequences_in_file = 0


    for record in SeqIO.parse(fastqfile, "fastq"):
        startlist.append(str(record.seq)) #or without 'str'
        
    #test = 'ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATCGTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAA'
    print('Total sequences in file are ',len(startlist))
    for sequence in range(0,len(startlist)):
    
        #nonsense sequence test: We try to look for a substring consisting  of 20 Ts or Ns
        #find() function returns -1 if not found
        if (startlist[sequence].find(nonsense_sequence_1) != -1 or startlist[sequence].find(nonsense_sequence_2) != -1) or (len(startlist[sequence]) > 1210) or (len(startlist[sequence]) < 1190):
            nonsense_sequence_count += 1
            continue
        
        #when the resultant sequences is of the same length as reference, we do not want
        #the alignment program inserting gaps to align better
        number_of_sequences_in_file += 1 
        if (len(startlist[sequence]) == len(reference)):
            alignment_score = int(pairwise2.align.globalxs(reference, startlist[sequence], -1.0, -0.1, score_only = True))
   #     print(alignment_score)
            if(alignment_score != len(reference) and alignment_score > alignment_condition):
                #print('Sequence is ', sequence)
                sequences_with_substitutions += 1
                mutationlist, temppositionlist = id_mutation(reference, startlist[sequence])
                positionlist = np.add(positionlist,temppositionlist)
                for index in range(0,len(mutationlist)):
                    countlist[mutationlist[index]] +=1
            elif (alignment_score == len(reference)):
                unchanged_sequence_count += 1
                
        elif (len(startlist[sequence]) < len(reference)):
            alignment_score2 = int(pairwise2.align.globalxs(reference, startlist[sequence], -1.0, -0.1, score_only = True))
   #     print(alignment_score)
            if (alignment_score2 != len(reference) and alignment_score2 > alignment_condition):
                
                sequences_with_deletions += 1
                alignment_result = pairwise2.align.globalxs(reference, startlist[sequence], -1.0, -0.1, one_alignment_only=True)
                #print(alignment_result)
                                
                first_alignment = alignment_result.pop()
                new_reference = first_alignment[0]
                new_sequence = first_alignment[1]
                #print('New Reference is ',new_reference)
                #print('New Sequence is ',new_sequence)
                mutationlist, temppositionlist = id_mutation(new_reference,new_sequence)
                positionlist = np.add(positionlist,temppositionlist)
                for index in range(0,len(mutationlist)):
                    countlist[mutationlist[index]] +=1    
            
        elif (len(startlist[sequence]) > len(reference)) :
            alignment_score3 = int(pairwise2.align.globalxs(reference, startlist[sequence], -1.0, -0.1, score_only = True))
   #     print(alignment_score)
            if (alignment_score3 != len(reference) and alignment_score3 > alignment_condition):
                
                sequences_with_additions += 1
                alignment_result = pairwise2.align.globalxs(reference, startlist[sequence], -1.0, -0.1, one_alignment_only=True)
                #print(alignment_result)
                                                
                first_alignment = alignment_result.pop()
                new_reference = first_alignment[0]
                new_sequence = first_alignment[1]
                # print('New Reference is ',new_reference)
                # print('New Sequence is ',new_sequence)
                mutationlist, temppositionlist = id_mutation(new_reference,new_sequence)
                positionlist = np.add(positionlist,temppositionlist)
                for index in range(0,len(mutationlist)):
                    countlist[mutationlist[index]] +=1
            
    total_number_of_mutations_in_file = np.sum(positionlist)
    avg_mutation_per_sequence = total_number_of_mutations_in_file / number_of_sequences_in_file
    print('Additions in ',sequences_with_additions,' sequences, deletions in ',sequences_with_deletions,' sequences & substitutions in ',sequences_with_substitutions)
    print('Number of nonsense sequences in file are ',nonsense_sequence_count)
    print('Number of unchanged sequences in file are ', unchanged_sequence_count)
    print('Number of sequences in file are ', number_of_sequences_in_file)
    print('Avg no. of mutations per sequence are ', avg_mutation_per_sequence)
    
    return(countlist, positionlist, sequences_with_additions, sequences_with_deletions, sequences_with_substitutions, unchanged_sequence_count, number_of_sequences_in_file, avg_mutation_per_sequence)

# print(countmut('WT_short.fastq'))

# iterates through all fastq files in folder 'input_mutations' and creates corresponding csv files with mutation counts into 'output_mutation' folder

for file in os.listdir('input_mutation'):
    if file.endswith(".fastq"):
        print()
        print(file)
        score_distribution, mutationpositions, number_of_additions, number_of_deletions, number_of_substitutions, unchanged_sequences, Total_sequences, avg_mutations_per_seq = (countmut('input_mutation/'+file))
        print(score_distribution)
        with open('output_mutation/'+file[:-5]+'csv', 'w') as f:
            for item in score_distribution:
                f.write("%s\n" % str(item+','+str(score_distribution[item])))
            f.write("%s\n" % str('Number of sequences with additions,'+ str(number_of_additions)))
            f.write("%s\n" % str('Number of sequences with deletions,'+ str(number_of_deletions)))
            f.write("%s\n" % str('Number of sequences with substitutions,'+ str(number_of_substitutions)))
            f.write("%s\n" % str('Number of unchanged sequences,'+ str(unchanged_sequences)))
            f.write("%s\n" % str('Number of Total sequences,'+ str(Total_sequences)))
            f.write("%s\n" % str('Avg. no. of mutations per sequence,'+ str(avg_mutations_per_seq)))
            for position in range(0, len(mutationpositions)):
                f.write("%s\n" % str(str(position) + ',' + str(mutationpositions[position])))
