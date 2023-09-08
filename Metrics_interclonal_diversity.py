import numpy as np
import random
import sys 
from Bio import SeqIO
import re
from ete3 import Tree 
import heapq, operator
import csv
import pandas as pd
from scipy.spatial import distance
import Metrics.py


# calculate metrics for VHC
dict_VHC = {}
for i in range (1, 10):
    file_fasta = 'VHC/VHC_dataset0' + str(i) + '_1_200_sequences.aln.fa'
    file_newick = 'VHC/VHC_dataset0' + str(i)+'_1_200_sequences.aln.fa.nk'
    dict_VHC[file_fasta[4:33]] = list(allMetrics(file_fasta, file_newick))
        
for i in range (10, 41):
    file_fasta = 'VHC/VHC_dataset' + str(i) + '_1_200_sequences.aln.fa'
    file_newick = 'VHC/VHC_dataset' + str(i)+'_1_200_sequences.aln.fa.nk'
    dict_VHC[file_fasta[4:33]] = list(allMetrics(file_fasta, file_newick))
    #print (i, allMetrics(file_fasta, file_newick))

# calculate metrics for MUT
MUT_fasta = []
with open ('MUTnames.txt', 'r') as file:
    for line in file:
        MUT_fasta.append(line.replace('\n', ''))

dict_MUT = {}
for seq in MUT_fasta:
    file_fasta = 'MUT/' + seq
    file_newick = 'MUT/' + seq + '.nk'
    dict_MUT[file_fasta[4:-7]] = list(allMetrics(file_fasta, file_newick))

# calculate metrics for LLC
dict_LLC = {}
for i in range (1, 10):
    file_fasta = 'LLC/LLC_dataset0' + str(i) + '_1_200_sequences.aln.fa'
    file_newick = 'LLC/LLC_dataset0' + str(i)+'_1_200_sequences.aln.fa.nk'
    dict_LLC[file_fasta[4:33]] = list(allMetrics(file_fasta, file_newick))
        
for i in range (10, 57):
    file_fasta = 'LLC/LLC_dataset' + str(i) + '_1_200_sequences.aln.fa'
    file_newick = 'LLC/LLC_dataset' + str(i)+'_1_200_sequences.aln.fa.nk'
    dict_LLC[file_fasta[4:33]] = list(allMetrics(file_fasta, file_newick))


# create csv with metrics
def write_metrics_in_csv(our_dict, filename):    
    for key in our_dict.keys():
        if filename == 'LLC' or filename == 'MUT':
            groupe = 'LLC + MUT'
        else:
            groupe = 'VHC'
        s = [key, groupe]
        for elem in our_dict[key]:
            s.append(elem)
        with open(groupe + '.csv', "a") as p:
            pr = csv.writer(p, delimiter = ";", lineterminator = '\n')
            pr.writerow(s)

# create csv with metrics for PCA (categorie 1- LLC+MUT categorie 2 - VHC)
def write_metrics_in_csv_for_PCA(our_dict, filename):    
    for key in our_dict.keys():
        if filename == 'LLC' or filename == 'MUT':
            name = 'LLC + MUT'
            groupe = 1
        else:
            name = 'VHC'
            groupe = 2
        s = [key, groupe]
        for elem in our_dict[key]:
            s.append(elem)
        with open(name + '.csv', "a") as p:
            pr = csv.writer(p, delimiter = ";", lineterminator = '\n')
            pr.writerow(s)


write_metrics_in_csv(dict_LLC, 'LLC')
write_metrics_in_csv(dict_MUT, 'MUT')
write_metrics_in_csv(dict_VHC, 'VHC')


write_metrics_in_csv_for_PCA(dict_LLC, 'LLC')
write_metrics_in_csv_for_PCA(dict_MUT, 'MUT')
write_metrics_in_csv_for_PCA(dict_VHC, 'VHC')


# datasets LLC+MUT and VHC
df_LLC_MUT = pd.read_csv('LLC + MUT.csv', delimiter=';', names=['Newick file', 'Dataset', 'Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree'])
df_VHC = pd.read_csv('VHC.csv', delimiter=';', names=['Newick file', 'Dataset', 'Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree'])


# percentage of H1, H2, H3 == 0
def h_count_zero(metric):
    print ('LLC + MUT ', metric, ' : ', round ((df_LLC_MUT[metric] == 0).sum(axis=0) / len(df_LLC_MUT), 1) * 100, '%\n')
    print ('VHC       ', metric, ' : ', (df_VHC[metric] == 0).sum(axis=0) / len(df_VHC) * 100, '%\n')

# percentage of D == i   
def d_count(metric, i):
    print ('LLC + MUT ', metric, ' = ', i, ' : ', round ((df_LLC_MUT[metric] == i).sum(axis=0) / len(df_LLC_MUT), 1) * 100, '%\n')
    #print ('VHC       ', metric, ' = ', i, ' : ', round ((df_VHC[metric] == i).sum(axis=0) / len(df_VHC), 1) * 100, '%\n')
    
metricsH = ['H1', 'H2', 'H3']
metricsD = [ 'D1', 'D2', 'D3']

for m in metricsH:
    h_count(m)

for m in metricsD:
    for i in range(1, 6):
        d_count(m, i)


pd.set_option('display.float_format', lambda x: '%.1f' % x)
df_VHC.describe()
df_LLC_MUT.describe()

df_all = df_LLC_MUT.append(df_VHC, ignore_index=True)

# boxplots for metrics
def made_boxplot(metric):
    boxplot = df_all.boxplot(by = 'Dataset', column =[metric], grid = True, rot=45, fontsize=10)

    boxplot.get_figure().gca().set_title("") 
    boxplot.get_figure().suptitle(metric)
    boxplot.plot()

    plt.show()
    
metrics = ['Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree']
for metric in metrics:
    made_boxplot(metric)

