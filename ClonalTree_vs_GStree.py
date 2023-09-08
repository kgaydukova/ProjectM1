import numpy as np
import random
import sys 
import re
from ete3 import Tree 
import csv
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt

pd.set_option('display.float_format', lambda x: '%.2f' % x)


# generate bash commands to launch clonaltree
def generate_command_for_bash(n, j):
    for i in range (1, n + 1):
        print ('python src/clonalTree.py  -i Data/Simulations/' + str(j) + '/' + str(j) +'_' + str(i) + '.fasta -o Data/Simulations/' + str(j) + '/' + str(j) +'_' + str(i) + '.clonalTree.nk -a 0 -r 0 -t 0')
nn = [50, 60, 70, 80, 90, 100, 200]
for number in nn:
    generate_command_for_bash(10, number)
    print ('fin')

generate_command_for_bash(10, 150)

dict = {} # {name of sequence : [all metrics]}

#function for metrics of ClonalTree and GCtree sequences from 30_, 40_, 50_, 60_, 70_, 80_, 90_, 100_
def run2arbres(j):    
    for i in range(1,11):
        fastaFile = 'Simulations/' + str(j) + '/' + str(j) + '_'  + str(i) + '.fasta'
        nClonalTree = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i) +'.clonalTree.nk'
        nGT = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i)+ '.GT.naive.nk'
        if j < 100:
            dict[nClonalTree[15:-3]] = list(allMetrics(fastaFile, nClonalTree))
            dict[nGT[15:-3]] = list(allMetrics(fastaFile, nGT))
        else:
            dict[nClonalTree[16:-3]] = list(allMetrics(fastaFile, nClonalTree))
            dict[nGT[16:-3]] = list(allMetrics(fastaFile, nGT))
            
#function for metrics of ClonalTree and GCtree sequences from 200_           
def run2arbres200(j):    
    for i in range(1,11):
        if i != 2:
            fastaFile = 'Simulations/' + str(j) + '/' + str(j) + '_'  + str(i) + '.fasta'
            nClonalTree = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i) +'.clonalTree.nk'
            nGT = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i)+ '.GT.naive.nk'
            dict[nClonalTree[16:-3]] = list(allMetrics(fastaFile, nClonalTree))
            dict[nGT[16:-3]] = list(allMetrics(fastaFile, nGT))

#function for metrics of ClonalTree and GCtree sequences from 150_  
def run2arbres150(j):    
    for i in [1, 2, 10]:
        fastaFile = 'Simulations/' + str(j) + '/' + str(j) + '_'  + str(i) + '.fasta'
        nClonalTree = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i) +'.clonalTree.nk'
        nGT = 'Simulations/' + str(j) + '/' + str(j) + '_' + str(i)+ '.GT.naive.nk'
        dict[nClonalTree[16:-3]] = list(allMetrics(fastaFile, nClonalTree))
        dict[nGT[16:-3]] = list(allMetrics(fastaFile, nGT))


numbers = [30, 40, 50, 60, 70, 80, 90, 100] 
for j in numbers:
    run2arbres(j)
    
run2arbres150 (150)   
run2arbres200 (200)


# function for score (Euclidean distance)
def score(a, b):
    mean = (a + b) / 2
    sd = np.sum(np.sqrt(np.square(a - mean) + np.square(b - mean)))
    sd = round (sd, 1)
    ed = np.sqrt(np.sum(np.square(a - b)))
    ed = round(ed, 3)
    
    met = []
    n = 0
    a_b = a - b 
    metrics = ['Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree']
    
    for i in range(len(metrics)):
        if a_b[i] != 0:
            met.append(metrics[i])
            
    return [sd, ed, met]


dict_euclidean = {} # {name of sequence : score}
dict_euclidean_m = {} # {name of sequence : array of absolute value difference for all metrics}
dict_non_zero = {}  #{name of sequence : score, [metrics with difference]}

for key in dict.keys():
    if 'clonalTree' in key:
        a = np.asarray(dict[key])
        ab = key
        ab = re.sub('.[a-z]*[A-Z][a-z]*', '', ab)
    if 'GT.naive' in key:
        b = np.asarray(dict[key])
        a_b = a - b
        if a_b.any() != 0:
            sc = score(a, b)
            dict_euclidean[ab] = sc
            #dict_euclidean[ab] = np.abs(a - b)
        else:
            dict_euclidean[ab] = [0, 0]
        if dict_euclidean[ab] != [0, 0]:
            dict_euclidean_m[ab] = np.absolute(a - b)
            dict_non_zero[ab] = dict_euclidean[ab]



dict_with_scores = {} # {name of sequence : score, groupe 30-100, 150, 200}
d_keys = []
d_score_sd = []
d_score_ed =[]
d_groupe = [] 
    
for key in dict_euclidean.keys():
    if key[0:2] == '10':
        groupe = 100
    elif key[0:2] == '15':
        groupe = 150
    elif key[0:2] == '20':
        groupe = 200
    else :
        groupe = int (key[0:2])
    dict_with_scores[key] = [dict_euclidean[key][0], dict_euclidean[key][1], groupe]
    d_keys.append(key)
    d_score_sd.append(dict_euclidean[key][0])
    d_score_ed.append (dict_euclidean[key][1])
    d_groupe.append(groupe)


# creating a Dataframe object from a list 
df_scores = pd.DataFrame({'Name':d_keys , 'Score_sd': d_score_sd,'Score_ed': d_score_ed,  'Groupe' : d_groupe})


# percentage of different scores 
print ((df_scores['Score_ed'] > 6).sum(axis=0) /92 * 100)
print ((df_scores['Score_ed'] > 4).sum(axis=0)  / 92 * 100)
print ((df_scores['Score_ed'] > 2).sum(axis=0) /92 * 100)
print ((df_scores['Score_ed'] > 0).sum(axis=0)  / 92 * 100)
print ((df_scores['Score_ed'] == 0).sum(axis=0)  / 92 * 100)


# box-plot for scores for different groups
boxplot = df_scores.boxplot(by = 'Groupe', column =['Score_ed'], grid = True, rot=45, fontsize=15)
boxplot.get_figure().gca().set_title("") 
boxplot.get_figure().suptitle('Scores for simulated lineages trees')
boxplot.get_figure().gca().set_ylabel("Score") 
boxplot.get_figure().gca().set_xlabel("Number of sequences in repertoire")
boxplot.plot()

plt.show()


# write in csv diffrence between ClonalTree and GStree trees
for key in dict_euclidean_m.keys():
    s = [key]
    for elem in dict_euclidean_m[key]:
        s.append(elem)
    with open("2arbresdifferent.csv", "a") as p:
        pr = csv.writer(p, delimiter = ";", lineterminator = '\n')
        pr.writerow(s)


df2 = pd.read_csv('2arbresdifferent.csv', delimiter=';', names=['Newick file', 'Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree'])


# count how many times metrics were different
metrics = ['Number of branches', 'PD', 'avPD', 'H1', 'H2', 'H3', 'D1', 'D2', 'D3', 'SizeTree']
for metric in metrics:
    #print (df2[metric].value_counts())
    print ((df2[metric] != 0).sum(axis=0), '\n')


# statistic about metrics difference 
df3 = df2
df3 = df3.replace(0, np.NaN)
df3.describe()


