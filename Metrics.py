import numpy as np
import random
import sys 
from Bio import SeqIO
import re
from ete3 import Tree 
import heapq, operator
import csv

#====================================================================
# function to read Fasta file
# input - fasta file
# output - array of sequnces names and array of sequences
#====================================================================

def readFasta_for_metrics(fastaFile):
    labels = []
    arraySeqs = []
    with open (fastaFile, 'r') as file:
        for line in file:
            if (line[0] == '>'):
                name = line[1:].replace ('\n', '')
                name = re.sub('@[0-9]*', '', name)
                labels.append(name)
            else:
                arraySeqs.append(line.replace('\n', ''))
    return labels, arraySeqs

#====================================================================
# function to calculate Hamming distance 
# input - two sequences
# output - Hamming distance 
#====================================================================

def hamming_distance_for_metrics(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


#====================================================================
# function to calculate:
# The number of branches 
# Sum of all branch lengths (PD) 
# Sum of all branch lengths divided by the number of clonotypes of the lineage (avPD) 
# input1 - ete3 tree 
# input2 - array of sequnces names 
# input3 - array of sequences
# output - number of branches, PD, avPD
#====================================================================

def metricsPD(tree, labels, arraySeqs):
    num_of_branches = 0
    PD = 0
    avPD = 0.
    
    for parent in tree.traverse("preorder"):
        children = parent.children
        for child in children:
            if child.name != 'naive':
                if child.name == '':
                    remember_parent = parent
                else:
                    if parent.name == '':
                        parent_name = remember_parent.name
                    else:
                        parent_name = parent.name
                    child_name = child.name
                    index_p = labels.index(parent_name)
                    index_c = labels.index(child_name)
                    ham = hamming_distance_for_metrics(arraySeqs[index_p], arraySeqs[index_c])
                    num_of_branches +=1
                    PD += ham
                    #print (parent_name, child_name, ham)
    avPD = PD / num_of_branches
    avPD = round(avPD , 2)
    return num_of_branches, PD, avPD


#====================================================================
# function to choose 3 most abundant sequences
# input - fasta file
# output - array of names for 3 most abundant sequences
#====================================================================

def get_3_most_abundance (file) :
    abundance = {}
    #abundance["total"] = 0
    with open(file) as f:
        for line in f :
            if(line[0] == ">"):
                seq = line.split("\n")[0]
                seq_info = seq.split("@")
                seq_name = seq_info[0][1:]
                if(len(seq_info) == 2) :
                    number_seq = int(seq_info[1])
                    abundance[seq_name] = number_seq
                    #abundance["total"] += number_seq
                else :
                    if seq_name in abundance :
                        abundance[seq_name] += 1
                        #abundance["total"] += 1
                    else :
                        abundance[seq_name] = 1
                        #abundance["total"] += 1
    abud3 = list(zip(*heapq.nlargest(3, abundance.items(), key=operator.itemgetter(1))))[0]
    return abud3


#====================================================================
# function that return:
# The path from nodeName to root 
# input1 - nodeName : name of first node of path,
# input2 - ete3 tree
# output - list of nodes in the path
#====================================================================

def pathToRoot(nodeName, tree):
	D = tree&nodeName
	# Get the path from B to the root
	node = D
	path = []
	while node.up:
		nn = node.name
		if node.name == '' or node.name.isdigit():
			nn = 'none'
		path.append(nn)
		node = node.up
	if path and path[len(path)-1] == 'none':
		path[len(path)-1] = 'naive'
	return path

#====================================================================
# function that return:
# The path from nodeName to nodeFin 
# input - nodeName : name of 1 node of path, nodeFin : name of last node in path (won't print in path), tree
# output - list of nodes in the path
#====================================================================

def pathToNode(nodeName, nodeFin, tree):
	D = tree&nodeName
	node = D
	path = []
	while node.up and node.name != nodeFin:
		nn = node.name
		if node.name == '' or node.name.isdigit():
			nn = 'none'
		path.append(nn)
		node = node.up
	return path

#====================================================================
# function that return:
# The height of a node : the number of edges on the longest path from the node to a leaf. A leaf node will have a height of 0.
# input - name of node that height we want to calculate, tree
# output1 - height
# output2 - flag = 1 if we had unoserved node on path 
#====================================================================

def height(nodeName, tree):
	leafs = []

	for node in tree.traverse("preorder"):
			if node.is_leaf() == True:
				leafs.append(node.name)

	path_depth = []			
	if nodeName in leafs:
		h = 0
	else:
		h = -1
		for leaf in leafs:
			path = pathToNode(leaf, nodeName, tree)
			len_path = len(path)
			if len_path > h and ('naive' not in path):
				h = len_path
				path_depth = path # for print path and how we received this height
                
	if 'none' in path_depth: #check for unobserved node
		flag = 1
	else:
		flag = 0

	return h, flag

#====================================================================
# function that return:
# The depth of a node : the number of edges from the node to a root. 
# input - name of node that depth we want to calculate, tree
# output1 - depth
# output2 - flag = 1 if we had unoserved node on path 
#====================================================================

def depth(nodeName, tree):
    path_depth = pathToRoot(nodeName, tree)
    d = len(path_depth) - 1
    
    if 'none' in path_depth: #check for unobserved node
        flag = 1
    else:
        flag = 0
        
    return d, flag

#====================================================================
# function that return:
# The overall depth of tree. 
# input - tree
# output - overall depth
#====================================================================

def TreeSize(tree):

	leafs = []

	for node in tree.traverse("preorder"):
			if node.is_leaf() == True:
				leafs.append(node.name)

	h = -1
	for leaf in leafs:
		path = pathToRoot(leaf, tree)
		len_path = len(path)
		if len_path > h:
			h = len_path

	return h


#===================================================================================================
# function that return:
# Array of depth and height 
# input1 - fasta file
# input2 - tree
# outputs - height array, depth array, height flags, depth flags, names of 3 most abundant sequences
#====================================================================================================

def metricsDH(fasta, tree):
    
    array3 = get_3_most_abundance(fasta)
    
    d_list = []
    d_list_flag = []
    h_list = []
    h_list_flag = []
    
    for node in array3:
        
        d, flag_d = depth(node, tree)
        d_list.append(d)
        d_list_flag.append(flag_d)
        
        h, flag_h = height(node, tree)
        h_list.append(h)
        h_list_flag.append(flag_h)
        
    return h_list, d_list, h_list_flag, d_list_flag, array3


#====================================================================
# function that return:
# all metrics 
# input1 - fasta file
# input2 - Newick file
# output - all metrics
#====================================================================

def allMetrics(input_fasta, input_newick):    
    labels1, arraySeqs1  = readFasta_for_metrics(input_fasta)
    with open (input_newick, 'r') as file:
        for line in file :
            t = Tree(line, format = 1)
    Num1, PD1, avPD1 = metricsPD(t, labels1, arraySeqs1)
    h1, d1, h1_flags, d1_flags, most_abund3 = metricsDH(input_fasta, t)
    size = TreeSize(t)
    return Num1, PD1, avPD1, h1[0], h1[1], h1[2], d1[0], d1[1], d1[2], size

#====================================================================
# function that write all metrics in file
# input1 - fasta file
# input2 - tree
# input3 - output file
#====================================================================

def all_metrics_write_in_file(input_fasta, tree, output_file):    
    labels1, arraySeqs1  = readFasta_for_metrics(input_fasta)
    t = tree
    Num1, PD1, avPD1 = metricsPD(t, labels1, arraySeqs1)
    h1, d1, h1_flags, d1_flags, most_abund3 = metricsDH(input_fasta, t)
    size = TreeSize(t)
    
    filename_txt = output_file.replace('.abRT.nk', '') + '_Metrics.txt'
    f = open(filename_txt, 'w')
    f.writelines ('\n SizeTree = ' + str(size) + '\n')
    f.write('\n Number of branches : '+ str(Num1) + '\n')
    f.write('\n Sum of all branches length (PD) : '+ str(PD1) + '\n')
    f.write('\n Sum of all branch length divided by the number of clonotypes of the lineage (avPD) : '+ str(avPD1) + '\n')
    
    f.writelines ('\n#==============================================================================================#\n\n')
    f.writelines (' Height\n')
    
    if h1_flags[0] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
        
    f.writelines('\n The most abundant clonotype - ' + most_abund3[0] + ' : height = ' + str(h1[0]) + s_un)
    
    if h1_flags[1] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
    
    f.writelines('\n The second most abundant clonotype - ' + most_abund3[1] + ' : height = ' + str(h1[1]) + s_un)
    
    if h1_flags[2] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
    
    f.writelines('\n The third most abundant clonotype - ' + most_abund3[2] + ' : height = ' + str(h1[2]) + s_un)
    
        
    f.writelines ('\n#==============================================================================================#\n\n')
    f.writelines (' Depth\n')
    
    if d1_flags[0] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
        
    f.writelines('\n The most abundant clonotype - ' + most_abund3[0] + ' : depth = ' + str(d1[0]) + s_un)
    
    if d1_flags[1] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
        
    f.writelines('\n The second most abundant clonotype - ' + most_abund3[1] + ' : depth = ' + str(d1[1]) + s_un)
    
    if d1_flags[2] == 1:
        s_un = ' (Unobserved node!) \n'
    else:
        s_un = ' \n'
        
    f.writelines('\n The third most abundant clonotype - ' + most_abund3[2] + ' : depth = ' + str(d1[2]) + s_un)
    
        
    f.close()


