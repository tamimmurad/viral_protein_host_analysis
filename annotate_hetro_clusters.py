#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:04:40 2023
Description:

This script takes as an input the a tsv file containing the hetro host clusters(sequence based) information and 
add two extra columns. The first column is the cluster size (number of members of in the cluster) and 
the second column is the NCBI annotation fetched directly by using Entrez. The results is written
to a new tsv file.  

Non standard Functions: NA

Procedure:
    
    1. Read the clusters file into a dataframe.
    2. Loop over the rows of the dataframe and determine the size of each cluster and 
        the annotation of the each member. Add the the results to ne columns.
    3. Sort the dataframe in a descending form.
    4. Write the resulting dataframe to a file.
    
Input: Hetro clusters tsve file and Entrez email address to be specified in the script.
Output: Tsv file of hetro host clusters with annotations.

Version 1.00
Date: July 3rd 2023.
Author: Tamim AlMurad

"""

import pandas as pd
from Bio import Entrez

# Read the hetro host clusters file.
clusters = pd.read_csv('hetro_cell_clusters.tsv',sep='\t').drop('Unnamed: 0',axis=1)

# debugging use
i = 0

Entrez.email = "ta4818al-s@student.lu.se"
Entrez.api_key = '63fcf49b0aa53822a393229e8373b128e608'

#Loop over the dataframe row by row.
for index, row in clusters.iterrows():
    #Debugging use to check progress.
    i+=1
    
    #cluster rep of the row.
    rep = row['cluster_rep_accn']
    
    #Determin the cluster size by counting the number of occureneces of the cluster rep.
    clusterSize = clusters['cluster_rep_accn'].value_counts()[rep]
    #Write the cluster size into a new column.
    clusters.at[index,'cluster_host_size'] = clusterSize
    
    #Take the member accn number.
    member = row['member_accn']
    
    #Fetch the entry of the cluster member.
    try:
        handle = Entrez.efetch(db="protein", id=member,rettype ='xml')
        record = Entrez.read(handle)
    
        #Get the annotation.
        annotation = record[0]['GBSeq_definition']
        #Write the annotation to the annotation column.
        clusters.at[index,'annotation'] = annotation
        #Show the progress.
        print('reaching  '+str(i))
    
    except:
        #In case the entry is not found.
        print(member + '   not available at index ' +str(i) )

# Sort the results so larger clusters are on top.
clusters.sort_values(by = ['cluster_host_size'], ascending = False , inplace = True)

# Write the output to a new tsv file.
clusters.to_csv('hetro_seq_clusters_annotations.tsv',sep='\t')
