#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:20:21 2023
Description:

In this script, the previously produced hetro_cell_clusters is loaded as a dataframe.
Then a new dataframe is produced where a row for each cluster is found.
The row contains columns that are : accn, cluster size, hetro level, and 4 more columns to
indicate the hosts cell types. Additionally, a column showing the hetro strength, the dominant 
life domain and the 2nd dominant.

Non standard Functions: NA

Procedure:
    1. Read the hetro clusters file into a dataframe.
    2. Determine host size of all cluster by counting the number of rows where a rep is repeated.
    3. The results from step 2 are stored in a dataframe with rows containing cluster rep accn and its host size.
    3. Initialize a new dataframe that only contains one row for each cluster.
    4. Iterated over clusters, determine related hosts, dominant and 2nd dominant life domains and its hetrogeniety 
        strength. The later is defined as the percentage of number of hosts from 2nd dominant life domain infected
        by proteins from viruses in the cluster.
        Example: Cluster 'A', has member proteins of viruses with 100 host-interactions.
        50 of these host-interactions are prokaryotes, 40  eukaryotes and 10 are Archaea. Then the following
        applies: Dominant host: Prokaryotes, 2nd Dominant host: Eukaryotes, Hetro Strength = 40%
    7. Note: a protein can belong to a virus that infects several hosts. e.g. mammalian viruses. A host can have 
        be infected by several viruses with proteins in the same cluster. Host-interactions are what we are considering.
Input: Hetro clusters tsve file.
Output: Hetro clusters hetro analysis (one row per cluster)
      
Version 1.00
Date: July 3rd 2023.
Author: Tamim AlMurad

"""

import pandas as pd

#Read hetro host clusters in terms of domains of life (cell types).
hetroCellClusters = pd.read_csv('hetro_cell_clusters.tsv',sep='\t',engine='python')
hetroCellClusters.drop('Unnamed: 0',axis ='columns', inplace = True)
#Counting the number of occurences of the cluster rep. This will show how many host-interactions a cluster have.
clusters = hetroCellClusters['cluster_rep_accn'].value_counts()
clusters = pd.DataFrame(clusters)
clusters = clusters.reset_index()

#Hetro clusters analysis dataframe
clustersHetro = pd.DataFrame(columns=['cluster_rep_accn','cluster_host_size','hetro_level','Bacteria', 'Archaea','Eukaryota' ,'Viruses',
                                      'hetro_strength','dominant','2nd_dominant'])
for index, row in clusters.iterrows():
    
    accn = row['cluster_rep_accn']
    #cluster size
    count = row['count']
    #Get cluster members info
    clusterData = hetroCellClusters[hetroCellClusters['cluster_rep_accn']==accn]
    
    #Life domains affected by the cluster.
    affectedCells = list(clusterData['host_cell']) 
    #Count the number of interactions of the cluster with every life domain.
    hetroCounts = [affectedCells.count('Bacteria'),affectedCells.count('Archaea'),
                          affectedCells.count('Eukaryota'),affectedCells.count('Viruses')]
    
    #Sort the number of interactions.
    hetroCounts.sort()
    #Most affected life domain.
    clusterDominant = clusterData['host_cell'].mode()[0]
    #2nd most affected life domain.
    cluster2ndDominant = clusterData['host_cell'].value_counts().index[1]
    #Create a row for the cluster and add it to the clusterHetro dataframe.
    newRow = {'cluster_rep_accn':accn,'cluster_host_size':count,'hetro_level':len(set(affectedCells)), 
              'Bacteria':affectedCells.count('Bacteria'), 'Archaea': affectedCells.count('Archaea'), 
              'Eukaryota': affectedCells.count('Eukaryota'), 'Viruses':affectedCells.count('Viruses'),
              'hetro_strength': 100*(hetroCounts[-2]/count),'dominant':clusterDominant,'2nd_dominant': cluster2ndDominant  }
    clustersHetro.loc[len(clustersHetro)] = newRow

# Write the results into a tsv file
clustersHetro.to_csv('hetro_cell_cluster_analysis.tsv', sep='\t')

#%%
# This is a function to get cluster members of a certain clusters and its information.
def get_hetro_cluster (accn, hetroCellClusters):
    return hetroCellClusters[hetroCellClusters['cluster_rep_accn']==accn]

#Examples
clusterToStudy1 = get_hetro_cluster('4GQH_A', hetroCellClusters)
clusterToStudy2 = get_hetro_cluster('YP_002265459.1', hetroCellClusters)









    

















