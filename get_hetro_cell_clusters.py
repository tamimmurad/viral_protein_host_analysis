#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:

In this script, the virus_host DB is used to fetch the taxonomy of the viral protein sequence
and its host. The code uses ACCN to find taxID then use this id to fetch taxonomy of the host.
Additionally, the script produces three output files.
 1- Clusters with members taxa and host taxa.
 2- One for clusters with members proteins of viruses having host-interactions with more than one life domain.
 3- One for clusters with hosts belong to more than one species.
 
Non standard Functions: NA

Procedure:
     
     1. Read the full viral taxID file keyed by viral protein accn into a dataframe.
     2. Read the clusters file into a dataframe.
     3. Read the virus-host-db into a dataframe.
     4. First merge the clusters file and the taxID file. So for every cluster member have associated taxID.  
     5. Merge the result of 4 with the virus-host-db so that for every cluster member we have the full tax
         of the virus and its host(s).
     6. Loop over clusters and extract clusters with members interacting with more than one life domains 
         and seperately extract clusters with more than one host.
     
 Input: taxIDs file keyed by protein accns, protein clusters (reps and members) file, virus-host-DB file.  
 Output: 3 TSV files as mention in the description.

Version 1.00
Date: July 7th 2023.
Author: Tamim AlMurad
"""
import pandas as pd


#Read all viral accns with taxid of viruses (downloaded & curated from NCBI on May31st 2023)
viralAccnsTaxa =pd.read_csv('accn_tax_corrected.tsv',sep='\t',engine='python',names=['ACCN','viral_tax_ID'])
#Read clusters members data into a dataframe. The dataframe will contain only two columns with cluster rep and member respectively.
vpClusters = pd.read_csv('vp_clusters_cluster.tsv',sep='\t', engine='python',names=['cluster_rep_accn','member_accn'])
#Read virus_host_db as a dataframe with selected columns.
fields = ['virus tax id','virus name','virus lineage','DISEASE','host tax id','host name','host lineage']
virusHostDB= pd.read_csv('virushostdb.tsv',sep='\t',engine='python',usecols=fields)
#Change column name of tax id to be matching used name throughout the code..
virusHostDB = virusHostDB.rename(columns = {'virus tax id': 'viral_tax_ID'})
#A new dataframe with by merging the clusters dataframe with the tax data frame. 
# Every row shows the cluster member taxID.
vpClustersTax = pd.merge(vpClusters, viralAccnsTaxa,left_on='member_accn',right_on='ACCN')
#Drop extra column.
vpClustersTax = vpClustersTax.drop(columns=['ACCN'])
# Final dataframe to contain for every member the tax data of its virus and the tax data of the host.
vpClustersHostTax = pd.merge(vpClustersTax, virusHostDB,on='viral_tax_ID')
vpClustersHostTax[['host_cell', 'host_lineage']] = vpClustersHostTax['host lineage'].str.split(';',n=1,expand =True)



#Make a list of all clusters' reps.
allReps = list(dict.fromkeys(list(vpClustersHostTax['cluster_rep_accn'])))

#Initialize DFs to be populated
hetroCellClusters = pd.DataFrame()
hetroHostClusters = pd.DataFrame()

#Iterate over clusters and get clusters data of clusters having members with hosts from different cell types.
#Also another dataframe will contain hetro clusters form more than one host.
for rep in allReps:
    #Get the cluster data by its rep accn.
    clusterData = vpClustersHostTax[vpClustersHostTax['cluster_rep_accn']==rep]
    #remove null entries where host is not identified. This happens in uncultured viral sequences.
    clusterData = clusterData[clusterData['host_cell'].isnull()==False]
    #Get clusters with hosts distributed over more than one life domain.
    if len(dict.fromkeys(list(clusterData['host_cell']))) > 1:
        hetroCellClusters = pd.concat([hetroCellClusters,clusterData],ignore_index=True)
    #Get clusters with hosts distributed over more more than one host.
    if len(dict.fromkeys(list(clusterData['host tax id']))) > 1:
        hetroHostClusters = pd.concat([hetroHostClusters,clusterData],ignore_index=True)
         


#Write clusters host interactions for all and then for Hetro clusters into files.
vpClustersHostTax.to_csv('all_clusters_host_tax.tsv', sep = '\t')
hetroCellClusters.to_csv('hetro_cell_clusters.tsv', sep='\t')
hetroHostClusters.to_csv('hetro_host_clusters.tsv', sep='\t')
