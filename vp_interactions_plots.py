#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:51:19 2023
Description:

The aim of this script is to provide statistical insights on viral proteins clusters and their host-interactions.
Non standard Functions: NA

Procedure:

      
Version 1.00
Date: July 4th 2023.
Author: Tamim AlMurad
"""
import pandas as pd
import matplotlib.pyplot as plt

#Load files into dataframes and remove unnamed columns.
vpHostInteractions = pd.read_csv('all_clusters_host_tax.tsv',sep='\t',engine='python')
hetroClustersMembers = pd.read_csv('hetro_cell_clusters.tsv',sep='\t',engine='python')
hetroClusters = pd.read_csv('hetro_cell_cluster_analysis.tsv',sep='\t',engine='python')
virusHostDB = pd.read_csv('virushostdb.tsv',sep='\t',engine='python')
vpHostInteractions.drop('Unnamed: 0',axis ='columns', inplace = True)
hetroClustersMembers.drop('Unnamed: 0',axis ='columns', inplace = True)
hetroClusters.drop('Unnamed: 0',axis ='columns', inplace = True)
#%%
# Plot viral_entries_host_distribution among life domains.
vpHostInteractions['host_cell'].value_counts().plot(kind = 'bar',rot = 0 ,grid=True,log=True,
                                                    title='Viral Entries counts According to host (log)',
                                                    xlabel='Host Life Domain',ylabel='# of Interactions').figure.savefig('viral_entries_host_distribution.png')
plt.close()

#%%
# Plot cluster a histogram for cluster interactions.
hetroClusters['cluster_host_size'].plot(kind = 'hist', rot =0, grid = True,log=True,bins =50, 
                                        title='Histogram of Number of Host-Interactions among Hetro Clusters (log)',
                                        xlabel='#of Interactions',ylabel='Frequency').figure.savefig('hetro_clus_interactions_hist.png')
plt.close()
#%%
# Plot a histogram of the heterogeneity level of clusters.
hetroClusters['hetro_strength'].plot(kind = 'hist', rot =0, grid = True,bins =50, 
                                        title='Host Hetrogeniety Histogram',
                                        xlabel='Hetro Strength % (percentage of 2nd Dominant Host)',ylabel='Frequency').figure.savefig('Host_Hetrogeniety_Hist.png')
plt.close()
#%%
#Print some statistics
print('Number of hosts from available viral protein data: %d'%len(set(list(vpHostInteractions['host tax id']))))
print('Number of virus species from available viral protein data: %d'%len(set(list(vpHostInteractions['viral_tax_ID']))))
print('Number of hosts species according to Virus-Host DB: %d'%len(set(list(virusHostDB['host tax id']))))
print('Number of viral species according to Virus-Host DB: %d'%len(set(list(virusHostDB  ['virus tax id']))))



