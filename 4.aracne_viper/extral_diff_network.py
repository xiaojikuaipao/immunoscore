# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 13:34:43 2023

@author: Administrator
"""
import pandas as pd

dfin = pd.read_table("TCGA-STAD_fpkm_cluster_diff.txt",delimiter="\t",index_col=0)
genes = list(dfin.index)

nfin = pd.read_table("network.txt",delimiter="\t")
filter_nfin = nfin[nfin["Target"].isin(genes)]
filter_nfin.to_csv("diff_network.txt",sep="\t",index=False)

filter_nfin.to_csv("diff_network_sub.txt",sep="\t",index=False,header=False)
