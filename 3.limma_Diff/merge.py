# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 09:39:08 2023

@author: Administrator
"""

import os


files = os.listdir(".")
files = [ f for f in files if f.startswith("immune") and f.endswith("_diff.xls")]

diff_compare = {}
compares = {}
genes = {}
all_result = []
for f in files:
    with open(f) as fin:
        h = fin.readline().strip().split("\t")
        compare = f.split("_")[2]
        compares[compare] = 1
        if compare not in diff_compare:
            diff_compare[compare] = {}
        for line in fin:
            content = line.strip().split("\t")
            ncontent = content+[compare]
            all_result.append(ncontent)
            genes[content[0]] = 1
            diff_compare[compare][content[0]] = content[1]
all_result = [h+["compare"]]+all_result

with open("TCGA-STAD.htseq_fpkm.tsv_format.txt_tumor_mRNA_exp.txt") as efin,\
    open("all_cluster_compare_result.txt","w") as rout1,\
    open("all_cluster_comapre_matrix.txt","w") as rout2,\
    open("TCGA-STAD_fpkm_cluster_diff.txt","w") as rout3:
        eh = efin.readline().strip()
        print(eh,file=rout3)
        for line in efin:
            content = line.strip().split("\t")
            if content[0] in genes:
                print(line.strip(),file=rout3)
        for result in all_result:
            print("\t".join(result),file=rout1)
        header = ["gene"] + [c for c in compares]
        print("\t".join(header),file=rout2)
        for g in genes:
            v = [g]
            for c in compares:
                v.append(diff_compare[c].get(g,"0"))
            vcount = len([m for m in v[1:] if m!="0"])
            print("\t".join(v),file=rout2)
            
            
            
            