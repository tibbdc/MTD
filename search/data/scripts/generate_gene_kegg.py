#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-05-10 06:09:43
LastEditors: wangruoyu wangry@tib.cas.cn
LastEditTime: 2023-09-06 05:54:37
Description: file content
FilePath: /mtd_cdk/lambda/mtd/data/0/generate_gene_kegg.py
'''
import os
import sys
import json
import pandas as pd

species_index = sys.argv[1]

with open(f"./{species_index}/gene_information.json",'r') as f:
    genes = json.load(f)

kegg_dict = {}
kegg_mtm_dict  = {}
with open(f'./{species_index}/2_KEGG_filter.tsv','w') as fs:
    fs.write(f"Gene\tPathway\n")
    for gene_key in genes:
        kegg_list = genes[gene_key]["KEGG Info"]
        gene = genes[gene_key]["ID"]
        for kegg in kegg_list:
            kegg_mtm = kegg["Pathway ID"]
            kegg_name = kegg["Pathway Name"]
            if kegg_name == "":
                continue
            
            if kegg_name not in kegg_dict:
                kegg_dict[kegg_name] = [gene]
            else:
                kegg_dict[kegg_name].append(gene)
                
            if kegg_mtm not in kegg_mtm_dict:
                kegg_mtm_dict[kegg_mtm] = kegg_name
    for i in kegg_dict:
        for j in kegg_dict[i]:
            fs.write(f"{j}\t{i}\n")
            
df = pd.read_csv(f"./{species_index}/2_KEGG_filter.tsv",sep="\t")
df.drop_duplicates(subset=["Pathway"],keep='first',inplace=True)
# print(df)
with open(f'./{species_index}/kegg_pathway.txt','w') as f:
    for i in list(df['Pathway']):
        f.write(f"{i}\n")
        
    
with open(f"./{species_index}/pathway_id.tsv",'w') as fs:
    for i in kegg_mtm_dict:
        fs.write(
            f"{i}\t{kegg_mtm_dict[i]}\n"
        )