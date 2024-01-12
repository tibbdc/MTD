#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-05-10 06:09:26
LastEditors: wangruoyu wangry@tib.cas.cn
LastEditTime: 2023-09-06 05:23:56
Description: file content
FilePath: /mtd_cdk/lambda/mtd/data/0/generate_gene_id.py
'''
import os
import sys
import json

species_index = sys.argv[1]

with open(f"./{species_index}/gene_information.json",'r') as f:
    genes = json.load(f)
    
with open(f'./{species_index}/gene_id.tsv','w') as fs:
    fs.write(f"ID\tGene\n")
    for gene in genes:
        fs.write(f"{gene}\t{genes[gene]['ID']}\n")