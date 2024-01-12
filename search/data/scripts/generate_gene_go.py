'''
Author: wangruoyu wangry@tib.cas.cn
Date: 2023-06-28 01:37:54
LastEditors: wangruoyu wangry@tib.cas.cn
LastEditTime: 2023-09-06 07:16:51
FilePath: /mtd_cdk/lambda/mtd/data/0/generate_gene_go.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''
import os
import sys
import json
import pandas as pd

species_index = sys.argv[1]
with open(f"./{species_index}/gene_information.json",'r') as f:
    genes = json.load(f)

kegg_dict = {}
with open(f'./{species_index}/1_Go.tsv','w') as fs:
    fs.write(f"Gene\tGOID\tONTOLOGY\tTERM\n")
    for gene_key in genes:
        kegg_list = genes[gene_key]["GO Info"]
        gene = genes[gene_key]["ID"]
        for kegg in kegg_list:
            go_id = kegg["ID"]
            go_term = kegg["ONTOLOGY"]
            go_ONTOLOGY = kegg["Type"]
            if go_id == "":
                continue
            fs.write(f"{gene}\t{go_id}\t{go_ONTOLOGY}\t{go_term}\n")
            
df = pd.read_csv(f"./{species_index}/1_Go.tsv",sep="\t")
df.drop_duplicates(subset=["GOID"],keep='first',inplace=True)
# print(df)
with open(f'./{species_index}/go.txt','w') as f:
    for i in list(df['GOID']):
        f.write(f"{i}\n")