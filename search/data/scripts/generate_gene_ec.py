import os
import sys
import json
import pandas as pd

species_index = sys.argv[1]

with open(f"./{species_index}/gene_information.json",'r') as f:
    genes = json.load(f)
    
    
    

with open(f"./{species_index}/Gene_Descripition_Symbol_EC.tsv",'w') as fs:
    fs.write(
        "Gene\tDescription\tSymbol\tEC\n"
    )
    
    for gene_key in genes:
        gene = genes[gene_key]["ID"]
        des = genes[gene_key]["Description"]
        symbol = genes[gene_key]["Name"]
        ec = genes[gene_key]["EC"]

        fs.write(
        f"{gene}\t{des}\t{symbol}\t{ec}\n"
    )
    