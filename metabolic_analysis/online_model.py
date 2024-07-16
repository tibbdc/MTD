import riptide
import sys
import os
from os.path import join  
import cobra
import json
import pandas as pd
import argparse
from script.ECMpy_function import *
from script.AutoPACMEN_function import *
import re


import matplotlib.pyplot as plt

"""
python3 online_model.py -w /Users/dongjiacheng/Desktop/Github/metabolic_analysis -m GEM -s R2399 -p R2435 -a pfba -o /Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/flux.csv
python3 online_model.py -w /Users/dongjiacheng/Desktop/Github/metabolic_analysis -m ecGEM -s R2399 -p R2435 -a pfba -o /Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/flux.csv
"""



def read_config():
    parser = argparse.ArgumentParser(description='This code is used to model to analyse flux')
    parser.add_argument('--workdir', '-w', required=True, help='workdir')
    parser.add_argument('--model', '-m', required=True, help='input model、substrate、product、product')
    parser.add_argument('--substrate', '-s', required=True, help='input model、substrate、product、product')
    parser.add_argument('--product', '-p', required=True, help='input model、substrate、product、product')
    parser.add_argument('--algorithm', '-a', required=True, help='input model、substrate、product、product')
    parser.add_argument('--output', '-o', required=False, default='output/output.csv ', help='tre file')
    arg = parser.parse_args()
    workdir = arg.workdir
    model = arg.model
    substrate = arg.substrate
    product=arg.product
    algorithm=arg.algorithm
    flux_output = arg.output
    return workdir,model,substrate,product,algorithm,flux_output

def Initial_parameter_selection(workdir, model,substrate,product,algorithm,flux_output):
    #model有GEM,ecGEM
    #substrate输入为反应id例如R2399
    #为产物id

    # 模型背景文件路径
    ecModel_output_file = os.path.join(workdir, 'ecMTM_TurNup_C13.json')
    model_file = os.path.join(workdir, 'model20230921.json')
    
    fluxes_outfile='flux.csv'
    if model=='GEM':#选择模型
        print('选择GEM')
        model=cobra.io.load_json_model(model_file)
        model.reactions.get_by_id("R2399").bounds=(0,0)
        model.reactions.get_by_id(substrate).bounds=(-10,1000)#输入底物
        model.objective=product#输入产物
        if algorithm=='fba':#选择算法
            print('选择fba')
            solution = model.optimize()
            reac=dict(solution.fluxes[abs(solution.fluxes) > 1e-10])
            
            genelist=[]
            reactionlist=[]#反应列表信息
            fluxlist=[]#通量列表信息
            for i in reac:
                reID=str(model.reactions.get_by_id(i))
                reactionlist.append(reID)
                fluxlist.append(reac[i])
                genelist.append(model.reactions.get_by_id(i).gene_reaction_rule)
            df=pd.DataFrame({'reaction':reactionlist,'gene':genelist,'flux':fluxlist})
            df.to_csv(flux_output,index=0)#输出文件
        else:
            print('选择pfba')
            solution=cobra.flux_analysis.pfba(model)
            reac=dict(solution.fluxes[abs(solution.fluxes) > 1e-10])
            print(reac)
            genelist=[]
            reactionlist=[]#反应列表信息
            fluxlist=[]#通量列表信息
            for i in reac:
                reID=str(model.reactions.get_by_id(i))
                fluxlist.append(reac[i])
                reactionlist.append(reID)
                genelist.append(model.reactions.get_by_id(i).gene_reaction_rule)
            df=pd.DataFrame({'reaction':reactionlist,'gene':genelist,'flux':fluxlist})
            df.to_csv(flux_output,index=0)#输出文件

    else:
        print('选择ecGEM')
        model=get_enzyme_constraint_model(ecModel_output_file)
        model.reactions.get_by_id("R2399_reverse").bounds=(0,0)
        model.reactions.get_by_id(substrate).bounds=(-10,1000)#输入底物
        model.objective=product#输入产物
        if algorithm=='fba':#选择算法
            print('选择fba')
            solution = model.optimize()
            reac=dict(solution.fluxes[abs(solution.fluxes) > 1e-10])
            genelist=[]
            reactionlist=[]#反应列表信息
            fluxlist=[]#通量列表信息
            for i in reac:
                reID=str(model.reactions.get_by_id(i))
                fluxlist.append(reac[i])
                reactionlist.append(reID)
                genelist.append(model.reactions.get_by_id(i).gene_reaction_rule)
            print(f'Length of reactionlist: {len(reactionlist)}')
            print(f'Length of genelist: {len(genelist)}')
            print(f'Length of fluxlist: {len(fluxlist)}')
            df=pd.DataFrame({'reaction':reactionlist,'gene':genelist,'flux':fluxlist})
            df.to_csv(flux_output,index=0)#输出文件
        else:
            print('选择pfba')
            solution=cobra.flux_analysis.pfba(model)
            solution = cobra.flux_analysis.pfba(model)
            solution = get_fluxes_detail_in_model(model,solution,flux_output,ecModel_output_file)


def main():

    # read args
   
    workdir,model,substrate,product,algorithm,flux_output=read_config()
    Initial_parameter_selection(workdir,model,substrate,product,algorithm,flux_output)
    


if __name__ == "__main__":
    main()