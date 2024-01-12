#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-03-29 08:09:56
LastEditors: wangruoyu wangry@tib.cas.cn
LastEditTime: 2023-09-06 06:14:47
Description: file content
FilePath: /mtd_cdk/lambda/mtd/config.py
'''
import os
import sys
import json
import subprocess
import uuid

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy


DATA = {
    "Myceliophthora thermophila":"0",
    "Aspergillus niger":"1",
    "Neurospora crassa":"2",
    "Trichoderma reesei":"3",
}


class Strain:

    def __init__(self,name:str) -> None:
        self.name = name
        self.index = self.check_strain()
        self.pathway = f"data/{self.index}/2_KEGG_filter.tsv"
        self.go = f"data/{self.index}/1_Go.tsv"
        self.gene = f"data/{self.index}/gene_information.json"
        self.gene_info = self.read_gene_information()
        self.sample = f"data/{self.index}/sample.tsv"
        self.expression = f"data/{self.index}/exp.tsv"
        self.kegg_id = f"data/{self.index}/pathway_id.tsv"
        self.sample_condition = f"data/{self.index}/sample_condition.tsv"
        self.gene_ec_info = f"data/{self.index}/Gene_Descripition_Symbol_EC.tsv"
    
    def read_gene_information(self):
        with open(self.gene,'r') as f:
            genes = json.load(f)
            return genes
    def check_strain(self):
        """
        通过strain名返回相关文件
        :param strain: _description_
        :type strain: str
        """
        if not DATA.get(self.name):
            raise ValueError(f"the strain {self.name} not support.")
        return DATA.get(self.name)
    
    def get_samples(self) -> dict:
        """_summary_

        Returns:
            dict: _description_
        """
        df = pd.read_csv(self.sample,sep="\t")
        return df
    def get_kegg_pathway(self,value:str = None) -> pd.DataFrame:
        """_summary_

        Returns:
            pd.DataFrame: _description_
        """
        df = pd.read_csv(self.pathway,sep="\t")
        if value:
            print("go search value: ", value)
            df = df[df["Pathway"].str.strip() == value.strip()]
        print('Pathway df Gene list: ', list(df["Gene"]))
        return df
    def get_go(self,value:str =None) -> pd.DataFrame:
        """_summary_

        Returns:
            pd.DataFrame: _description_
        """
        df = pd.read_csv(self.go,sep="\t")
        if value:

            df = df[df["GOID"].str.strip() == value.strip()]
        print('GO df Gene list : ',list(df["Gene"]))
        return df
    def get_expression(self,gsm_list:list) -> pd.DataFrame:
        """_summary_

        Returns:
            dict: _description_
        """
        headers= ["Gene id"] + gsm_list
        df = pd.read_csv(self.expression,sep="\t")
        df = df[headers]
        return df

    def get_expression_by_gene(self,gene_list:list) -> pd.DataFrame:
        """
        """    
        gene_list = list(set(gene_list))
        print('get expression by gene input gene list :',gene_list)
        df = pd.read_csv(self.expression,sep="\t")
        # print(df)
        df = df[df["Gene id"].isin(gene_list)]
        print('get expression by gene: ', list(df["Gene id"]))
        return df
    
    def get_gene_expression_condition(self,gsm_list:list):
        df = pd.read_csv(self.sample_condition,sep="\t")[["Sample ID","Study Title","Condition","Series Accession","Related Sample"]]
        title_list = []
        condition_list = []
        accession_list = []
        for i  in gsm_list:
            tmp_df = df[df["Sample ID"]==i]
            if not tmp_df.empty:
                # title
                value = tmp_df["Study Title"].values[0]
                title_list.append(value)
                # condition
                print(tmp_df)
                c_value = f"{tmp_df['Condition'].values[0]} {tmp_df['Sample ID'].values[0]} {tmp_df['Related Sample'].values[0]}"
                condition_list.append(c_value)
                
                accession_value = f"{tmp_df['Series Accession'].values[0]}"
                accession_list.append(accession_value)
                
        return title_list,condition_list,accession_list

        
    def plot(self,gsm_list:list,gene_list:list,plot_file:str) -> str:
        """_summary_

        Args:
            gsm_list (list): _description_
            gene_list (list): _description_
            output_file (str): _description_

        Returns:
            str: heatmap png filepath
        """
        df = self.get_expression(gsm_list)
        df = df[df["Gene id"].isin(gene_list)]
        df.set_index('Gene id',inplace=True)
        print(df)
        # plot
        index = True
        if len(gene_list) >= 15 or len(gsm_list) > 15:
            index = False
        cluster_fig = sns.clustermap(df,annot=index,cmap="RdBu_r")
        # cluster_fig.set_xlabel('Sample id')
        cluster_fig.savefig(plot_file,dpi=400,bbox_inches="tight")
        plt.close()
        return plot_file

    def remove_duplicate(self,input_df:pd.DataFrame,stype:str) -> pd.DataFrame:
        """

        :param input_df: _description_
        :type input_df: pd.DataFrame
        :param stype: _description_
        :type stype: str
        :return: _description_
        :rtype: pd.DataFrame
        """
        dup_items = {}
        # gene重复，第二个往后的
        dup_df = input_df[input_df.duplicated(subset="Gene")].reset_index(drop=True)
        
        no_dup_df = input_df.drop_duplicates(subset="Gene").reset_index(drop=True)
        
        if stype == "go":
            no_dup_df["gene_ontology_1"] = no_dup_df['GOID'] + "," + no_dup_df['ONTOLOGY'] + "," + no_dup_df['TERM'] +";"
            
            print('no dup df : ',no_dup_df)

            for i in range(len(dup_df)):
                gene = dup_df.iloc[i]["Gene"]
                goid = dup_df.iloc[i]['GOID']
                ontology = dup_df.iloc[i]['ONTOLOGY']
                term = dup_df.iloc[i]['TERM']
                
                v = f"{goid},{ontology},{term};"
                if gene not in dup_items:
                    dup_items[gene] = {
                        "gene":gene,
                        "gene_ontology_2":v
                    }
                else:
                    dup_items[gene]["gene_ontology_2"] += v
            
            new_dup_df = pd.DataFrame(dup_items.values())
            print('new duplicate df: ',new_dup_df)
            
            if not new_dup_df.empty:
                df = pd.merge(
                    no_dup_df,
                    new_dup_df,
                    left_on="Gene",
                    right_on = "gene",
                    how="left"
                )
                df[["gene_ontology_2"]] = df[["gene_ontology_2"]].fillna("")
                df["GENE_ONTOLOGY"] = df["gene_ontology_1"] + df["gene_ontology_2"]
            else:
                df = no_dup_df
                df["GENE_ONTOLOGY"] = df["gene_ontology_1"]
            print('new merge df: ',df)
            return df[["Gene","GENE_ONTOLOGY"]]
        elif stype == "pathway":
            
            print('no dup df : ',no_dup_df)

            for i in range(len(dup_df)):
                gene = dup_df.iloc[i]["Gene"]
                pathway = f"{dup_df.iloc[i]['Pathway']};"
                
                if gene not in dup_items:
                    dup_items[gene] = {
                        "gene":gene,
                        "pathway":pathway
                    }
                else:
                    dup_items[gene]["pathway"] += pathway
            
            new_dup_df = pd.DataFrame(dup_items.values())
            print('new duplicate df: ',new_dup_df)
            
            if not new_dup_df.empty:
                df = pd.merge(
                    no_dup_df,
                    new_dup_df,
                    left_on="Gene",
                    right_on = "gene",
                    how="left"
                )
                df[["pathway"]] = df[["pathway"]].fillna("")
                df["Pathway"] = df["Pathway"] + ";" + df["pathway"]
            else:
                df = no_dup_df
            print('new merge df: ',df)
            return df[["Gene","Pathway"]]
            
    def search(self,database:str,value:str) -> dict:
        """_summary_

        Args:
            database (str): _description_
            value (str): _description_

        Returns:
            dict: _description_
        """
        d_value =  getattr(self,database)
        
        if database == "gene":
            if value.startswith('MYCTH'):
                pass
            elif value.startswith('ANI'): 
                pass
            elif value.startswith('NCU'): 
                pass
            elif value.startswith('TRIREDRAFT'):
                pass
            else:
                raise ValueError('Search Gene please use correct ID.')
                
            gene_info_json = self.gene_info[f"KEY:{value}"]
            gene_value_df = self.get_expression_by_gene([value])
            gene_value_list = gene_value_df.values.tolist()[0]
            print('gene value list',gene_value_list)
            samples = gene_value_df.columns[1:]
            name_list,condition_list,pmid_list = self.get_gene_expression_condition(samples)
            data = {
                "search":gene_info_json,
                "expression":{
                    "pmid":pmid_list,
                    "full_name":name_list,
                    "condition":condition_list,
                    "value":gene_value_list[1:]
                }
            }
        else:   
            if database == "pathway":
                p_df = self.get_kegg_pathway(value)
                go_df = self.get_go()
                direction = 'left'
            elif database == "go":
                p_df = self.get_kegg_pathway()
                go_df = self.get_go(value)
                direction = 'right'
            go_df = self.remove_duplicate(go_df,'go')
            p_df = self.remove_duplicate(p_df,'pathway')
            print(f"search {database} database, value is {value}")
            df = pd.merge(
                p_df,
                go_df,
                left_on="Gene",
                right_on = "Gene",
                how=direction
            )
            # df = df.dropna(axis=0,how='any')
            print('Pathway GO merge df Gene list: ', list(df["Gene"]))
            # search_json_string = df.to_json(orient = 'records')
            expression_json_string = self.get_expression_by_gene(list(df["Gene"])).to_json(orient = 'records')
        
            # pathway 返回 svg
            data = {
                "search":{
                    "url": f"https://{os.environ['s3Result']}.s3.amazonaws.com/svg/{self.index}/{self.get_pathway_id_by_pathway(value)}.svg",
                    "genes":json.loads(self.get_pathway_contains_genes(value).to_json(orient = 'records'))
                    },
                "expression":json.loads(expression_json_string)
            }
        return data

    def s_expression(self,gsm_list:list,gene_list:list) -> dict:
        df = self.get_expression(gsm_list)
        df = df[df["Gene id"].isin(gene_list)]
        json_string = df.to_json(orient = 'records')
        return json.loads(json_string)
    
    def get_pathway_id_by_pathway(self,pathway_name: str) -> str:
        """根据pathway的名字，依据kegg_id表，返回在相应pathway的id，生成pathway的通路图，保存在kegg-viewers目录下

        Args:
        * workdir: 工作目录
        * species_name: 菌种的名字
        * pathway_name: pathway的名字
        ----------
        Returns:
            pathway_id: pathway的id

        Usage: 
        安装：pip install kegg_viewer
        :   生成pathway的通路图的命令为：kegg_viewer -p pathway_id -O output -c cache
        :   pathway_id为pathway在kegg数据库中的id，output_dir为输出的通路图的目录，cache为缓存的目录
        :   通路图的svg矢量图输出目录为：workdir/species_name/kegg-viewers
        :   缓存的目录为：workdir/species_name/cache
        """
 
        # 读取kegg_id表
        df_kegg_id = pd.read_csv(self.kegg_id, sep='\t', names=['Id', 'Pathway'] )

        df = df_kegg_id[df_kegg_id['Pathway'].str.contains(pathway_name)]
        if df.empty:
            raise ValueError(f"The Pathway {pathway_name} can not be found in pathway id file")
        pathway_id = df['Id'].tolist()[0]
        print("pathway_id",pathway_id)
        return pathway_id
    
    def get_pathway_contains_genes(self,pathway_name):
        """_summary_

        Args:
            pathway_name (_type_): _description_
        """
        p_df = self.get_kegg_pathway(pathway_name)
        gene_list = list(p_df["Gene"])
        
        df = pd.read_csv(self.gene_ec_info,sep="\t")
        df = df[df["Gene"].isin(gene_list)]
        print('get gene ec info by gene: ', list(df["Gene"]))
        return df