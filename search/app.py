#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-03-29 08:45:30
LastEditors: wangruoyu wangry@tib.cas.cn
LastEditTime: 2023-06-28 05:30:41
Description: file content
FilePath: /mtd_cdk/lambda/mtd/app.py
'''
# coding:utf-8

import os
import sys
import json
import uuid
import boto3



try:
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import scipy
except Exception as e:
    s3 = boto3.resource('s3')
    object_name = "py38-seaborn-scipy-kegg-viewer.tar.gz"
    local_package = os.path.join("/tmp/", object_name)
    s3.Object(os.environ["s3Reference"], f"lambda-layers/{object_name}").download_file(local_package)
    cmd = "tar zxf "+local_package+" -C /tmp/"
    os.system(cmd)
    sys.path.insert(0,"/tmp/python/")
    print(sys.path)
    cmd = "rm " + local_package
    os.system(cmd)

    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import scipy

from config import  Strain

args_dict = {
            'ACL': "public-read",
            'ContentType':"image/svg+xml"
            }


def api_response(statuscode,data):
    # lambda response to api gateway proxy lambda
    response_body = {
        "statusCode":statuscode,
        "data" : data
    }
    response = {
        "statusCode": 200,
        "headers":{
            'Access-Control-Allow-Headers': 'Content-Type,X-Amz-Date,Authorization,X-Api-Key',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': 'OPTIONS,POST,GET'
        },
        "body": json.dumps(response_body)
    }
    return response

def api_url_check(path):
    paths = ['search','expression','plot']
    avaliable_paths = [f"/api/{i}" for i in paths]
    print(f'url path: {path}')
    if path not in avaliable_paths:
        raise ValueError(f"url path {path} not support")

def upload_svg(plot_file:str):
    # upload svg
    key = f"plot/{plot_file.split('/')[-1]}" 
    s3.meta.client.upload_file(plot_file, os.environ["s3Result"], key,ExtraArgs=args_dict)

    url = f"https://{os.environ['s3Result']}.s3.amazonaws.com/{key}"
    return url
    
def lambda_handler(event,context):
    print(event)
    try:
        path = event["path"]
        api_url_check(path)
        event = json.loads(event['body'])
        print(event)
        if path == "/api/search":
            name = event["strain_name"]
            database = event["database"]
            value = event["value"]
            s = Strain(name)
            print('Strain instance success.')
            data = s.search(database,value)
            # if database == "pathway":
            #     data["search"] = upload_svg(data["search"])
        elif path == "/api/expression":
            name = event["strain_name"]
            values = event["gsm_list"]
            gene_values = event["gene_list"]
            s = Strain(name)
            data = s.s_expression(values,gene_values)
        elif path == "/api/plot":
            name = event["strain_name"]
            gsm_values = event["gsm_list"]
            gene_values = event["gene_list"]
            s = Strain(name)
            plot_file = f"/tmp/{uuid.uuid4()}.svg"
            s.plot(gsm_values,gene_values,plot_file)
            
            url = upload_svg(plot_file)
            data = {
                "url":url
            }
        return api_response(200,data)
    except Exception as e:
        print(str(e))

        code = 500
        data = {
            "code":code,
            "msg":str(e)
        } 
        return api_response(code,data)



if __name__ == "__main__":
    
    
    
    
    event = {
        "path":"/api/expression",
        "body":"{\"strain_name\":\"Myceliophthora thermophila\",\"gsm_list\":[\"GSM4074527\",\"GSM4074528\"]}"
    }
    
    event = {
        "path":"/api/plot",
        "body":"{\"strain_name\":\"Myceliophthora thermophila\",\"gsm_list\":[\"GSM4074527\",\"GSM4074528\",\"GSM4074531\",\"GSM4074532\",\"GSM4074533\"]}"
    }

    event = {
        "path":"/api/search",
        "body":"{\"strain_name\":\"Myceliophthora thermophila\",\"database\":\"go\",\"value\":\"GO:0004865\"}"
    }

    event = {
        "path":"/api/search",
        "body":"{\"strain_name\":\"Myceliophthora thermophila\",\"database\":\"pathway\",\"value\":\"Glycolysis / Gluconeogenesis\"}"
    }
    
    event = {
        "path":"/api/search",
        "body":"{\"strain_name\":\"Myceliophthora thermophila\",\"database\":\"gene\",\"value\":\"MYCTH_100011\"}"
    }
    print(lambda_handler(event,{}))