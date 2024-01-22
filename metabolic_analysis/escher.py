import pandas as pd
import subprocess
import os


"""命令行运行示例

# 用户提交一个tsv文件
python GEM_model_transform_1.py --workdir /Users/dongjiacheng/Desktop/Github/metabolic_analysis --input model_input_path/model_input.tsv --output model/model_output.tsv

# 用户提交两个tsv文件
python GEM_model_transform_2.py --workdir /Users/dongjiacheng/Desktop/Github/metabolic_analysis --input_c model_input_path/model_input_control.tsv --input_t model_input_path/model_input_treatment.tsv --output model/model_output_difference.tsv

"""


def run_mt_model_1(workdir, model_input_path, model_output_path):
    """根据用户提交的转录组数据tsv文件, 运行模型, 预测通量，生成结果tsv文件

    Args:
        workdir (str): 工作目录
        model_input_path (str): 输入文件路径
        model (str): 输出文件路径

    Returns: 返回一个字典, key为反应ID, value为通量的值
    [dict]
        {
            'R00001': 2.0,
            'R00002': 0.3,
            'R00003': 0.03,
            ...
        }
    """
    # 脚本文件路径
    script_path = os.path.join(workdir, 'GEM_model_transform_1.py')

    # 运行脚本
    cmd = [
        'python', script_path,
        '--workdir', workdir,
        '--input', model_input_path,
        '--output', model_output_path,
    ]
    # 运行
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        return e.stderr
    
    # 读取模型的输出文件
    df_reaction = pd.read_csv(model_output_path,sep='\t')
    df_reaction['flux'] = df_reaction['flux'].round(6)
    df_reaction = df_reaction.iloc[:, [0, 2]]

    data = df_reaction.set_index('reactionid')['flux'].to_dict()
    # print(data)

    return data
    

def run_mt_model_2(workdir, model_input_path1, model_input_path2, model_output_path):
    """根据用户提交的两个转录组数据tsv文件, 运行模型, 预测通量，生成结果tsv文件，展示两个结果的差异

    Args:
        workdir (str): 工作目录
        model_input_path1 (str): 输入文件路径1
        model_input_path2 (str): 输入文件路径2
        model_output_path (str): 输出文件路径

    Returns: 返回一个字典, key为反应ID, value为通量差异的值
    [dict]
        {
            'R00001': 2.0,
            'R00002': -1.0,
            'R00003': 0.03,
            ...
        }

    """    
    # 脚本文件路径
    script_path = os.path.join(workdir, 'GEM_model_transform_2.py')

    # 运行脚本
    cmd = [
        'python', script_path,
        '--workdir', workdir,
        '--input_c', model_input_path1,
        '--input_t', model_input_path2,
        '--output', model_output_path,
    ]
    # 运行
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        return e.stderr
    
    # 读取模型的输出文件
    df_reaction = pd.read_csv(model_output_path,sep='\t')

    # 计算两种条件下通量的差异
    df_reaction['flux_difference'] = df_reaction['flux_df2'] - df_reaction['flux_df1']
    df_reaction['flux_difference'] = df_reaction['flux_difference'].round(6)

    # 构造字典
    data = df_reaction.set_index('reactionid')['flux_difference'].to_dict()
    # print(data)

    # 将equation_df1列改名为equation
    df_reaction.rename(columns={'equation_df1': 'equation'}, inplace=True)

    # 保存reactionid、equation、flux_difference三列
    df_reaction = df_reaction[['reactionid', 'equation', 'flux_difference']]
    df_reaction.to_csv(model_output_path,sep='\t',index=False)

    return data



if __name__ == '__main__':

    import time
    # 记录运算开始时间
    start_time = time.time()
    # run_mt_model_1('/Users/dongjiacheng/Desktop/Github/metabolic_analysis', 
    #                '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input.tsv', 
    #                '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/model_output.tsv')
    # 记录运算结束时间
    end_time1 = time.time()
    # 计算运算时间
    total_time = end_time1 - start_time
    print(f"模型1运行耗时: {total_time:.2f} 秒")

    run_mt_model_2('/Users/dongjiacheng/Desktop/Github/metabolic_analysis',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input_control.tsv',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input_treatment.tsv',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/model_output_difference.tsv')
    # 记录运算结束时间
    end_time2 = time.time()
    # 计算运算时间
    total_time2 = end_time2 - end_time1
    print(f"模型2运行耗时: {total_time2:.2f} 秒")