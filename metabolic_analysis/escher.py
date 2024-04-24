import pandas as pd
import subprocess
import os


"""命令行运行示例

# 用户提交一个tsv文件
python GEM_model_transform_1.py --workdir /Users/dongjiacheng/Desktop/Github/metabolic_analysis --input model_input_path/model_input.tsv --output model/model_output.tsv

# 用户提交两个tsv文件
python GEM_model_transform_2.py --workdir /Users/dongjiacheng/Desktop/Github/metabolic_analysis --input_c model_input_path/model_input_control.tsv --input_t model_input_path/model_input_treatment.tsv --output model/model_output_difference.tsv

"""

"""模型运行时间
- 提交单个转录组数据，运行45秒左右
- 提交两个转录组数据，运行80秒左右
"""


def run_mt_model_1(workdir, model_input_path, model_output_path):
    """根据用户提交的转录组数据csv文件, 运行模型, 预测通量，生成结果tsv文件

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

    # # 将输入的csv文件转换为tsv文件
    # df = pd.read_csv(model_input_path, header=0)
    # model_input_dir = os.path.join(workdir, 'input_file')
    # os.makedirs(model_input_dir, exist_ok=True)
    # model_input_path = os.path.join(model_input_dir, 'model_input.tsv')
    # df.to_csv(model_input_path, sep='\t', header=False, index=False)

    # PS:代谢模型运行需要tsv文件,但是用户一般不懂啥是tsv文件，所以需要将用户提交的csv文件转换为tsv文件
    # 获取原始输入csv文件路径，更改扩展名为.tsv并保存
    df = pd.read_csv(model_input_path, header=0)
    model_input_tsv_path = os.path.splitext(model_input_path)[0] + '.tsv'
    df.to_csv(model_input_tsv_path, sep='\t', header=False, index=False)

    # 运行脚本
    cmd = [
        'python', script_path,
        '--workdir', workdir,
        '--input', model_input_tsv_path,
        '--output', model_output_path,
    ]
    # 运行
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        return e.stderr
    
    # 读取模型的输出文件
    df_reaction = pd.read_csv(model_output_path,sep='\t')
    df_reaction.rename(columns={'reactionid': 'ReactionID', 'flux': 'Flux', 'equation': 'Equation'}, inplace=True)
    df_reaction['Flux'] = df_reaction['Flux'].round(6)
    df_reaction['Flux'] = df_reaction['Flux'].apply(lambda x: '{:.6f}'.format(x)) # 将flux列转换为不使用科学记数法的字符串

    df_reaction['Equation'] = df_reaction['Equation'].apply(lambda x: x.split(':', 1)[1])

    data = df_reaction.set_index('ReactionID')['Flux'].to_dict()
    print(data)

    df_reaction['Flux'] = df_reaction['Flux'].astype(float)
    df_reaction = df_reaction.sort_values(by='Flux', ascending=False)
    df_reaction.to_csv(model_output_path,sep='\t',index=False)

    # df_reaction = df_reaction.iloc[:, [0, 2]]

    return data
    


def run_mt_model_2(workdir, model_input_treatment, model_input_control, model_output_path):
    """根据用户提交的两个转录组数据csv文件, 运行模型, 预测通量，生成结果tsv文件，展示两个结果的差异

    Args:
        workdir (str): 工作目录
        model_input_treatment (str): 输入实验组数据
        model_input_control (str): 输入对照组数据
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

    # # 将输入的csv文件转换为tsv文件
    # df1 = pd.read_csv(model_input_treatment, header=0)
    # model_input_dir = os.path.join(workdir, 'input_file')
    # os.makedirs(model_input_dir, exist_ok=True)
    # model_input_treatment = os.path.join(model_input_dir, 'model_input_control.tsv')
    # df1.to_csv(model_input_treatment, sep='\t', header=False, index=False)

    # df2 = pd.read_csv(model_input_control, header=0)
    # model_input_control = os.path.join(model_input_dir, 'model_input_treatment.tsv')
    # df2.to_csv(model_input_control, sep='\t', header=False, index=False)

    # PS:代谢模型运行需要tsv文件,但是用户一般不懂啥是tsv文件，所以需要将用户提交的csv文件转换为tsv文件
    # 获取原始输入csv文件路径，更改扩展名为.tsv并保存
    df1 = pd.read_csv(model_input_treatment, header=0)
    model_input_tsv_treatment = os.path.splitext(model_input_treatment)[0] + '.tsv'
    df1.to_csv(model_input_tsv_treatment, sep='\t', header=False, index=False)

    df2 = pd.read_csv(model_input_control, header=0)
    model_input_tsv_control = os.path.splitext(model_input_control)[0] + '.tsv'
    df2.to_csv(model_input_tsv_control, sep='\t', header=False, index=False)

    # 运行脚本
    cmd = [
        'python', script_path,
        '--workdir', workdir,
        '--input_c', model_input_tsv_treatment,  # 这里input_c 对应的是treatment，不影响结果
        '--input_t', model_input_tsv_control,
        '--output', model_output_path,
    ]
    # 运行
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        return e.stderr
    
    # 读取模型的输出文件
    df_reaction = pd.read_csv(model_output_path,sep='\t')
    df_reaction.rename(columns={'reactionid': 'ReactionID', 
                                'equation_df1': 'Equation_Treatment',
                                'flux_df1': 'Flux_Treatment',
                                'equation_df2': 'Equation_Control',
                                'flux_df2': 'Flux_Control',
                                }, inplace=True)

    
    # 计算两种条件下通量的差异
    df_reaction['Flux_Difference'] = df_reaction['Flux_Treatment'] - df_reaction['Flux_Control']
    df_reaction['Flux_Difference'] = df_reaction['Flux_Difference'].round(6)
    df_reaction['Flux_Difference'] = df_reaction['Flux_Difference'].apply(lambda x: '{:.6f}'.format(x)) # 将flux列转换为不使用科学记数法的字符串

    # 构造字典
    data = df_reaction.set_index('ReactionID')['Flux_Difference'].to_dict()
    print(data)

    df_reaction['Equation'] = df_reaction['Equation_Treatment']
    df_reaction['Equation'] = df_reaction['Equation'].apply(lambda x: x.split(':', 1)[1])

    df_reaction['Flux_Treatment'] = df_reaction['Flux_Treatment'].apply(lambda x: '{:.6f}'.format(x))
    df_reaction['Flux_Control'] = df_reaction['Flux_Control'].apply(lambda x: '{:.6f}'.format(x))

    df_reaction['Flux_Difference'] = df_reaction['Flux_Difference'].astype(float)
    df_reaction['Flux_Treatment'] = df_reaction['Flux_Treatment'].astype(float)
    df_reaction['Flux_Control'] = df_reaction['Flux_Control'].astype(float)
    
    df_reaction = df_reaction.sort_values(by='Flux_Difference', ascending=False)
    df_reaction_save = df_reaction[['ReactionID', 'Equation', 'Flux_Treatment', 'Flux_Control','Flux_Difference']]
    df_reaction_save.to_csv(model_output_path,sep='\t',index=False)

    return data




if __name__ == '__main__':

    import time
    # 记录运算开始时间
    start_time = time.time()
    run_mt_model_1('/Users/dongjiacheng/Desktop/Github/metabolic_analysis', 
                   '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input.csv', 
                   '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/model_output.tsv')
    
    # 记录运算结束时间
    end_time1 = time.time()
    # 计算运算时间
    total_time = end_time1 - start_time
    print(f"模型1运行耗时: {total_time:.2f} 秒")

    run_mt_model_2('/Users/dongjiacheng/Desktop/Github/metabolic_analysis',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input_control.csv',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/input_file/model_input_treatment.csv',
                     '/Users/dongjiacheng/Desktop/Github/metabolic_analysis/output_file/model_output_difference.tsv')
    
    # 记录运算结束时间
    end_time2 = time.time()
    # 计算运算时间
    total_time2 = end_time2 - end_time1
    print(f"模型2运行耗时: {total_time2:.2f} 秒")
