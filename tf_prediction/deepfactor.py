import subprocess
import os
import pandas as pd
import time


"""命令行运行示例
python tf_running.py -i ./Dataset/tf.fasta -o ./result -g cpu
python tf_running.py -i ./Dataset/tf.fasta -o ./result -g cuda:1

- 可以自定义目录

CPU运行: M1芯片
- 运行2个蛋白序列, 用时约2秒
- 运行20个蛋白序列, 用时约16秒

GPU运行不成功
"""


def deepfactor_predict(workdir, input_fasta_path, output_predict_result_path, gpu=False):
    """运行deepfactor, 并返回预测结果

    Args:
        workdir (str): deepfactor工作目录
        input_fasta_path (str): 输入fasta文件路径
        output_predict_result_path (str): 输出预测结果文件目录
        gpu (bool, optional): 是否使用gpu. Defaults to False.

    Returns:
        str: 预测结果json格式
    """

    script_path = os.path.join(workdir, 'tf_running.py')
    start_time = time.time()  # 记录开始时间

    if gpu:
        cmd = [
            'python', script_path,
            '-i', input_fasta_path,
            '-o', output_predict_result_path,
            '-g', 'cuda:1',
        ]
    else: 
        cmd = [
            'python', script_path,
            '-i', input_fasta_path,
            '-o', output_predict_result_path,
            '-g', 'cpu',
        ]

    # 运行deepfactor
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        return e.stderr  
    
    # 读取结果文件, 生成dataframe, 并重命名列名
    output_predict_result = os.path.join(output_predict_result_path, 'prediction_result.txt')
    df_tf = pd.read_csv(output_predict_result, sep='\t', header=None)
    df_tf.columns = ['Sequence_ID', 'Result', 'Score']

    # 转为json格式
    df_tf_json = df_tf.to_json(orient='records')

    end_time = time.time()  # 记录结束时间
    total_time = end_time - start_time  # 计算总运行时间
    print(f"运行耗时: {total_time:.2f} 秒")  # 输出运行时间

    # print(df_tf_json)
    return df_tf_json


if __name__ == '__main__':

    deepfactor_predict('/Users/dongjiacheng/Desktop/Github/tf_prediction/',
                        '/Users/dongjiacheng/Desktop/Github/tf_prediction/Dataset/tf.fasta',
                        '/Users/dongjiacheng/Desktop/Github/tf_prediction/result')