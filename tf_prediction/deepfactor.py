import subprocess
import os
import pandas as pd


workdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'deepfactor')

"""
python tf_running.py -i ./Dataset/tf.fasta -o ./result -g cpu
python tf_running.py -i ./Dataset/tf.fasta -o ./result -g cuda:1
"""
def deepfactor_predict(workdir):

    script_path = os.path.join(workdir, 'tf_running.py')
    input_fasta_path = os.path.join(workdir, 'Dataset', 'tf.fasta')
    output_predict_result_path = os.path.join(workdir, 'result')

    # python tf_running.py -i ./Dataset/tf.fasta -o ./result -g cpu
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
        return e.stderr  # 正确使用 stderr 而不是 stderrs
    
    # 读取结果文件, 生成dataframe, 并重命名列名
    df_tf = pd.read_csv(output_predict_result_path, sep='\t', header=None)
    df_tf.columns = ['Sequence_ID', 'Result', 'Score']

    # 转为json格式
    df_tf_json = df_tf.to_json(orient='records')
    return df_tf_json


if __name__ == '__main__':

    deepfactor_predict('/Users/dongjiacheng/Desktop/Github/tf_prediction/')