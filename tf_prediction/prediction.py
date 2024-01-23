import subprocess
import os
import pandas as pd


"""命令行运行示例
python tf_running.py -i ./input_file/tf.fasta -o ./output_file -g cpu
python tf_running.py -i ./input_file/tf.fasta -o ./output_file -g cuda:1

- 可以自定义目录

CPU运行: M1芯片
- 运行2个蛋白序列, 用时约2秒
- 运行20个蛋白序列, 用时约16秒

GPU运行不成功
"""

def deepfactor_predict(workdir, input_fasta_path, output_predict_result_path, gpu=False):
    """运行deepfactor, 并返回预测结果
    
    Args:
        workdir: 工作目录
        input_fasta_path: 输入FASTA文件路径
        output_predict_result_path: 输出预测结果路径
        gpu: 是否使用GPU运算, 默认为False
    Returns:
        json格式的预测结果
    """

    # 确保路径正确
    script_path = os.path.join(workdir, 'tf_running.py')
    cmd = [
        'python', script_path,
        '-i', input_fasta_path,
        '-o', output_predict_result_path,
        '-g', 'cuda:1' if gpu else 'cpu',
    ]

    # 检查FASTA文件中的序列数量
    max_sequences = 30
    try:
        with open(input_fasta_path, 'r') as file:
            sequence_count = sum(1 for line in file if line.startswith('>'))
            if sequence_count > max_sequences:
                return f"Error: The input FASTA file contains {sequence_count} sequences, which exceeds the maximum allowed number of {max_sequences}."

        # 运行deepfactor
        print(f"Input FASTA file contains {sequence_count} sequences.")
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except subprocess.CalledProcessError as e:
        return f"An error occurred while running deepfactor: {e.stderr}"
    except Exception as e:
        return f"An unexpected error occurred: {e}"

    output_predict_result = os.path.join(output_predict_result_path, 'prediction_result.txt')
    df_tf = pd.read_csv(output_predict_result, sep='\t', header=None)
    df_tf.columns = ['Sequence_ID', 'Result', 'Score']

    # 转为json格式并返回
    return df_tf.to_json(orient='records')



if __name__ == '__main__':

    import time
    # 记录运算开始时间
    start_time = time.time()

    deepfactor_predict('/Users/dongjiacheng/Desktop/Github/tf_prediction/',
                        '/Users/dongjiacheng/Desktop/Github/tf_prediction/input_file/tf.fasta',
                        '/Users/dongjiacheng/Desktop/Github/tf_prediction/output_file')
    
    # 记录运算结束时间  
    end_time1 = time.time()
    # 计算运算时间
    total_time = end_time1 - start_time
    print(f"运行耗时: {total_time:.2f} 秒")