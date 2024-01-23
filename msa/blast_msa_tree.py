
import sys
import os
import pandas as pd
import subprocess



def run_blastp(blast_input_path, blast_output_path, db_prot_name_path, evalue=1e-6):
    """
    根据输入序列信息，对蛋白库进行blastp比对，得到比对结果。
    
    Args:
        blast_input_path (str): 输入txt文件的路径。
        blast_output_path (str): 输出txt文件的路径。
        db_prot_name_path (str): 使用的蛋白数据库路径
        evalue (float): evalue值。
    """
    # 构建命令行参数
    cmd = [
        'blastp',
        '-query', blast_input_path,
        '-out', blast_output_path,
        '-db', db_prot_name_path,
        '-outfmt', '6',
        '-evalue', str(evalue)
    ]

    # 执行命令并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


def run_mafft(workdir, blast_result_path, mafft_result_path):
    """使用mafft进行多序列比对
    
    Args:
        blast_result_path: 输入文件路径
        mafft_result_path: 输出文件路径
    """
    mafft_path = os.path.join(workdir, 'mafft-mac/mafft.bat')
    # print(mafft_path)

    os.system(mafft_path+" --auto "+blast_result_path+" > "+mafft_result_path)
    
# 示例调用
# output_mafft = run_mafft("output_file/blast_seq.fasta", "output_file/mafft_result.fasta")
    

if __name__ == '__main__':

    # 运行blastp
    run_blastp("/Users/dongjiacheng/Desktop/Github/msa/input_file/blast_input.txt",
               "")