#-*- coding:utf-8 -*-

import subprocess
import os

def run_pca(workdir, input_count_path, input_sample_path, output_hmtl_2d_path, output_hmtl_3d_path):
    """
    根据输入的表达量矩阵和样本信息表，运行R脚本，生成PCA分析的结果

    Args:
        workdir: 工作目录
        input_count_path: 输入的表达量矩阵路径
        input_sample_path: 输入的样本信息表路径
        output_hmtl_2d_path: 输出的2D PCA分析结果路径
        output_hmtl_3d_path: 输出的3D PCA分析结果路径
    """

    # R脚本的路径，需要师哥你改路径
    script_path = os.path.join(workdir, 'PCA.R')

    cmd = [
        'Rscript', script_path,
        '--input_count', input_count_path,
        '--input_sample', input_sample_path,
        '--output_html_2d', output_hmtl_2d_path,
        '--output_html_3d', output_hmtl_3d_path,
    ]

    # 执行R脚本并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


if __name__ == "__main__":

    # 示例调用
    run_pca(
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/input_file/expression_matrix.csv',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/input_file/sample_info.csv',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/output_file/pca_2d.html',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/output_file/pca_3d.html'
    )