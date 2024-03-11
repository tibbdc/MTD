


### Notion
https://dongjc.notion.site/dongjc/MTD-e39ef7b0e627460eab37b35ee6700590


### 2024.1.23
- 更新了omic_analysis/heatmap模块的代码，修复了一个错误
- 添加了metabolic_analysis与tf_prediction的内容：代码、文档、界面


### 2024.1.24
- 添加了msa（蛋白序列比对、多序列比对、进化树构建）的内容：代码、文档、界面


### 2024.1.25
- 优化了msa中的代码


### 2024.1.29
- 更新了metabolic_analysis/escher.py的代码
    1. 代码生成的文件中存在科学记数法的数字，不利于前端展示，已修改，并更新了生成的表格
    2. 添加了Python调用Escher的代码，供参考（实际环境后端只需要返回具体的值，由前端展示）


### 2024.2.4
- 更新了omic_analysis/pca_analysis/pca.py与PCA.R的代码
    1. 添加了2D PCA图的显示
    2. 删除了静态PCA图的生成


### 2024.2.27
- 更新了msa/blast_msa_tree.py的代码，在linux环境进行了测试
    1. 运行Blastp
    2. 运行mafft与fasttree，即多序列比对与蛋白进化树构建


### 2024.3.11
- 更新了tf_prediction/prediction.py的代码，删除了部分无意义内容