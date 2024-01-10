
# 在终端中运行此脚本，安装所需的R包： Rscript package_install.R
# 如果没有权限，可以运行 sudo Rscript package_install.R

if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

BiocManager::install(version = "3.18")

BiocManager::install("DESeq2")
BiocManager::install("optparse")
BiocManager::install("clusterProfiler")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("plotly")


# 如果安装失败，可以在R命令行界面中用常规方法安装
# install.packages("XXXXX")