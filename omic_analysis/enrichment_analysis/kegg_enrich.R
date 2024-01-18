#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(dplyr)


## 测试
#  Rscript kegg_enrich.R -w /Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis -s "Myceliophthora thermophila" -p 0.05 -i "input_file/gene_list.txt" -o "output_file/kegg.tsv"


# 设置命令行选项
option_list <- list(
  make_option(c("-w", "--workdir"), type = "character", default = NULL, help = "工作目录路径", metavar = "workdir"),
  make_option(c("-s", "--species"), type = "character", default = "", help = "菌种名称", metavar = "species"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, help = "pvalue阈值", metavar = "pvalue"),
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件路径", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "输出文件路径", metavar = "file")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 设置工作目录
if (!is.null(opt$workdir)) {
  setwd(opt$workdir)
}

background_dir <- "background_file"

# 根据不同的物种读取不同的文件
if (opt$species == "Aspergillus niger") {
  ko2name_path <- file.path(background_dir, "kegg_name_an.txt")
  ko2gene_path <- file.path(background_dir, "kegg_gene_an.txt")
} else if (opt$species == "Myceliophthora thermophila") {
  ko2name_path <- file.path(background_dir, "kegg_name_mt.txt")
  ko2gene_path <- file.path(background_dir, "kegg_gene_mt.txt")
} else if (opt$species == "Trichoderma reesei") {
  ko2name_path <- file.path(background_dir, "kegg_name_tr.txt")
  ko2gene_path <- file.path(background_dir, "kegg_gene_tr.txt")
} else if (opt$species == "Neurospora crassa") {
  ko2name_path <- file.path(background_dir, "kegg_name_nc.txt")
  ko2gene_path <- file.path(background_dir, "kegg_gene_nc.txt")
}
  ko2gene <- read.delim(ko2gene_path, stringsAsFactors = FALSE)
  ko2name <- read.delim(ko2name_path, stringsAsFactors = FALSE)


# 读取gene list
genelist <- read.csv(file = opt$input, row.names = 1)
genelist <- as.character(rownames(genelist))

# 如果genelist中有重复的基因，只保留一个
genelist <- unique(genelist)

# kegg
enrich_kegg <- enricher(genelist,
  TERM2GENE = ko2gene,
  TERM2NAME = ko2name,
  pAdjustMethod = "BH", # 使用FDR进行校正
  pvalueCutoff = opt$pvalue,
  qvalueCutoff = 1
)

# 保存结果到指定的输出文件
write.table(enrich_kegg, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
