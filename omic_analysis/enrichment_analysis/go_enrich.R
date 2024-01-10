#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(dplyr)


## 测试
#  Rscript go_enrich.R -s 'Myceliophthora thermophila' -p 0.05 -i 'input_file/gene_list.txt' -o 'output_file/enrich_kegg.tsv'

# 设置命令行选项
option_list <- list(
  make_option(c("-s", "--species"), type = "character", default = "", help = "菌种名称", metavar = "species"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, help = "pvalue阈值", metavar = "pvalue"),
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件路径", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "输出文件路径", metavar = "file")
)

# 解析命令行参数
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# 检查输入文件是否存在
if (!file.exists(opt$input)) {
  stop("输入文件不存在: ", opt$input)
}

# 背景文件，需要师哥你修改一下路径
go2name <- read.delim('/Users/dongjiacheng/Desktop/coder/mtd/code/omic_analysis/enrichment_analysis/background_file/go2name.txt', stringsAsFactors=FALSE)
if (opt$species == "Aspergillus niger") {
  go2gene <- read.delim('/Users/dongjiacheng/Desktop/coder/mtd/code/omic_analysis/enrichment_analysis/background_file/go_gene_an.txt', stringsAsFactors=FALSE)
} else if (opt$species == "Myceliophthora thermophila") {
  go2gene <- read.delim('/Users/dongjiacheng/Desktop/coder/mtd/code/omic_analysis/enrichment_analysis/background_file/go_gene_mt.txt', stringsAsFactors=FALSE)
} else if (opt$species == "Trichoderma reesei") {
  go2gene <- read.delim('/Users/dongjiacheng/Desktop/coder/mtd/code/omic_analysis/enrichment_analysis/background_file/go_gene_tr.txt', stringsAsFactors=FALSE)
} else if (opt$species == "Neurospora crassa") {
  go2gene <- read.delim('/Users/dongjiacheng/Desktop/coder/mtd/code/omic_analysis/enrichment_analysis/background_file/go_gene_nc.txt', stringsAsFactors=FALSE)
}

# 读取gene list
genelist <- read.csv(file = opt$input, row.names = 1)
genelist <- as.character(rownames(genelist))

# 如果genelist中有重复的基因，只保留一个
genelist <- unique(genelist)

# GO
go2gene = split(go2gene , with(go2gene , CLASS))

enrich_MF = enricher(genelist,
                    TERM2GENE=go2gene [['MF']][c(1,2)],
                    TERM2NAME=go2name,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

enrich_BP = enricher(genelist,
                    TERM2GENE=go2gene [['BP']][c(1,2)],
                    TERM2NAME=go2name,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

enrich_CC = enricher(genelist,
                    TERM2GENE=go2gene [['CC']][c(1,2)],
                    TERM2NAME=go2name,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

# 合并三个表
enrich_MF = as.data.frame(enrich_MF)
enrich_BP = as.data.frame(enrich_BP)
enrich_CC = as.data.frame(enrich_CC)
enrich_all <- bind_rows(
  mutate(enrich_MF, category = "MF"),
  mutate(enrich_BP, category = "BP"),
  mutate(enrich_CC, category = "CC")
)

# 保存结果到指定的输出文件
write.table(enrich_all, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)