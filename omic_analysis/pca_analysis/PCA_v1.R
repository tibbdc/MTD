#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
library(plotly)
library(optparse)

## 测试
# Rscript PCA.R -c 'input_file/expression_matrix.csv' -s 'input_file/sample_info.csv' -o 'output_file/pca.png' -j 'output_file/pca.json'


pdf(file = NULL)

option_list <- list(
  make_option(c("-c", "--input_count"), type = "character", default = "", help = "基因表达谱(Count)"), 
  make_option(c("-s", "--input_sample"), type = "character", default = "", help = "样本分组信息"),
  make_option(c("-o", "--output_png"), type = "character", default = "", help = "输出的PCA静态图"),
  # make_option(c("-d", "--output_3d_png"), type = "character", default = "", help = "输出的PCA交互式图的高分辨图片"),
  make_option(c("-j", "--output_json"), type = "character", default = "", help = "输出的PCA交互式图的json文件")  
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查输入文件是否存在
if (!file.exists(opt$input_count)) {
  stop("输入文件不存在: ", opt$input_count)
}
if (!file.exists(opt$input_sample)) {
  stop("输入文件不存在: ", opt$input_sample)
}

# 读取基因表达谱
count_df <- read.csv(file = opt$input_count, header = TRUE, row.names = 1)
# 读取样本分组信息
sample_df <- read.csv(file = opt$input_sample, header = TRUE, row.names = 1)

# 对样本信息进行预处理
rownames(count_df) <- gsub("-", ".", rownames(count_df))
rownames(sample_df) <- gsub("-", ".", rownames(sample_df))
sample_df$Group <- gsub("-", ".", sample_df$Group)

# 创建dds对象用于差异表达分析
sample_df$Group <- factor(sample_df$Group)
deseq2.obj1 <- DESeqDataSetFromMatrix(
    countData = count_df,
    colData = DataFrame(sample_df),
    design = ~Group
)

# 执行DESeq分析
deseq2.obj2 <- DESeq(deseq2.obj1)

# 准备数据进行PCA分析
dds <- deseq2.obj2
dds <- estimateSizeFactors(dds)
rld <- rlog(dds)

# 绘制PCA图
pcaData <- plotPCA(rld, intgroup="Group", returnData=TRUE)
percentVar <- attr(pcaData, "percentVar")

ggplot(pcaData, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=6) +  # 调整点的大小为6，可以根据需要调整
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed() +
  ggtitle("") +
  theme_bw() +  # 使用theme_bw作为基础主题
  theme(
    panel.border = element_rect(colour = "#000000b0", fill=NA, linewidth=0.6)  # 设置边框颜色和线宽
  )

# 保存图片
ggsave(opt$output_png, width=6, height=4, dpi=300)

##  3D PCA
pca_data <- prcomp(t(assay(rld)))
plot_df <- as.data.frame(pca_data$x)

# 计算方差百分比
pca_var <- (pca_data$sdev^2 / sum(pca_data$sdev^2))[1:3] * 100

fig_3d <- plot_ly(
  data = plot_df, x = ~PC1, y = ~PC2, z = ~PC3, text = rownames(plot_df), 
  color = sample_df$Group,
  width = 900, height = 600,
  marker = list(size = 6, line = list(width = 1, color = 'DarkSlateGray'))
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = paste0('PC1: ', round(pca_var[1], 2), '%')),
      yaxis = list(title = paste0('PC2: ', round(pca_var[2], 2), '%')),
      zaxis = list(title = paste0('PC3: ', round(pca_var[3], 2), '%')),
      aspectmode = 'cube',  # 设置坐标轴比例为立方体
      aspectratio = list(x = 1, y = 1, z = 1)  # 设置三个轴的比例都为1
    ),
    template = "simple_white"
  )

# 保存为png，scale设置为4
# kaleido(fig_3d, file = opt$output_3d_png, scale = 4)

# 将Plotly图转换为JSON
json_data <- plotly::plotly_json(fig_3d, jsonedit = FALSE) 

# 将JSON数据写入文件
json_file_path <- opt$output_json
write(json_data, file = json_file_path)