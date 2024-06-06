library(optparse)
library(pheatmap)
library(dplyr)
library(ggplot2)

# Define command line options
option_list <- list(
    make_option(c("--input"), type = "character", default = NULL, help = "Input file"),
    make_option(c("--output"), type = "character", default = NULL, help = "Output file"),
    make_option(c("--color_up"), type = "character", default = "default", help = "High expression color"),
    make_option(c("--color_down"), type = "character", default = "default", help = "Low expression color"),
    make_option(c("--color_mid"), type = "character", default = "default", help = "Medium expression color"),
    make_option(c("--show_border"), type = "logical", default = FALSE, help = "Show border"),
    make_option(c("--cluster_rows"), type = "logical", default = TRUE, help = "Cluster rows"),
    make_option(c("--cluster_cols"), type = "logical", default = FALSE, help = "Cluster columns"),
    make_option(c("--scale"), type = "character", default = "log2", help = "Normalization method"),
    make_option(c("--cellwidth"), type = "numeric", default = 20, help = "Cell width"),
    make_option(c("--cellheight"), type = "numeric", default = 20, help = "Cell height"),
    make_option(c("--fontsize"), type = "numeric", default = 10, help = "Font size")
)

# Parse command line options
args <- parse_args(OptionParser(option_list = option_list))

row_zscores <- function(x) {
  row_means <- rowMeans(x)
  row_sds <- apply(x, 1, sd)
  scale(x, center = row_means, scale = row_sds)
}

clean_colnames <- function(colnames) {
  gsub("^X", "", colnames)
}

draw_heatmap <- function(args) {
    data <- read.csv(args$input) # Read CSV file
    data <- as.data.frame(data) # Convert to data frame
    row.names(data) <- make.unique(data[,1]) # Ensure unique row names
    data <- data[,-1]
    colnames(data) <- clean_colnames(colnames(data)) # Clean column names
    
    if(args$color_up == "default" & args$color_down == "default" & args$color_mid == "default") {
      p <- pheatmap(data,
          border_color = ifelse(args$show_border, "black", NA),
          scale = "row",
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
    )
    }
    else {
      my_colors <- colorRampPalette(c(args$color_down, args$color_mid, args$color_up))(n = 256)
      p <- pheatmap(data,
        color = my_colors,
        border_color = ifelse(args$show_border, "black", NA),
        scale = "row",
        cluster_rows = args$cluster_rows,
        cluster_cols = args$cluster_cols,
        cellwidth = args$cellwidth, 
        cellheight = args$cellheight,
        fontsize = args$fontsize,
        legend = TRUE,
        silent = TRUE
    )
    }
    return(p)
}

draw_heatmap_log <- function(args) {
    data <- read.csv(args$input) # Read CSV file
    data <- as.data.frame(data) # Convert to data frame
    row.names(data) <- make.unique(data[,1]) # Ensure unique row names
    data <- data[,-1]
    colnames(data) <- clean_colnames(colnames(data)) # Clean column names
    data <- log2(data+1) # Log2 transformation

    if(args$color_up == "default" & args$color_down == "default" & args$color_mid == "default") {
      p <- pheatmap(data,
          border_color = ifelse(args$show_border, "black", NA),
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
    )
    return(p)
    }
    else {
      my_colors <- colorRampPalette(c(args$color_down, args$color_mid, args$color_up))(n = 256)
      p <- pheatmap(data,
          color = my_colors,
          border_color = ifelse(args$show_border, "black", NA),
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
      )
    return(p)
}
}

if(args$scale == "row")
    p <- draw_heatmap(args)

if(args$scale == "log2")
    p <- draw_heatmap_log(args)

ggsave(args$output, plot = p, width = 10, height = 10, dpi = 300)
