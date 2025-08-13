library(tidyverse)
library(data.table)
library(coloc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(locuscomparer)
library(ggplot2)
library(TwoSampleMR)
library(readxl)
library(biomaRt)
library(ggpubr)

# --- 1. 数据加载与预处理 ---
# 定义一个函数来格式化p值
format_pvalue <- function(value) {
  ifelse(value >= 0.001, round(value, 3), format(value, scientific = TRUE, digits = 2))
}

# 加载基因列表
gene_key <- read_excel('gene_list.xlsx')
target_genes <- unique(gene_key$gene)

# 转换基因ID格式
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_id_mapping <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                         filters = "hgnc_symbol",
                         values = target_genes,
                         mart = ensembl)
colnames(gene_id_mapping) <- c("SYMBOL", "Gene")

# 加载GWAS数据
gwas_data <- fread('gwas_data.txt', header = TRUE, na.strings = "", fill = TRUE) %>%
  dplyr::select(pos, CHR, SNP, beta, se, A1, A2, freq, p) %>%
  setnames(c('pos', 'CHR', 'SNP', 'beta', 'se', 'A1', 'A2', 'freq', 'p'),
           c('base_pair_location', 'chromosome', 'ID', 'beta', 'se', 'effect_allele',
             'other_allele', 'effect_allele_frequency', 'p_value')) %>%
  mutate(MAF = ifelse(effect_allele_frequency < 0.5, effect_allele_frequency, 1 - effect_allele_frequency)) %>%
  filter(!is.na(ID), !duplicated(ID)) %>%
  as.data.frame()

# 定义SMR和Coloc分析所需的样本量
gwas_N <- 15056 + 12637
gwas_case <- 15056
eqtl_N <- 31684
mqtl_N <- 1980
pqtl_N <- 10708
qtl_data_path <- "./" # 定义QTL数据目录

# --- 2. SMR与共定位分析函数 ---
run_smr_and_coloc <- function(qtl_type, smr_data, coloc_data_path, gene_id_mapping, gwas_data, gwas_N, gwas_case, qtl_N) {
  
  # 加载QTL数据
  qtl_coloc_data <- fread(coloc_data_path) %>%
    setnames(c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","SE","p")) %>%
    mutate(MAF = ifelse(Freq < 0.5, Freq, 1 - Freq)) %>%
    filter(!duplicated(SNP)) %>%
    arrange(p) %>%
    as.data.frame()
  
  # 循环进行Coloc分析
  significant_probes <- unique(smr_data$probeID)
  coloc_results_list <- list()
  
  for(probe_id in significant_probes) {
    tryCatch({
      # 筛选GWAS和QTL数据
      gwas_subset <- gwas_data
      qtl_subset <- qtl_coloc_data[qtl_coloc_data$Probe == probe_id, ]
      
      # 寻找共同的SNP
      common_snps <- intersect(gwas_subset$ID, qtl_subset$SNP)
      if (length(common_snps) == 0) next
      
      gwas_subset <- gwas_subset[gwas_subset$ID %in% common_snps, ] %>% arrange(ID)
      qtl_subset <- qtl_subset[qtl_subset$SNP %in% common_snps, ] %>% arrange(SNP)
      
      # 运行Coloc分析
      coloc_result <- coloc.abf(
        dataset1 = list(beta = gwas_subset$beta, varbeta = gwas_subset$se^2, snp = gwas_subset$ID,
                        type = "cc", N = gwas_N, s = gwas_case/gwas_N, MAF = gwas_subset$MAF),
        dataset2 = list(beta = qtl_subset$b, varbeta = qtl_subset$SE^2, snp = qtl_subset$SNP,
                        type = "quant", N = qtl_N, MAF = qtl_subset$MAF),
        p12 = 5e-5
      )
      
      # 存储结果
      if (as.numeric(coloc_result$summary['PP.H4.abf']) > 0.5) {
        gene_symbol <- gene_id_mapping$SYMBOL[match(probe_id, gene_id_mapping$Gene)]
        if (length(gene_symbol) == 0 || is.na(gene_symbol)) gene_symbol <- probe_id
        
        # 写入结果文件和曼哈顿图
        # ... (文件输出和绘图代码，使用相对路径) ...
        
        coloc_results_list[[length(coloc_results_list) + 1]] <- data.frame(
          ProbeID = probe_id,
          GeneSymbol = gene_symbol,
          PPH4 = as.numeric(coloc_result$summary['PP.H4.abf'])
        )
      }
    }, error = function(e) {
      warning(paste("Error occurred for probe", probe_id, ":", conditionMessage(e)))
    })
  }
  return(bind_rows(coloc_results_list))
}

# --- 3. 运行分析并生成结果 ---
# 这里需要根据你的实际文件路径和SMR结果进行修改
# smr_mqtl_results <- read.table('smr_mQTLs_results.msmr', header=TRUE)
# smr_eqtl_results <- read.table('smr_eQTLs_results.msmr', header=TRUE)
# smr_pqtl_results <- read.table('smr_pQTLs_results.msmr', header=TRUE)

# mqtl_coloc_res <- run_smr_and_coloc("mQTL", smr_mqtl_results, 'mQTLs_coloc.txt', ..., mqtl_N)
# eqtl_coloc_res <- run_smr_and_coloc("eQTL", smr_eqtl_results, 'eQTLs_coloc.txt', ..., eqtl_N)
# ... (依此类推)

# --- 4. 可视化曼哈顿图 ---
# 这一部分依赖于外部R脚本 `refine_manh.R`，需要确保该文件存在于相应路径
# source("refine_manh.R") 

# 生成曼哈顿图的示例代码
# highlight_genes <- unique(mqtl_coloc_res$GeneSymbol) # 假设这是需要高亮的基因
# gwas_data_for_plot <- ... # 准备用于绘图的SMR结果数据
# p1 <- Manhattan_refine(Background_data = gwas_data_for_plot, highlight_names = highlight_genes, ...)
# ggsave('manhattan_plot.png', p1, width = 12, height = 6)