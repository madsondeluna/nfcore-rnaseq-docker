# ============================================================
# Analise de Expressao Diferencial com DESeq2
# Entrada: matrizes de contagem do nf-core/rnaseq
# ============================================================

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)

# --------------------------------------------------
# 1. Carregar a matriz de contagem
# --------------------------------------------------
counts_file <- "results/star_salmon/salmon.merged.gene_counts.tsv"
counts_raw <- read_tsv(counts_file)

# A primeira coluna e o gene_id e a segunda e gene_name
gene_info <- counts_raw[, c("gene_id", "gene_name")]

# Montar a matriz de contagem (somente colunas numericas)
count_matrix <- as.matrix(counts_raw[, -c(1, 2)])
rownames(count_matrix) <- counts_raw$gene_id

# Arredondar para inteiros (necessario para DESeq2)
count_matrix <- round(count_matrix)

# --------------------------------------------------
# 2. Definir as condicoes experimentais
# --------------------------------------------------
# IMPORTANTE: ajustar conforme os nomes das suas colunas
# O padrao abaixo assume que amostras com "CONTROLE" no nome
# pertencem ao grupo controle, e as demais ao grupo tratamento.
sample_names <- colnames(count_matrix)

condition <- ifelse(grepl("CONTROLE", sample_names), "controle", "tratamento")

col_data <- data.frame(
    row.names = sample_names,
    condition = factor(condition, levels = c("controle", "tratamento"))
)

cat("Amostras e condicoes:\n")
print(col_data)

# --------------------------------------------------
# 3. Criar o objeto DESeq2 e executar a analise
# --------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition
)

# Filtrar genes com baixa contagem
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes apos filtragem:", nrow(dds), "\n")

# Executar a analise diferencial
dds <- DESeq(dds)

# Extrair os resultados
res <- results(dds, contrast = c("condition", "tratamento", "controle"))
res <- res[order(res$padj), ]

# --------------------------------------------------
# 4. Filtrar DEGs significativos
# --------------------------------------------------
degs <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

degs$gene_id <- rownames(degs)
degs <- merge(degs, gene_info, by = "gene_id", all.x = TRUE)
degs$regulacao <- ifelse(degs$log2FoldChange > 0, "UP", "DOWN")

cat("\n========================================\n")
cat("Total de DEGs encontrados:", nrow(degs), "\n")
cat("UP-regulados:", sum(degs$regulacao == "UP"), "\n")
cat("DOWN-regulados:", sum(degs$regulacao == "DOWN"), "\n")
cat("========================================\n")

# --------------------------------------------------
# 5. Salvar resultados
# --------------------------------------------------
dir.create("results/degs", showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(res), "results/degs/resultados_completos.csv")
write.csv(degs, "results/degs/degs_significativos.csv", row.names = FALSE)
write.csv(
    degs %>% filter(regulacao == "UP"),
    "results/degs/degs_upregulados.csv",
    row.names = FALSE
)
write.csv(
    degs %>% filter(regulacao == "DOWN"),
    "results/degs/degs_downregulados.csv",
    row.names = FALSE
)

# --------------------------------------------------
# 6. Visualizacoes
# --------------------------------------------------

# -- Volcano Plot --
pdf("results/degs/volcano_plot.pdf", width = 10, height = 8)
res_df <- as.data.frame(res) %>%
    mutate(
        significant = ifelse(
            !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
            ifelse(log2FoldChange > 0, "UP", "DOWN"),
            "NS"
        )
    )

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("UP" = "#D73027", "DOWN" = "#4575B4", "NS" = "grey60")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10 p-value",
        color = "Regulacao"
    ) +
    theme_minimal()
dev.off()

# -- PCA --
pdf("results/degs/pca_plot.pdf", width = 8, height = 6)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") +
    theme_minimal() +
    labs(title = "PCA - Amostras")
dev.off()

# -- Heatmap dos top 50 DEGs --
if (nrow(degs) > 0) {
    n_top <- min(50, nrow(degs))
    pdf("results/degs/heatmap_top_degs.pdf", width = 10, height = 12)
    top_genes <- head(degs$gene_id, n_top)
    mat <- assay(vsd)[top_genes, ]
    mat <- mat - rowMeans(mat)
    pheatmap(
        mat,
        annotation_col = as.data.frame(colData(dds)["condition"]),
        scale = "row",
        show_rownames = TRUE,
        main = paste("Top", n_top, "DEGs")
    )
    dev.off()
}

cat("\nAnalise concluida. Resultados salvos em results/degs/\n")
cat("Arquivos gerados:\n")
cat("  - results/degs/resultados_completos.csv\n")
cat("  - results/degs/degs_significativos.csv\n")
cat("  - results/degs/degs_upregulados.csv\n")
cat("  - results/degs/degs_downregulados.csv\n")
cat("  - results/degs/volcano_plot.pdf\n")
cat("  - results/degs/pca_plot.pdf\n")
cat("  - results/degs/heatmap_top_degs.pdf\n")
