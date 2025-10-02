
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)


sample_ids <- gsub("d0$|d1$", "", grep("lane\\d+", names(Test_count_annot), value = TRUE))
sample_ids <- unique(sample_ids)

analyze_sample <- function(sample_id) {
  d0_col <- paste0(sample_id, "d0")
  d1_col <- paste0(sample_id, "d1")
  
  delta <- Test_count_annot[[d1_col]] - Test_count_annot[[d0_col]]
  
  category <- case_when(
    delta > 1 ~ "Upregulated",
    delta < -1 ~ "Downregulated",
    TRUE ~ "No change"
  )
  

  as.data.frame(table(category)) %>%
    rename(count = Freq) %>%
    mutate(sample = sample_id)
}


all_samples_summary <- lapply(sample_ids, analyze_sample) %>%
  bind_rows()


ranking <- all_samples_summary %>%
  filter(category == "Upregulated") %>%
  arrange(desc(count)) %>%
  pull(sample)


all_samples_summary$sample <- factor(all_samples_summary$sample, levels = ranking)
##Streght filtration column add##
Consensus_H <- Consensus_H %>%
  mutate(`All sample <10` = if_else(
    rowSums(select(., matches("^Ind\\d+d[01]")) > 10) == 
      length(select(., matches("^Ind\\d+d[01]"))),
    "Yes", "No"
  ))

##Deseq2##

library(DESeq2)
library(dplyr)

# 1. Dodaj kolumnę peak_id - załóżmy, że po prostu indeks wiersza
Consensus_H <- Consensus_H %>%
  mutate(peak_id = row_number())

# 2. Wybierz kolumny z counts i stwórz macierz counts
# Załóżmy, że counts to wszystkie kolumny z nazwami zawierającymi "Ind" i "d0" lub "d1"
count_cols <- grep("^Ind.*d[01]$", colnames(Consensus_H), value = TRUE)

count_data <- Consensus_H %>% select(all_of(count_cols))
count_data <- round(count_data)
rownames(count_data) <- Consensus_H$peak_id

# 3. Stwórz tabelę metadanych (colData) - grupa to 'd0' lub 'd1'
# Na podstawie nazwy kolumny wyciągamy dzień (ostatni znak kolumny)
coldata <- data.frame(
  sample = colnames(count_data),
  condition = ifelse(grepl("d0$", colnames(count_data)), "d0", "d1")
)
rownames(coldata) <- coldata$sample

# 4. Stwórz obiekt DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = coldata,
                              design = ~ condition)

# 5. Uruchom DESeq2
dds <- DESeq(dds)

# 6. Pobierz wyniki porównania d1 vs d0
res <- results(dds, contrast = c("condition", "d1", "d0"))

# 7. Dodaj wyniki do oryginalnej ramki
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

Consensus_H <- left_join(Consensus_H, res_df, by = "peak_id")

# Teraz Consensus_H ma kolumny z wynikami testu (log2FoldChange, pvalue, padj, itd.)

ggplot(all_samples_summary, aes(x = sample, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Liczba regulowanych regionów ATAC-seq w każdej próbce",
    x = "Próbka",
    y = "Liczba regionów"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_manual(values = c(
    "Upregulated" = "firebrick",
    "Downregulated" = "steelblue",
    "No change" = "grey80"
  ))
## RowSum > 90##
gr500kp <- GRanges(
  seqnames = RanjiniRNA$chromosome_name,
  ranges = IRanges(
    start = pmax(RanjiniRNA$start_position - 500000, 1),
    end = RanjiniRNA$end_position + 500000
  ))
gr500bp <- GRanges(
  seqnames = RanjiniRNA$chromosome_name,
  ranges = IRanges(
    start = pmax(RanjiniRNA$start_position - 500, 1),
    end = RanjiniRNA$end_position + 500
  ))
###Overlap Analysis with RNA###
gr_H <- GRanges(
  seqnames = Consensus_H$seqnames,
  ranges = IRanges(start = Consensus_H$start, end = Consensus_H$end)
)
hits_merged <- findOverlaps(gr_H, gr500bp)
Consensus_H$gr500bp <- FALSE
Consensus_H$gr500bp[unique(queryHits(hits_merged))] <- TRUE

# Sprawdzenie overlapu gr_H z gr_merged2
hits_merged2 <- findOverlaps(gr_H, gr500kp)
Consensus_H$gr500kp <- FALSE
Consensus_H$gr500kp[unique(queryHits(hits_merged2))] <- TRUE


colnames(consensus_matrix)
R
# Filtrujemy kolumny, które zawierają "laneXXXd0" lub "laneXXXd1"
# Załóżmy, że wszystkie próbki mają nazwę w formacie "laneXXXd0" lub "laneXXXd1"
consensus_matrix_filtered <- consensus_matrix[, grepl("lane\\d+d[01]", colnames(consensus_matrix))]

# Sprawdzamy, czy pozostały tylko poprawne próbki
colnames(consensus_matrix_filtered)

# Transponowanie danych
consensus_t <- t(consensus_matrix_filtered)

# Przygotowanie PCA
pca_res <- prcomp(consensus_t, scale. = TRUE)

# Przygotowanie metadanych
sample_names <- rownames(consensus_t)
condition <- ifelse(grepl("d0", sample_names), "d0", "d1")
pca_df <- data.frame(pca_res$x[, 1:2], condition = factor(condition))

# Wykres PCA
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = rownames(pca_df))) +
  geom_point(size = 3) +                            # Rysowanie punktów
  geom_text(size = 3, vjust = -1) +                 # Etykiety próbek (opcjonalnie)
  labs(title = "PCA ATAC-seq", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "top")                    # Dodanie legendy u góry



library(ggplot2)
library(dplyr)

# --- 2. Wybierz kolumny d0 i d1 ---
d0_cols <- grep("d0$", colnames(Consensus_F), value = TRUE)
d1_cols <- gsub("d0$", "d1", d0_cols)  # zakładamy, że pary istnieją

# Sprawdź poprawność
stopifnot(all(d1_cols %in% colnames(Consensus_F)))

# --- 3. Log2-transformacja ---
log_d0 <- log2(Consensus_F[, d0_cols] + 1)
log_d1 <- log2(Consensus_F[, d1_cols] + 1)

# --- 4. Oblicz różnice d1 - d0 ---
delta_mat <- log_d1 - log_d0

# --- 5. Ustaw uproszczone nazwy kolumn: laneXXX (bez d0/d1) ---
sample_names <- gsub("d0$", "", d0_cols)
colnames(delta_mat) <- sample_names

# --- 6. Standaryzacja po genach ---
scaled_delta <- t(scale(t(delta_mat)))  # Z-score po wierszach

# --- 7. Klasteryzacja hierarchiczna ---
dist_mat <- dist(t(scaled_delta))  # próbki w kolumnach
hc <- hclust(dist_mat, method = "ward.D2")

# Dendrogram
plot(hc, main = "Clustering based on expression changes (d1 - d0)", xlab = "", sub = "")

# --- 8. K-means clustering ---
set.seed(123)
k <- 3
kmeans_res <- kmeans(t(scaled_delta), centers = k)
cluster_labels <- factor(kmeans_res$cluster)

# --- 9. PCA ---
pca <- prcomp(t(scaled_delta), scale. = FALSE)
pca_df <- data.frame(pca$x[, 1:2])
pca_df$sample <- sample_names
pca_df$cluster <- cluster_labels

# --- 10. PCA wykres ---
ggplot(pca_df, aes(PC1, PC2, color = cluster, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA of Expression Changes (d1 - d0)",
       x = "PC1", y = "PC2") +
  theme_minimal()

###Cohort_figures###

##Manhatan Plot###

# 🔹 Kolumny dnia 0
day0_cols <- grep("d0$", colnames(Consensus_high), value = TRUE)
Consensus_high$avg_d0 <- rowMeans(Consensus_high[, day0_cols], na.rm = TRUE)

# 🔹 Mapowanie chromosomów
chrom_map <- setNames(1:24, c(paste0("chr", 1:22), "chrX", "chrY"))

# 🔹 Dane do wykresu
manhattan_df_day0 <- Consensus_high %>%
  mutate(
    SNP = paste0(seqnames, ":", start),
    CHR = chrom_map[seqnames],
    BP = start,
    P = avg_d0
  ) %>%
  filter(!is.na(CHR))

# 🔹 Kolory
my_colors <- rep(c("#0096d6", "#d76067"), length.out=length(unique(manhattan_df_day0$CHR)))

# 🔹 Styl graficzny – pogrubienie osi
par(lwd = 2)  # pogrubienie linii osi

# 🔹 Wykres
manhattan(manhattan_df_day0, chr="CHR", bp="BP", p="P", snp="SNP", logp=FALSE,
          ylab="Average Counts Day 0", xlab="Chromosome", ylim=c(0, 500), col=my_colors)
### Venn###
library(VennDiagram)
library(grid)

# Wersja z overlaps
overlaps <- findOverlaps(gr_merged, gr_consensus_f)

area1 <- length(gr_merged)
area2 <- length(gr_consensus)
cross_area <- length(unique(queryHits(overlaps)))

venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross_area,
  category = c("Merged", "Consensus_F"),
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cat.cex = 1.5,
  scaled = TRUE
)

grid.draw(venn.plot)


##Venn##
library(GenomicRanges)
library(VennDiagram)

# Zamiana tabel w GRanges
gr_L <- GRanges(seqnames = Consensus_L_filtered$seqnames,
                ranges = IRanges(start = Consensus_L_filtered$start,
                                 end = Consensus_L_filtered$end))

gr_M <- GRanges(seqnames = Consensus_M_filtered$seqnames,
                ranges = IRanges(start = Consensus_M_filtered$start,
                                 end = Consensus_M_filtered$end))

gr_H <- GRanges(seqnames = Consensus_H_filtered$seqnames,
                ranges = IRanges(start = Consensus_H_filtered$start,
                                 end = Consensus_H_filtered$end))

# Liczenie nakładania się
n_L <- length(gr_L)
n_M <- length(gr_M)
n_H <- length(gr_H)

n_LM <- length(findOverlaps(gr_L, gr_M, ignore.strand=TRUE))
n_LH <- length(findOverlaps(gr_L, gr_H, ignore.strand=TRUE))
n_MH <- length(findOverlaps(gr_M, gr_H, ignore.strand=TRUE))

# Nakładanie się wszystkich trzech
# findOverlaps zwraca obiekty typu Hits, więc trzeba policzyć unikalne L, które nakładają się na M i H
hits_LM <- findOverlaps(gr_L, gr_M, ignore.strand=TRUE)
hits_LH <- findOverlaps(gr_L, gr_H, ignore.strand=TRUE)
hits_L_M_H <- intersect(queryHits(hits_LM), queryHits(hits_LH))
n_LMH <- length(hits_L_M_H)

# Rysowanie Venn diagramu
venn.plot <- draw.triple.venn(
  area1 = n_L,
  area2 = n_M,
  area3 = n_H,
  n12 = n_LM,
  n23 = n_MH,
  n13 = n_LH,
  n123 = n_LMH,
  category = c("L", "M", "H"),
  fill = c("red", "green", "blue"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = 0
)

#Clustering
library(ggplot2)

# Przygotowanie danych
pca_plot_df <- as.data.frame(pca_res$x)
rownames(pca_plot_df) <- rownames(diff_samples)
pca_plot_df$cluster <- factor(kmeans_res$cluster, levels = 1:3, labels = c("high","mid","low"))
pca_plot_df$sample <- rownames(pca_plot_df)  # nazwy próbek do label

ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = sample), vjust = -0.6, size = 3, show.legend = FALSE) +  # show.legend=FALSE usuwa "a"
  scale_color_manual(
    values = c("high" = "#66c2a5", "mid" = "#fc8d62", "low" = "#8da0cb")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.8, color = "black")
  ) +
  ggtitle("kmeans clustering")
###RNA Analysis##
# Nazwy próbek
sample_names <- colnames(count_data)

# Warunek (d0 vs d1)
condition <- ifelse(grepl("^d0_", sample_names), "d0", "d1")

# Osoba (Ind01, Ind02 itd.)
individual <- gsub(".*Ind", "Ind", sample_names)

col_data <- data.frame(
  row.names   = sample_names,
  condition   = factor(condition, levels = c("d0", "d1")),
  individual  = factor(individual)  # pozwala na analizę sparowaną
)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData   = col_data,
                              design    = ~ individual + condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "d1", "d0"))

# shrink log2FC dla lepszej interpretacji
res <- lfcShrink(dds, coef = "condition_d1_vs_d0", type = "apeglm")

# sortowanie po padj
resOrdered <- res[order(res$padj), ]
# konwersja na dataframe z kolumną "Gene"
res_RNA_df <- as.data.frame(resOrdered)
res_RNA_df$Gene <- rownames(resOrdered)
res_RNA_df <- res_RNA_df %>% select(Gene, everything())

# zapis do Excela
write.xlsx(res_RNA_df, file = "res_RNA.xlsx", rowNames = FALSE)

RNA_M2 <- merge(RNA_M2, RNAres_M, 
                    by.x = "hgnc_symbol", 
                    by.y = "Gene", 
                    all.x = TRUE)
RNA_H2 <- merge(RNA_H2, RNAres_H, 
                    by.x = "hgnc_symbol", 
                    by.y = "Gene", 
                    all.x = TRUE)
RNA_L2 <- merge(RNA_L2, RNAres_L, 
                    by.x = "hgnc_symbol", 
                    by.y = "Gene", 
                    all.x = TRUE)
##Venn diagram##
library(VennDiagram)
library(grid)

# Filtrujemy geny po padj < 0.05
genes_L2 <- RNA_L2$hgnc_symbol[RNA_L2$padj < 0.05]
genes_M2 <- RNA_M2$hgnc_symbol[RNA_M2$padj < 0.05]
genes_H2  <- RNA_H2$hgnc_symbol[RNA_H2$padj < 0.05]

# Tworzymy listę dla VennDiagram
gene_list <- list(L2 = genes_L2, M2 = genes_M2, H2 = genes_H2)

# Tworzymy Venn diagram z poprawionymi kolorami i pozycją etykiet
venn.plot <- venn.diagram(
  x = gene_list,
  filename = NULL,
  fill = c("blue", "green", "red"),  # L2 - niebieski, M2 - zielony, H - czerwony
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("blue", "green", "red"),
  cat.pos = c(-20, 20, 180),        # przesunięcie etykiet aby nie nachodziły
  cat.dist = c(0.05, 0.05, 0.05),   # dystans od koła
  margin = 0.1
)

# Wyświetlamy diagram
grid.newpage()
grid.draw(venn.plot)



