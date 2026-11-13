# STEP 1

install.packages(c("data.table", "R.utils"))   # run once
library(data.table)

clin <- fread("TCGA-BRCA.clinical.tsv.gz")
expr <- fread("TCGA-BRCA.star_fpkm.tsv.gz")    # large file: may take 1-2 minutes
map  <- fread("gencode.v36.annotation.gtf.gene.probemap")

dim(clin); dim(expr); dim(map)   # rows and columns for each
names(map)                       # mapping column names
head(map)                        # view structure
names(clin)                      # all clinical columns
expr[1:5, 1:4]                   # corner of the expression matrix

--------------------------------------------------------------------------------

# STEP 2

genes <- c("MKI67","BRCA1","TP53","CCNB1","CDK1","ERBB2","ESR1","FOXA1","EGFR","PGR")

# A) find IDs for the 10 genes, extract their rows and transpose the table
sel <- map[gene %in% genes, .(id, gene)]
print(sel)                          # must have 10 rows, one per gene

x <- as.data.frame(expr[Ensembl_ID %in% sel$id])
rownames(x) <- x$Ensembl_ID
x$Ensembl_ID <- NULL
m <- t(as.matrix(x))                # now: samples on rows, genes on columns
colnames(m) <- sel$gene[match(colnames(m), sel$id)]
dat <- data.frame(sample = rownames(m), m, row.names = NULL)
dim(dat); head(dat, 3)

# B) clinical quality checks
sum(duplicated(clin$sample))                  # how many rows repeat a sample
sum(dat$sample %in% clin$sample)              # how many dat samples are found in clinical data
table(clin$sample_type.samples, useNA = "ifany")   # sample types
summary(clin$age_at_index.demographic)        # age: minimum, maximum, missing values

--------------------------------------------------------------------------------

# STEP 3

# merge age and sample type with gene expression
clin2 <- as.data.frame(clin[, .(sample,
                                  sample_type_gdc = sample_type.samples,
                                  age = age_at_index.demographic)])
d <- merge(dat, clin2, by = "sample")

# patient = first 12 characters of the sample barcode
d$patient <- substr(d$sample, 1, 12)

# keep only primary tumors and normal tissue, and create group variable
d <- d[d$sample_type_gdc %in% c("Primary Tumor", "Solid Tissue Normal"), ]
d$sample_type <- ifelse(d$sample_type_gdc == "Primary Tumor", "Tumor", "Normal")

# checks
dim(d)
table(d$sample_type)
sum(duplicated(paste(d$patient, d$sample_type)))   # patients with multiple samples in the same group
length(intersect(d$patient[d$sample_type == "Tumor"],
                 d$patient[d$sample_type == "Normal"]))   # patients with both (paired test)

# if there are duplicates, keep a single sample per patient and group
d <- d[!duplicated(paste(d$patient, d$sample_type)), ]
table(d$sample_type)

# save the clean dataset
write.csv(d, "brca_clean.csv", row.names = FALSE)

# example table for slide 5 (4 tumors and 4 normals)
tab <- d[, c("patient", "sample_type", "age", "MKI67", "ESR1", "FOXA1", "EGFR")]
tab[, 4:7] <- round(tab[, 4:7], 2)
rbind(head(tab[tab$sample_type == "Tumor", ], 4),
      head(tab[tab$sample_type == "Normal", ], 4))

--------------------------------------------------------------------------------

# STEP 4

# if R was closed, reload: d <- read.csv("brca_clean.csv")
genes <- c("MKI67","BRCA1","TP53","CCNB1","CDK1","ERBB2","ESR1","FOXA1","EGFR","PGR")

# function computing a statistic per group and rounding to 2 decimal places
st <- function(f) {
  r <- aggregate(d[, genes], by = list(group = d$sample_type), FUN = f)
  r[, -1] <- round(r[, -1], 2)
  r
}
st(mean); st(median); st(sd); st(IQR)

# histograms for each gene (all samples)
png("hist_genes.png", width = 1800, height = 900, res = 150)
par(mfrow = c(2, 5))
for (g in genes) hist(d[[g]], main = g, xlab = "log2(FPKM+1)", col = "grey80", breaks = 30)
dev.off()

# tumor vs normal boxplots for each gene
png("box_tumor_normal.png", width = 1800, height = 900, res = 150)
par(mfrow = c(2, 5))
for (g in genes) boxplot(d[[g]] ~ d$sample_type, main = g, xlab = "",
                         ylab = "log2(FPKM+1)", col = c("lightblue", "salmon"))
dev.off()

--------------------------------------------------------------------------------
  
# STEP 5

tum <- d[d$sample_type == "Tumor", ]

# correlation among the 10 genes
cm <- cor(tum[, genes])
round(cm, 2)
png("corr_heatmap.png", width = 1200, height = 1000, res = 150)
heatmap(cm, symm = TRUE, margins = c(8, 8))
dev.off()

# correlation of each gene with age
sapply(genes, function(g) round(cor(tum[[g]], tum$age, use = "complete.obs"), 2))

# MKI67 against other genes (excluding CCNB1 and CDK1) and against age
pred <- c("BRCA1","TP53","ERBB2","ESR1","FOXA1","EGFR","PGR","age")
png("scatter_mki67.png", width = 1800, height = 900, res = 150)
par(mfrow = c(2, 4))
for (p in pred) plot(tum[[p]], tum$MKI67, main = p, xlab = p,
                     ylab = "MKI67", pch = 16, cex = 0.4)
dev.off()
  
--------------------------------------------------------------------------------

# STEP 6

fmt <- function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))
tab_desc <- data.frame(
  gene   = genes,
  Normal = sapply(genes, function(g) fmt(d[[g]][d$sample_type == "Normal"])),
  Tumor  = sapply(genes, function(g) fmt(d[[g]][d$sample_type == "Tumor"])),
  row.names = NULL)
tab_desc
write.csv(tab_desc, "table_descriptive.csv", row.names = FALSE)
