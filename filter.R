Pachythelia <- read.csv("mkt_Pachythelia.tsv", header = TRUE, sep = "\t")
filtered_Pachythelia <- subset(Pachythelia, NI > 1 & alpha < 0 & pvalue < 0.05)

Taleporia <- read.csv("mkt_Taleporia.tsv", header = TRUE, sep = "\t")
filtered_Taleporia <- subset(Taleporia, NI > 1 & alpha < 0 & pvalue < 0.05)

Tischeria <- read.csv("mkt_Tischeria.tsv", header = TRUE, sep = "\t")
filtered_Tischeria <- subset(Tischeria, NI > 1 & alpha < 0 & pvalue < 0.05)

Bombyx <- read.csv("mkt_Bombyx.tsv", header = TRUE, sep = "\t")
filtered_Bombyx <- subset(Bombyx, NI > 1 & alpha < 0 & pvalue < 0.05)


common_genes <- Reduce(intersect, list(
  filtered_Pachythelia$file,
  filtered_Taleporia$file,
  filtered_Tischeria$file,
  filtered_Bombyx$file
))

all_genes <- unique(c(
  Pachythelia$file,
  Taleporia$file,
  Tischeria$file,
  Bombyx$file
))

writeLines(common_genes, "common_genes.txt")
writeLines(all_genes, "all_genes.txt")

Pachythelia <- read.csv("mkt_Pachythelia.tsv", header = TRUE, sep = "\t")
filtered_Pachythelia <- subset(Pachythelia, NI < 1 & alpha > 0 & pvalue < 0.05)

Taleporia <- read.csv("mkt_Taleporia.tsv", header = TRUE, sep = "\t")
filtered_Taleporia <- subset(Taleporia, NI < 1 & alpha > 0 & pvalue < 0.05)

Tischeria <- read.csv("mkt_Tischeria.tsv", header = TRUE, sep = "\t")
filtered_Tischeria <- subset(Tischeria, NI < 1 & alpha > 0 & pvalue < 0.05)

Bombyx <- read.csv("mkt_Bombyx.tsv", header = TRUE, sep = "\t")
filtered_Bombyx <- subset(Bombyx, NI < 1 & alpha > 0 & pvalue < 0.05)


common_genes1 <- Reduce(intersect, list(
  filtered_Pachythelia$file,
  filtered_Taleporia$file,
  filtered_Tischeria$file,
  filtered_Bombyx$file
))

all_genes <- unique(c(
  Pachythelia$file,
  Taleporia$file,
  Tischeria$file,
  Bombyx$file
))

writeLines(common_genes1, "common_genes1.txt")
writeLines(all_genes, "all_genes.txt")