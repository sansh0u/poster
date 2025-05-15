library(topGO)
library(ggplot2)

gene2go <- read.table("gene2go.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene2go) <- c("gene", "go")

# 构造成 list: gene -> vector of GO terms
geneID2GO <- by(gene2go$go, gene2go$gene, function(x) as.character(x))

# === 2. 读取背景基因和目标基因 ===
background_genes <- readLines("all_genes_filtered.txt")
target_genes <- readLines("common_genes_filtered.txt")

# === 3. 构建 geneList（命名因子）===
geneList <- factor(as.integer(background_genes %in% target_genes))
names(geneList) <- background_genes

# === 4. 构建 topGO 对象 ===
GOdata <- new("topGOdata",
              ontology = "BP",  # 可选 "BP", "MF", "CC"
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

# === 5. 运行 Fisher 检验 ===
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# === 6. 输出显著 GO term ===
allRes <- GenTable(GOdata,
                   classicFisher = resultFisher,
                   orderBy = "classicFisher",
                   topNodes = 50)

write.table(allRes, file = "topGO_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 可选：只打印 p < 0.05 的
sigRes <- subset(allRes, as.numeric(classicFisher) < 0.05)
write.table(sigRes, file = "topGO_significant.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


