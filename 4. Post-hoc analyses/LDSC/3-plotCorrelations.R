# read in corr matrix
library(tidyverse)
library(corrplot)

cor.dir <- ""

cor.mat <- read_tsv(file.path(cor.dir, "corr-matrix_r3.1.txt"))
pval.mat <- read_tsv(file.path(cor.dir, "corr-matrix_p3.1.txt"))
  
names <- colnames(cor.mat)
names <- colnames(pval.mat)

row.names(cor.mat) <- names
as.matrix(cor.mat)

row.names(pval.mat) <- names
as.matrix(pval.mat)

corrplot(as.matrix(cor.mat), p.mat = as.matrix(pval.mat), 
         method = "square", type = "upper", diag = FALSE, tl.col = 'black', 
         tl.cex = 1.3,
         #tl.srt = 45,
         sig.level = c(0.0024, 0.05), pch.cex = 2,
         insig = 'label_sig', pch.col = 'white',
         col = COL2('RdBu', 10)) 
