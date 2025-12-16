setwd('./Mendelian')

mendelian_gene = read.csv('mendelian_gene.csv')
mendelian_gene = unique(mendelian_gene) 
cigma = read.csv('cis.jk.csv')
cigma$spe = cigma$"free.v" / (cigma$"free.v" + cigma$"free.hom_g2")

cigma_men = cigma[which(cigma$gene %in% mendelian_gene$gene),]
cigma_nonmen = cigma[which(!cigma$gene %in% mendelian_gene$gene),]
t.test(cigma_men$"free.v", cigma_nonmen$"free.v") 
t.test(cigma_men$"free.hom_g2", cigma_nonmen$"free.hom_g2")
t.test(cigma_men$"spe", cigma_nonmen$"spe") 

# remove MHC genes (the same as AMM)
cigma_noMHC = fread('cigma_gene_noMHC10178.txt')
cigma_men = cigma_men[which(cigma_men$gene %in% cigma_noMHC$x),]
cigma_nonmen = cigma_nonmen[which(cigma_nonmen$gene %in% cigma_noMHC$x),]

t.test(cigma_men$"free.v", cigma_nonmen$"free.v") 
t.test(cigma_men$"free.hom_g2", cigma_nonmen$"free.hom_g2")
t.test(cigma_men$"spe", cigma_nonmen$"spe") 