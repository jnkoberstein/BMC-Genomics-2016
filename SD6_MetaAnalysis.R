library(GEOquery)
library(RColorBrewer)
library(RUVSeq)
library(mouse4302.db)
library(mogene21sttranscriptcluster.db)
library(biomaRt)
library(stringr)
library(xlsx)
library(plotrix)

notUnique <- function(x) x %in% x[duplicated(x)]
map <- "ensembl_gene_id"

load("expression.rda")
gerstner <- eData[,grep("CC6|SD", colnames(eData))]
colnames(gerstner) <- c(paste0("CC6_Gerstner",1:7), paste0("SD6_Gerstner",1:7))
gerstner <- as.matrix(gerstner)
mogene2 <- rownames(gerstner)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
mapping_mogene2 <- getBM(attributes = c("affy_mogene_2_1_st_v1", map), filters = "affy_mogene_2_1_st_v1",
                         values = mogene2, mart = mart)
temp <- merge(mapping_mogene2, gerstner, by.x = "affy_mogene_2_1_st_v1", by.y = "row.names")
temp <- temp[which(!notUnique(temp[,1])),]
temp <- temp[which(!notUnique(temp[,2])),]

gse9444 <- getGEO('GSE9444', GSEMatrix = T)[[1]]
maret <- gse9444[,intersect(grep("ZT 6", pData(phenoData(gse9444))[,1]), grep("B6", pData(phenoData(gse9444))[,1]))]
maret <- as.matrix(maret[,-grep("L_", pData(phenoData(maret))[,1])])

gse6514 <- getGEO('GSE6514', GSEMatrix = T)[[1]]
mackiewicz <- gse6514[,grep("cerebral_cortex_6hrs", pData(phenoData(gse6514))[,1])]
mackiewicz <- as.matrix(log2(exprs(mackiewicz)))

gse33302 <- getGEO('GSE33302', GSEMatrix = T)[[1]]
peixoto <- exprs(gse33302)

combined <- combine(maret, mackiewicz, peixoto)
colnames(combined) <- c(paste0("CC6_Maret1",1:3), paste0("SD6_Maret1",1:3), paste0("CC6_Maret2",1:3), paste0("SD6_Maret2",1:3), 
                        paste0("CC6_Mack",1:5), paste0("SD6_Mack",1:5), paste0("CC6_Peixoto",1:9), paste0("SD6_Peixoto",1:8))
combined <- combined[,c(1:3, 7:9, 13:17, 23:31, 4:6, 10:12, 18:22, 32:39)]
combined <- as.matrix(combined)

mo430 <- rownames(combined)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
mapping_mo430 <- getBM(attributes = c("affy_mouse430_2", map), filters = "affy_mouse430_2", 
                       values = mo430, mart = mart)
temp2 <- merge(mapping_mo430, combined, by.x = "affy_mouse430_2", by.y = "row.names")
temp2 <- temp2[-grep("a_at|x_at|s_at", temp2$affy_mouse430_2),]
temp2 <- temp2[which(!notUnique(temp2[,1])),]
temp2 <- temp2[which(!notUnique(temp2[,2])),]
merged <- merge(temp, temp2, by.x = map, by.y = map, all.x = FALSE, all.y = FALSE)
merged <- merged[,-grep("affy", colnames(merged))]
rownames(merged) <- merged[,1]
merged <- merged[,-1]
merged <- merged[,sort(colnames(merged))]
#SANITY CHECK
dim(merged[which(!notUnique(rownames(merged))),]) == dim(merged)

background <- rownames(merged)
write.table(background, file = "background.txt", sep = "\t", quote = F, row.names = F, col.names = F)

merged <- as.matrix(merged)
filter <- apply(merged, 1, function(x) length(x[x < 4]) <= dim(merged)[2]/2)
merged <- merged[filter,]
x <- as.factor(str_split_fixed(colnames(merged), "_", n = 2)[,1])
names <- levels(x)
sd <- as.matrix(merged)
sd <- normalizeBetweenArrays(sd, method = "quantile")

symbols <- c(21,22,24,23)
lab <- factor(c(rep("Gerstner",7), rep("Mackiewicz",5), rep("Maret",6), rep("Peixoto",9),rep("Gerstner",7), 
                 rep("Mackiewicz",5), rep("Maret",6), rep("Peixoto",8)))
colors <- brewer.pal(3, "Set2")
par(mar = c(7,4,1,1))
plotRLE(sd, outline = FALSE, col = colors[x], las = 2, ylim = c(-0.5, 0.5))
svg(filename = "Figures/Figure1_RMA_SD6M_PCA_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
plotPCA(sd, col = c("red","black","black","black")[lab], lwd = 2, bg = colors[x], isLog = T, labels = F, pch = symbols[lab], cex = 3, cex.lab = 2.5, cex.axis = 1.5)
legend("bottomleft", c(levels(lab), "Sleep Deprivation", "Circadian Control"), pch = c(21,22,24,23,19,19), 
       col = c(rep("black", 4), "#66C2A5", "#FC8D62"), pt.cex = 1.5)
dev.off()

des <- model.matrix(~0 + x)
colnames(des) <- names
f <- lmFit(sd, des)
contrast.matrix <- makeContrasts(SD6 - CC6, levels = des)
f <- contrasts.fit(f, contrast.matrix)
f <- eBayes(f)
preRUV <- topTable(f, coef = 1,  adjust.method = "fdr", p.value = 1, number = Inf)
all.p <- preRUV$P.Val
hist(all.p, breaks = 1000, xlab = "p-value", main = NULL, xlim = c(0,1), ylim = c(0,500))

makeRow <- function(y) { c(grep(y, x), rep(-1, 27 - length(grep(y, x)))) }
groups <- matrix(data = unlist(lapply(names, makeRow)), nrow = length(names), byrow = TRUE)

n = 5
controls = rownames(sd)
s <- RUVs(sd, cIdx = controls, k = n, scIdx = groups, round = FALSE, isLog = T)

par(mar = c(7,4,1,1))
plotRLE(s$normalizedCounts, outline = FALSE, col = colors[x], las = 2)
svg(filename = "Figures/Figure1_RUV_SD6M_PCA_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
plotPCA(s$normalizedCounts, col = c("red","black","black","black")[lab], lwd = 2, bg = colors[x], isLog = T, labels = F, pch = symbols[lab], cex = 3, cex.lab = 2.5, cex.axis = 1.5)
legend("topleft", c(levels(lab), "Sleep Deprivation", "Circadian Control"), pch = c(21,22,24,23,19,19), 
       col = c(rep("black", 4), "#66C2A5", "#FC8D62"), pt.cex = 1.5)
dev.off()

##################################################################################################
#library(pca3d)
# library(rgl)
# pc <- prcomp(t(sd))
# lab <- factor(c(rep("Mack",5), rep("Maret1",3), rep("Maret2",3), rep("gerstner",7), rep("Mack",5), 
#                 rep("Maret1",3), rep("Maret2",3), rep("gerstner",7)))
# 
# pca3d(pc, comp = 1:3, show.labels = rep("",36), show.plane = T, group = lab, radius = 0.5, legend = "topleft",
#       show.shadows = T, show.centroids = T, show.group.labels = F, show.shapes = T, show.ellipses = T)
# rgl.snapshot("preRUV_3dPCA.png", fmt = "png", top = TRUE )
# 
# pc2 <- prcomp(t(s$normalizedCounts))
# trt <- factor(c(rep("CC",18),rep("SD",18)))
# 
# pca3d(pc2, comp = 1:3, show.labels = rep("",36), show.plane = F, group = trt, radius = 0.5, legend = "topleft",
#       show.shadows = T, show.centroids = T, show.group.labels = F, show.ellipses = T, show.shapes = T, new = T)
# snapshotPCA3d("postRUV_3dPCA.png")
##################################################################################################

des <- model.matrix(~0 + x + s$W)
colnames(des) <- c(names, paste0("W", 1:n))
f <- lmFit(sd, des)
contrast.matrix <- makeContrasts(SD6 - CC6, levels = des)
f <- contrasts.fit(f, contrast.matrix)
f <- eBayes(f)
postRUV <- topTable(f, coef = 1,  adjust.method = "fdr", p.value = 1, number = Inf)
all.p <- postRUV$P.Val
hist(all.p, breaks = 1000, xlab = "p-value", main = NULL, xlim = c(0,1), ylim = c(0,800), cex.lab = 2.5, cex.axis = 2)

pos.control.meta <- c('Arc', 'Dusp1', 'Egr1', 'Egr2', 'Homer1', 'Hspa5', 'Klf10', 'Nr4a1', 'Nr4a3', 'Ppm2c', 'Ptgs2', 'Rbm3',
                        'Sult1a1', 'Crh')

pos.control.maret <- c('Arc','Calr','Dio2','Fos','Hsp5a','Hspa1b','Npas2','Nr4a1','P4ha1','Per1','Per2','Rbm3','Sult1a1',
                       'Tfrc','Tipin','Dbp','Bdnf','Cirbp','Crh','Egr1','Egr3','Nptx2','Sgk','Slc2a1','Tmem10',
                       'Tsc22d3','Homer1')

pos.control <- union(pos.control.meta, pos.control.maret)
pos.control.genes <- unname(unlist(getBM(filters = "mgi_symbol", attributes = map, 
                                  values = pos.control, mart = mart)))

###PLOT
svg(filename = "Figures/Figure1_RMA_SD6M_volcano_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(5,5,3,1))
with(preRUV, plot(logFC, -log10(P.Value), pch=20, ylab="-log10(p-value)", xlab="logFC", cex.lab = 2.5, cex.axis = 1.5, ylim = c(0,18), cex = 2))
with(subset(preRUV, adj.P.Val < 0.01 ), points(logFC, -log10(P.Value), pch=20, col="dodgerblue", cex = 2))
with(subset(preRUV, rownames(preRUV) %in% pos.control.genes), points(logFC, -log10(P.Value), pch=1, lwd=4, col="red", cex = 1.5))
title(main = "13% Positive Controls", cex.main = 2)
dev.off()

svg(filename = "Figures/Figure1_RUV_SD6M_volcano_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(5,5,3,1))
with(postRUV, plot(logFC, -log10(P.Value), pch=20, ylab="-log10(p-value)", xlab="logFC", cex.lab = 2.5, cex.axis = 1.5, ylim = c(0,18), cex = 2))
with(subset(postRUV, adj.P.Val < 0.01 ), points(logFC, -log10(P.Value), pch=20, col="dodgerblue", cex = 2))
with(subset(postRUV, rownames(postRUV) %in% pos.control.genes), points(logFC, -log10(P.Value), pch=1, lwd=4, col="red", cex = 1.5))
title(main = "100% Positive Controls", cex.main = 2)
dev.off()


###Broken axis p-value plot
##################################################################################################################
svg(filename = "Figures/Figure1_RMA_SD6M_pvalue_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
h <- hist(preRUV$P.Value, breaks = 100)
h$counts[h$counts > 400] <- h$counts[h$counts > 400] - 900
plot(h, ylim = c(0,500), yaxt = "n", xlab = "p-value", main = NULL, xlim = c(0,1), cex.lab = 2.5, cex.axis = 1.5, col = "grey")
axis(side = 2, at = c(0,100,200,300,400,500), labels = c(0,100,200,300,1300,1400), cex.axis = 1.5)
axis.break(axis = 2, breakpos = 375, style = "slash", brw = 0.03)
dev.off()

svg(filename = "Figures/Figure1_RUV_SD6M_pvalue_hippocampus.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
h <- hist(postRUV$P.Value, breaks = 100)
h$counts[h$counts > 400] <- h$counts[h$counts > 400] - 900
plot(h, ylim = c(0,500), yaxt = "n", xlab = "p-value", main = NULL, xlim = c(0,1), cex.lab = 2.5, cex.axis = 1.5, col = "grey")
axis(side = 2, at = c(0,100,200,300,400,500), labels = c(0,100,200,300,1300,1400), cex.axis = 1.5)
axis.break(axis = 2, breakpos = 375, style = "slash", brw = 0.03)
dev.off()

###Supplementary Table
tab <- topTable(f, coef = 1,  adjust.method = "fdr", p.value = 1, number = Inf)
gene <- rownames(tab)
nt <- getBM(filters = map, attributes = c(map, "mgi_symbol"), values= gene, mart= mart)
tab <- merge(tab, nt, by.x = "row.names", by.y = map)
tab$mogene2_id <- mapping_mogene2$affy_mogene_2_1_st_v1[match(tab$Row.names, mapping_mogene2$ensembl_gene_id)]
tab$mo430_id <- mapping_mo430$affy_mouse430_2[match(tab$Row.names, mapping_mo430$ensembl_gene_id)]
tab <- tab[, c(1, 9, 10, 8, 2:7)]
colnames(tab)[1:4] <- c("Ensembl ID", "Mouse Gene 2.1 ST probeset ID", "Mouse Genome 430 2.0 probeset ID", "MGI Symbol")
write.xlsx(tab, file = "Tables/Additional File 2 full.xlsx", sheetName = "Maret_Mackiewicz_Gerstner_Peixoto", row.names = F)

