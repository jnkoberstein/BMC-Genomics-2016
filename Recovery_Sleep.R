library(RColorBrewer)
library(RUVSeq)
library(mogene21sttranscriptcluster.db)
library(stringr)
library(gplots)
library(plotrix)
library(xlsx)
library(dendextend)
library(ggplot2)
library(reshape2)
load("expression.rda")
notUnique <- function(x) x %in% x[duplicated(x)]

###ID MAP
#################################################################################################################
x <- mogene21sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
map_symbol <- as.data.frame(x[mapped_probes])
colnames(map_symbol) <- c("affy_mogene_2_1_st_v1", "symbol")

x <- mogene21sttranscriptclusterENSEMBL
mapped_probes <- mappedkeys(x)
map_ensembl <- as.data.frame(x[mapped_probes])
colnames(map_ensembl) <- c("affy_mogene_2_1_st_v1", "ensembl")

###POSITIVE CONTROLS
##################################################################################################################
pos.control.meta <- c('Arc','Dusp1','Egr1','Egr2','Homer1','Klf10','Nr4a1','Nr4a3','Pdp1','Ptgs2','Rbm3','Sult1a1','Crh')
pos.control.maret <- c('Arc','Calr','Dio2','Fos','Hspa5','Hspa1b','Npas2','Nr4a1','P4ha1','Per1','Per2','Rbm3','Sult1a1',
                       'Tfrc','Tipin','Dbp','Bdnf','Cirbp','Crh','Egr1','Egr3','Nptx2','Sgk1','Slc2a1','Opalin',
                       'Tsc22d3','Homer1')
pos.sd <- union(pos.control.meta, pos.control.maret)
pos.sd.genes <- map_symbol[match(pos.sd, map_symbol[,2]),]

pos.rs <- c("Bdnf","Egr3","Fosl2","Hsp90b1","Hspa5","Junb","Nptx2")
pos.rs.genes <- map_symbol[match(pos.rs, map_symbol[,2]),]

###DATA PREP
##################################################################################################################
eData <- as.matrix(eData)
eData <- eData[,-grep("RS1_3", colnames(eData))]
eData <- eData[,-grep("CC7_3", colnames(eData))]
eData <- eData[,c(1:34, 49:65, 35:48)]

filter <- apply(eData, 1, function(x) length(x[x < 4]) <= dim(eData)[2]/2)
filtered <- eData[filter,]

x <- as.factor(str_split_fixed(colnames(filtered), "_", n = 2)[,1])
x <- factor(x,levels(x)[c(1,3,4,5,2,10,6,7,8,9)])
names <- levels(x)
dat <- as.matrix(filtered)

symbols <- c(21,22,24)
symbols2 <- c(19,17,15)
trt <- factor(c(rep("CC", 34), rep("SD", 7), rep("RS", 24)))
colfunc.sd <- colorRampPalette(c("red", "darkred"))
colfunc.rs <- colorRampPalette(c("dodgerblue", "darkblue"))
colfunc.cc <- colorRampPalette(c("green", "darkgreen"))
colors <- c(colfunc.cc(5), colfunc.sd(1), colfunc.rs(4))
background <- rownames(dat)
background <- map_ensembl$ensembl_id[match(background, map_ensembl$probe_id)]
write.table(background, file = "Tables/background.txt", sep = "\t", quote = F, row.names = F, col.names = F)

###PRE RUV TEST
##################################################################################################################
par(mar = c(7,4,1,1))
plotRLE(dat, outline = FALSE, col = colors[x], las = 2, ylim = c(-0.5, 0.5))

svg(filename = "Figures/Figure2_RMA_PCA.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
plotPCA(dat, bg = colors[x], isLog = T, labels = F, pch = symbols[trt], cex = 3, cex.lab = 2.5, cex.axis = 1.5)
legend("bottomright", levels(x), col = colors, pch = symbols2[c(rep(1,5), 2, rep(3,4))],  pt.cex = 1.5, bty = "n")
dev.off()

des1 <- model.matrix(~0 + x)
colnames(des1) <- names
f1 <- lmFit(dat, des1)
contrast.matrix <- makeContrasts(SD - CC6, RS1 - CC7, RS2 - CC8, RS3 - CC8, RS6 - CC11, levels = des1)
f1 <- contrasts.fit(f1, contrast.matrix)
f1 <- eBayes(f1)

###RUV
##################################################################################################################
makeRow <- function(y) { c(grep(y, x), rep(-1, 7 - length(grep(y, x)))) }
groups <- matrix(data = unlist(lapply(names, makeRow)), nrow = length(names), byrow = TRUE)

k = 6
controls = rownames(dat)
s <- RUVs(dat, cIdx = controls, k = k, scIdx = groups, round = FALSE, isLog = T)

###POST RUV TEST
##################################################################################################################
par(mar = c(7,4,1,1))
plotRLE(s$normalizedCounts, outline = FALSE, ylim = c(-0.5, 0.5), col = colors[x], las = 2)

svg(filename = "Figures/Figure2_RUV_PCA.svg", width = 4, height = 4, pointsize = 7)
par(mar = c(6,6,1,1))
plotPCA(s$normalizedCounts, bg = colors[x], isLog = T, labels = F, pch = symbols[trt], cex = 3, cex.lab = 2.5, cex.axis = 1.5)
legend("topleft", levels(x), col = colors, pch = symbols2[c(rep(1,5), 2, rep(3,4))],  pt.cex = 1.5, bty = "n")
dev.off()

des2 <- model.matrix(~0 + x + s$W)
colnames(des2) <- c(names, paste0("W", 1:k))
f2 <- lmFit(dat, des2)
contrast.matrix <- makeContrasts(SD - CC6, RS1 - CC7, RS2 - CC8, RS3 - CC8, RS6 - CC11, levels = des2)
f2 <- contrasts.fit(f2, contrast.matrix)
f2 <- eBayes(f2)

###PLOT FUNCTIONS
##################################################################################################################
plot_volcano <- function(tests, pos, main = NULL, y) {
  par(mar = c(5,5,3,1))
  with(tests, plot(logFC, -log10(P.Value), pch=20, ylab="-log10(p-value)", xlab="logFC", ylim = c(0,y), cex.lab = 2.5, cex.axis = 1.5, cex = 2))
  with(subset(tests, adj.P.Val<0.01 ), points(logFC, -log10(P.Value), pch=20, col="dodgerblue", cex = 2))
  with(subset(tests, rownames(tests) %in% pos), points(logFC, -log10(P.Value), pch=1, lwd=4, cex = 1.5, col="red"))
  title(main = main, cex.main = 2)
}

plot_pvalue <- function(tests, pos, main = NULL, brk, max) { 
  par(mar = c(6,6,1,1))
  h <- hist(tests$P.Value, breaks = 100, ylab = "")
  h$counts[h$counts > brk] <- h$counts[h$counts > brk] - ((h$counts[h$counts > brk] - brk) - ((h$counts[h$counts > brk] - brk) %% 100))
  plot(h, ylim = c(0, brk + 100), yaxt = "n", xlab = "p-value", ylab = "", main = NULL, xlim = c(0,1), cex.lab = 2.5, cex.axis = 1.5, col = "grey")
  axis(side = 2, at = c(seq(0, brk, by = 200), brk, brk + 100), labels = c(seq(0, brk - 100, by = 200), max - 100, max), las = 1, cex.axis = 1.5, ylab = "")
  title(ylab = "Frequency", mgp=c(4,1,0), cex.lab = 2.5)
  axis.break(axis = 2, breakpos = brk - 50, style = "slash", brw = 0.03)
}

###RMA PLOTS
##################################################################################################################
toptable <- function(x, f, p) { topTable(f, coef = x, adjust.method = "fdr", p.value = p, number = Inf, confint = F) }
tests <- c("SD" = 1, "RS1" = 2, "RS2" = 3, "RS3" = 4, "RS6" = 5)
RMA.test.all <- lapply(tests, toptable, f = f1, p = 1)
RUV.test.all <- lapply(tests, toptable, f = f2, p = 1)

list2env(RMA.test.all, .GlobalEnv)
svg(filename = "Figures/AdditionalFile3_RMA_SD_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(SD, pos.sd.genes[,1], brk = 2700, max = 3300)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_SD_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(SD, pos.sd.genes[,1], main = "72% Positive Controls", y = 15)
dev.off()
svg(filename = "Figures/Figure2_RMA_RS1_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS1, pos.rs.genes[,1], main = "29% Positive Controls", y = 12)
dev.off()
svg(filename = "Figures/Figure2_RMA_RS1_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS1, pos.rs.genes[,1], brk = 700, max = 1300)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS2_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS2, pos.rs.genes[,1], main = "29% Positive Controls", y = 12)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS2_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS2, pos.rs.genes[,1], brk = 900, max = 1600)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS3_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS3, pos.rs.genes[,1], main = "29% Positive Controls", y = 15)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS3_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS3, pos.rs.genes[,1], brk = 900, max = 1500)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS6_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS6, pos.rs.genes[,1], main = "0% Positive Controls", y = 10)
dev.off()
svg(filename = "Figures/AdditionalFile3_RMA_RS6_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS6, pos.rs.genes[,1], brk = 1300, max = 1900)
dev.off()

###RUV PLOTS
########################################################################################################################
list2env(RUV.test.all, .GlobalEnv)
svg(filename = "Figures/AdditionalFile3_RUV_SD_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(SD, pos.sd.genes[,1], brk = 2700, max = 3300)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_SD_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(SD, pos.sd.genes[,1], main = "91% Positive Controls", y = 15)
dev.off()
svg(filename = "Figures/Figure2_RUV_RS1_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS1, pos.rs.genes[,1], main = "86% Positive Controls", y = 12)
dev.off()
svg(filename = "Figures/Figure2_RUV_RS1_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS1, pos.rs.genes[,1],  brk = 700, max = 1300)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS2_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS2, pos.rs.genes[,1], main = "86% Positive Controls", y = 15)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS2_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS2, pos.rs.genes[,1], brk = 900, max = 1600)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS3_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS3, pos.rs.genes[,1], main = "57% Positive Controls", y = 15)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS3_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS3, pos.rs.genes[,1], brk = 900, max = 1500)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS6_volcano.svg", width = 4, height = 4, pointsize = 7)
plot_volcano(RS6, pos.rs.genes[,1], main = "29% Positive Controls", y = 10)
dev.off()
svg(filename = "Figures/AdditionalFile3_RUV_RS6_pvalue.svg", width = 4, height = 4, pointsize = 7)
plot_pvalue(RS6, pos.rs.genes[,1], brk = 1300, max = 1900)
dev.off()

###Supplementary Table 2
########################################################################################################################
RUV.test.sig <- lapply(tests, toptable, f = f2, p = 0.01)
addcolumn <- function(dat) { cbind(dat, "MGI Symbol" = map_symbol[match(rownames(dat), map_symbol[,1]),2], 
                                   "Ensembl ID" = map_ensembl[match(rownames(dat), map_ensembl[,1]),2])[,c(8,7,1:6)]}
RUV.test.sig.anno <- lapply(RUV.test.sig, addcolumn)
list2env(RUV.test.sig.anno, .GlobalEnv)

write.xlsx(SD, file = "Tables/Additional File 4.xlsx", sheetName = "SD", row.names = T)
write.xlsx(RS1, file = "Tables/Additional File 4.xlsx", sheetName = "RS1", append = TRUE, row.names = T)
write.xlsx(RS2, file = "Tables/Additional File 4.xlsx", sheetName = "RS2", append = TRUE, row.names = T)
write.xlsx(RS3, file = "Tables/Additional File 4.xlsx", sheetName = "RS3", append = TRUE, row.names = T)
write.xlsx(RS6, file = "Tables/Additional File 4.xlsx", sheetName = "RS6", append = TRUE, row.names = T)

###DATA FOR FIGURE 3
########################################################################################################################
RMA.test.sig <- lapply(tests, toptable, f = f1, p = 0.01)
RMA.updown <- lapply(RMA.test.sig, function(dat) { if(dim(dat)[1] > 0) {c(length(which(dat[,3] > 0)), length(which(dat[,3] < 0))) } 
  else{c(0,0)}})
RMA.updown  <-  as.data.frame(matrix(unlist(RMA.updown), nrow = length(unlist(RMA.updown[1]))), row.names = c("RMA.Up","RMA.Down"))
colnames(RMA.updown) <- names(RMA.test.sig)

RUV.updown <- lapply(RUV.test.sig, function(dat) { c(length(which(dat[,3] > 0)), length(which(dat[,3] < 0))) })
RUV.updown  <-  as.data.frame(matrix(unlist(RUV.updown), nrow = length(unlist(RUV.updown[1]))), row.names = c("RUV.Up","RUV.Down"))
colnames(RUV.updown) <- names(RUV.test.sig)

updown <- t(rbind(RMA.updown, RUV.updown))
write.table(updown, "Tables/Figure 3 Data.txt", sep = "\t", quote = F, col.names = T, row.names = T)

###HEATMAP
##################################################################################################################
signames <- unique(unlist(lapply(RUV.test.sig, rownames)))
logratio <- s$normalizedCounts
colnames(logratio) <- x
logratio <- as.data.frame(avearrays(logratio))
logratio <- cbind(logratio$SD - logratio$CC6, logratio$RS1 - logratio$CC7, logratio$RS2 - logratio$CC8, logratio$RS3 - logratio$CC8, logratio$RS6 - logratio$CC11)
rownames(logratio) <- rownames(s$normalizedCounts) 
colnames(logratio) <- c("SD", "RS1", "RS2", "RS3", "RS6")
sigs <- logratio[rownames(logratio) %in% signames,]

breaks <- seq(-1.0, 1.0, by = 0.0001)
breaks <- append(breaks, 1.5)
breaks <- append(breaks, -1.5, after = 0)
mycol <- colorpanel(n = length(breaks) -1 , low = "green", mid = "black", high = "red")
c <- 7
color <- brewer.pal(c+1, "Set3")[c(1,3,4,8,6,7,5)]

# Rowv <- sigs %>% dist %>% hclust %>% as.dendrogram %>% set("branches_k_color", k = c, value = color)
# nRowv <- Rowv
# nRowv <- click_rotate(nRowv, horiz = T, continue = T)
load("ordered_dendrogram.rda")

hi <- cutree(nRowv, k = c)
labs <- nRowv %>% labels
col <- get_leaves_branches_col(nRowv)
col <- col[order(order.dendrogram(nRowv))]
cmap <- col[match(1:c, hi)]

svg(filename = "Figures/Figure3_heatmap.svg", width = 4, height = 4, pointsize = 8)
h <- heatmap.2(sigs, col = mycol, breaks = breaks, scale = "none", key = F, density.info = "none", 
               key.xlab = "log2(treatment/circadian control)", trace = "none", margins = c(8,1.5), 
               ylab = "", key.title = "", cexCol = 2.5, Colv = F, 
               dendrogram = "row", labRow = "", Rowv = nRowv, cex.lab = 1.5)
dev.off()

########################################################################################################################
#supplemental table 3
########################################################################################################################
cluster_data <- lapply(1:c, function(x) {as.data.frame(sigs[rownames(sigs) %in% rownames(sigs[hi == x,]),])})
addSymbol <- function(dat) {merge(dat, map_symbol, by.x = "row.names", by.y = "affy_mogene_2_1_st_v1", all.x = T)[,c(1,7,2:6)]}
cluster_data <- lapply(cluster_data, addSymbol)
addEnsembl <- function(dat) {merge(dat, map_ensembl, by.x = "Row.names", by.y = "affy_mogene_2_1_st_v1", all.x = T)[,c(1,8,2:7)]}
cluster_data <- lapply(cluster_data, addEnsembl)
names(cluster_data) <- paste0("Cluster", 1:c)
list2env(cluster_data, .GlobalEnv)

write.xlsx(Cluster1, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster1", row.names = F)
write.xlsx(Cluster2, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster2", append = TRUE, row.names = F)
write.xlsx(Cluster3, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster3", append = TRUE, row.names = F)
write.xlsx(Cluster4, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster4", append = TRUE, row.names = F)
write.xlsx(Cluster5, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster5", append = TRUE, row.names = F)
write.xlsx(Cluster6, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster6", append = TRUE, row.names = F)
write.xlsx(Cluster7, file = "Tables/Additional File 6.xlsx", sheetName = "Cluster7", append = TRUE, row.names = F)


###Time Series Plots
########################################################################################################################
mytheme <- theme(axis.text.y = element_text(size = 24),
                 axis.text.x = element_text(size = 24),
                 axis.title.y = element_text(size = 36, face = "bold"),
                 axis.title.x = element_text(size = 36, face = "bold"),
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill = NA, size = 2)) 


timeConstant <- function(cluster, genes, ki, a = 1) {
num <- cluster
cluster <- get(paste0("Cluster", cluster))
time <- c(0,1,2,3,6)
tmp <- cluster[,4:8]
id <- cluster[,3]
ind <- grep(genes, id)
n <- dim(tmp)[1]
if(mean(tmp[[1]]) < 0){
  tmp <- -tmp
}
m <- colMeans(tmp)
sd <- apply(tmp, 2, sd)
melted <- melt(tmp)
l <- dim(melted)[1]
melted[,1] <- rep(time, each = n)
colnames(melted) <- c("time", "logFC")
plot(logFC ~ time, data = melted)
abline(h = 0)
lines(time, m, col = cmap[num], lwd = 3)
i <- which(m == min(m))
fit <- nls(logFC ~ m[i] + (a - m[i]) * exp(-time/k), melted, start = c(a = a, k = ki)) 
lines(melted$time, predict(fit), col = cmap[num], lwd = 6)
k <- coef(fit)

list(melted, num, ind, m, sd, k)
}

timeSeriesUp <- function(dat1, dat2) {
  melted <- dat1[[1]]
  num <- dat1[[2]]
  ind <- dat1[[3]]
  m <- dat1[[4]]
  sd <- dat1[[5]]
  l <- dim(melted)[1]
  n <- l/5
  
  melted2 <- dat2[[1]]
  num2 <- dat2[[2]]
  ind2 <- dat2[[3]]
  m2 <- dat2[[4]]
  sd2 <- dat2[[5]]
  l2 <- dim(melted2)[1]
  n2 <- l2/5
  
  ps <- 7
  ls <- 1
  times <- c(0,1,2,3,6)
  
  p <- ggplot(melted[seq(ind[1],l,n),], aes(x = time, y = logFC))
  p + geom_point(color = cmap[num], size = ps) + geom_line(color = cmap[num], linetype = 2, size = ls) + 
    geom_line(aes(x = times, y = m), color = cmap[num], size = 3) +
    geom_point(data = melted[seq(ind[2],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[2],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_point(data = melted[seq(ind[3],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[3],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_ribbon(aes(ymin = m-sd, ymax = m+sd), alpha = 0.2, fill = cmap[num]) +
    geom_hline(aes(yintercept = 0), color = "grey", size = 2) +
    geom_line(aes(x = times, y = m2), color = cmap[num2], size = 3) +
    geom_point(data = melted2[seq(ind2[1],l2,n2),], aes(x = time, y = logFC), color = cmap[num2], size = ps) + 
    geom_line(data = melted2[seq(ind2[1],l2,n2),], aes(x = time, y = logFC), color = cmap[num2], linetype = 2, size = ls) +
    geom_point(data = melted2[seq(ind2[2],l2,n2),], aes(x = time, y = logFC), color = cmap[num2], size = ps) + 
    geom_line(data = melted2[seq(ind2[2],l2,n2),], aes(x = time, y = logFC), color = cmap[num2], linetype = 2, size = ls) +
    geom_ribbon(aes(ymin = m2-sd2, ymax = m2+sd2), alpha = 0.2, fill = cmap[num2]) +
    labs(x = "Treatment", y = "logFC") + 
    scale_x_continuous(breaks = times, labels=c("SD", "RS1", "RS2", "RS3", "RS6")) +
    scale_y_continuous(limits = c(-0.7, 1.5)) +
    mytheme
}

timeSeriesDown <- function(dat1) {
  dat1[[1]][,2] <- -dat1[[1]][,2]
  melted <- dat1[[1]]
  num <- dat1[[2]]
  ind <- dat1[[3]]
  if(length(ind) < 6){
  ind <- rep(ind, 6/length(ind))}
  m <- -dat1[[4]]
  sd <- -dat1[[5]]
  l <- dim(melted)[1]
  n <- l/5
  
  ps <- 7
  ls <- 1
  times <- c(0,1,2,3,6)
  
  p <- ggplot(melted[seq(ind[1],l,n),], aes(x = time, y = logFC))
  p + geom_point(color = cmap[num], size = ps) + geom_line(color = cmap[num], linetype = 2, size = ls) + 
    geom_line(aes(x = times, y = m), color = cmap[num], size = 3) +
    geom_point(data = melted[seq(ind[2],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[2],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_point(data = melted[seq(ind[3],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[3],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_point(data = melted[seq(ind[4],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[4],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_point(data = melted[seq(ind[5],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[5],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_point(data = melted[seq(ind[6],l,n),], aes(x = time, y = logFC), color = cmap[num], size = ps) + 
    geom_line(data = melted[seq(ind[6],l,n),], aes(x = time, y = logFC), color = cmap[num], linetype = 2, size = ls) +
    geom_ribbon(aes(ymin = m-sd, ymax = m+sd), alpha = 0.2, fill = cmap[num]) +
    geom_hline(aes(yintercept = 0), color = "grey", size = 2) +
    labs(x = "Treatment", y = "logFC") + 
    scale_x_continuous(breaks = times, labels=c("SD", "RS1", "RS2", "RS3", "RS6")) +
    scale_y_continuous(limits = c(-0.7, 1.5)) +
    mytheme
}

c4.data <- timeConstant(4, c("Arc|Per1|Per2"), 1)
c7.data <- timeConstant(7, c("Egr1|Egr2"), 1)
svg(filename = "Figures/Figure4_UpFast.svg", width = 8, height = 8, pointsize = 1)
timeSeriesUp(c4.data, c7.data)
dev.off()

c3.data <- timeConstant(3, c("Homer1|Bdnf|Fosb"), 1)
c1.data <- timeConstant(1, c("Hspa5|Npas2"), 1)
svg(filename = "Figures/Figure4_UpSlow.svg", width = 8, height = 8, pointsize = 1)
timeSeriesUp(c3.data, c1.data)
dev.off()

c2.data <- timeConstant(2, c("Eif3f|Hdac9|Usp2|Usp28|Sin3a|Nlgn1"), 1)
svg(filename = "Figures/Figure4_DownFast.svg", width = 8, height = 8, pointsize = 1)
timeSeriesDown(c2.data)
dev.off()

c5.data <- timeConstant(5, c("Cirbp|Dbp"), 4)
svg(filename = "Figures/Figure4_DownSlow.svg", width = 8, height = 8, pointsize = 1)
timeSeriesDown(c5.data)
dev.off()

###Circadian
########################################################################################################################
des2 <- model.matrix(~0 + x + s$W)
colnames(des2) <- c(names, paste0("W", 1:k))
f3 <- lmFit(dat, des2)
contrast.matrix <- makeContrasts(CC6-CC0, CC11-CC0, CC11-CC6, levels = des2)
f3 <- contrasts.fit(f3, contrast.matrix)
f3 <- eBayes(f3)

toptable <- function(x, f, p) { topTable(f, coef = x, adjust.method = "fdr", p.value = p, number = Inf, confint = F) }
tests <- c("CC6v0" = 1, "CC11v0" = 2, "CC11v6" = 3)
tests <- lapply(tests, toptable, f = f3, p = 1)
tests <- lapply(tests, addcolumn)
list2env(tests, .GlobalEnv)

write.xlsx(CC6v0, file = "Tables/Additional File 5 full.xlsx", sheetName = "CC6vCC0", row.names = F)
write.xlsx(CC11v0, file = "Tables/Additional File 5 full2.xlsx", sheetName = "CC11vCC0", row.names = F)#, append = TRUE)
write.xlsx(CC11v6, file = "Tables/Additional File 5 full3.xlsx", sheetName = "CC11vCC6", row.names = F)#, append = TRUE)
