# created by Xioakang Liu on 21/06/2026 => xkliu@cug.edu.cn
library(geomorph)
library(dispRity)
library(factoextra)

df = read.csv("C:\\Users\\xiaokang\\Jupyter_code\\GastropodMorph&MMR\\gastropod_landmarks_twosets.csv",header=T,row.names = 1)
ammonoid.links <- as.matrix(read.table("C:\\Users\\xiaokang\\Jupyter_code\\GastropodMorph&MMR\\gastropod_links_combined.txt"))
curves<-as.matrix(read.csv("C:\\Users\\xiaokang\\Jupyter_code\\GastropodMorph&MMR\\gastropods_curveslide_combined.csv",header = T))
num_sp = length(row.names(df))
AA <- array(data = unlist(t(df[, 1:336])), 
            dim = c(2, 168, num_sp), 
            dimnames = list(NULL, NULL, rownames(df)))
BB <- aperm(AA,c(2,1,3))
X.gpa <- gpagen(BB, curves=curves,PrinAxes = FALSE,Proj=T)#PrinAxes = FALSE,Proj=T, approxBE=T)
#X.gpa <- gpagen(AA,PrinAxes = FALSE,ProcD = TRUE)
summary(X.gpa)
#plot(X.gpa)


pdf("twosets_Xgpa.pdf", width = 8, height = 6)
plot(X.gpa)
dev.off()

average.shape <- mshape(X.gpa$coords)#average shape points
PCA_data<-gm.prcomp(X.gpa$coords)
rownames(PCA_data$x) <- rownames(df)
pca_scores <- PCA_data$x[, 1:2]
pdf("twosets_PCA_Plot.pdf", width = 12, height = 10)
plot(pca_scores[,1], pca_scores[,2],
     #pch = '.',           
     #cex = 2,              
     xlab = "PC1", ylab = "PC2",
     main = "PCA of gastropods landmarks")
dev.off()

saveRDS(PCA_data, file = "./ModernGastropods/PCA_data.rds")