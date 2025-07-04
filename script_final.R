library(FactoMineR)
library(factoextra)

##Proteomica Día Corto

# Exploración de Datos
library(readr)
sd.swath.raw.data <- read_tsv("report.pg_matrix_SD.tsv")

nrow(sd.swath.raw.data)
ncol(sd.swath.raw.data)
head(sd.swath.raw.data)

# Renombrar las columnas de la tabla
colnames(sd.swath.raw.data) <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description","N.Sequences","N.Sequences",
                                    "ZT0_1", "ZT4_1", "ZT8_1", "ZT12_1", "ZT16_1", "ZT20_1",
                                    "ZT0_2", "ZT4_2", "ZT8_2", "ZT12_2", "ZT16_2", "ZT20_2",
                                    "ZT0_3", "ZT4_3", "ZT8_3", "ZT12_3", "ZT16_3", "ZT20_3")

# Verificar los nombres de las columnas
colnames(sd.swath.raw.data)

# Especificar el orden deseado de las columnas
ordered_columns_SD <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description","N.Sequences","N.Sequences",
                     "ZT0_1", "ZT0_2", "ZT0_3", "ZT4_1", "ZT4_2", "ZT4_3",
                     "ZT8_1", "ZT8_2", "ZT8_3", "ZT12_1", "ZT12_2", "ZT12_3",
                     "ZT16_1", "ZT16_2", "ZT16_3", "ZT20_1", "ZT20_2", "ZT20_3")

# Reorganizar el data frame según el orden especificado
sd.swath.raw.data.ordered <- sd.swath.raw.data[, ordered_columns_SD]

# Verificar el resultado
head(sd.swath.raw.data.ordered)

boxplot(sd.swath.raw.data.ordered[,7:24],outline=F,las=2,col=rep(rainbow(6),each=3),
        main="SD Proteomics Data Before Normalization",cex.main=2)


design <- data.frame(sample=colnames(sd.swath.raw.data.ordered[,7:24]),
                     group=c(rep("zt0",3),rep("zt4",3),rep("zt8",3),
                             rep("zt12",3),rep("zt16",3),rep("zt20",3)))

write.table(x = design,file = "SD_3days_design.tsv",quote = F,row.names = F,
            sep = "\t")


library("NormalyzerDE")

write.table(sd.swath.raw.data.ordered, file = "SD_data.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


normalyzer(jobName = "SD_normalization",designPath = "SD_3days_design.tsv",
           dataPath = "SD_data.tsv",outputDir = "./")

sd.normalized.proteomic.data <- read.table(
  file = "SD_normalization/Quantile-normalized.txt",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  fill = TRUE,
  stringsAsFactors = FALSE
)

head(sd.normalized.proteomic.data)

boxplot(sd.normalized.proteomic.data[,7:24],las=2,col=rep(rainbow(6),each=3),
        outline = F,main="SD Proteomics Data After Normalization",cex.main=2)

# Imputacion de Datos
sd.normalized.proteomic.data <- sd.normalized.proteomic.data[, -c(2:6)]
sum(is.na(sd.normalized.proteomic.data))
for(i in 1:nrow(sd.normalized.proteomic.data))
{
  if(sum(is.na(sd.normalized.proteomic.data[i,])) != 0)
  { 
    na.points <- colnames(sd.normalized.proteomic.data)[which(is.na(sd.normalized.proteomic.data[i,]))]
    for(j in 1:length(na.points))
    {
      zt <- strsplit(na.points[j],split="_")[[1]][1]  
      
      imputed.value <- mean(as.numeric(sd.normalized.proteomic.data[i,paste(zt,1:3,sep="_")]),na.rm = T)
      if(!is.nan(imputed.value))
      {
        sd.normalized.proteomic.data[i,na.points[j]] <- imputed.value
      }
    }
  }
}

sum(is.na(sd.normalized.proteomic.data))

for(i in 1:nrow(sd.normalized.proteomic.data))
{
  if(sum(is.na(sd.normalized.proteomic.data[i,])) != 0)
  {
    na.points <- colnames(sd.normalized.proteomic.data)[which(is.na(sd.normalized.proteomic.data[i,]))]
    
    zts <- unique(sapply(X = strsplit(na.points,split="_"), FUN = function(x){ return(x[1])}))
    
    for(j in 1:length(zts))
    {
      current.time.point <- as.numeric(strsplit(zts[j],split="T")[[1]][2])
      
      if(current.time.point > 0)
      {
        previous.time.point <- paste("ZT",current.time.point-4,sep="")
      } else
      {
        previous.time.point <- "ZT20"
      }
      
      if(current.time.point < 20)
      {
        next.time.point <- paste("ZT",current.time.point+4,sep="")
      } else
      {
        next.time.point <- "ZT0"
      }
      
      sd.normalized.proteomic.data[i,paste(zts[j],1:3,sep="_")] <-  mean(as.numeric(c(sd.normalized.proteomic.data[i,paste(previous.time.point,1:3,sep="_")],
                                                                                      sd.normalized.proteomic.data[i,paste(next.time.point,1:3,sep="_")])))
    }
  }
}

sum(is.na(sd.normalized.proteomic.data))

sd.normalized.proteomic.data[is.na(sd.normalized.proteomic.data)] <- min(sd.normalized.proteomic.data[,2:19],na.rm = T)
head(sd.normalized.proteomic.data)

sum(is.na(sd.normalized.proteomic.data))

# Matriz de abundancia
sd.protein.ids <- sd.normalized.proteomic.data$Protein.Group

sd.normalized.proteomic.data <- as.matrix(sd.normalized.proteomic.data[,2:19])
head(sd.normalized.proteomic.data)
rownames(sd.normalized.proteomic.data) <- sd.protein.ids

#sd.zt0 <- rowMeans(sd.normalized.proteomic.data[, c("ZT0_1", "ZT0_2", "ZT0_3")])
#sd.zt4 <- rowMeans(sd.normalized.proteomic.data[, c("ZT4_1", "ZT4_2", "ZT4_3")])
#sd.zt8 <- rowMeans(sd.normalized.proteomic.data[, c("ZT8_1", "ZT8_2", "ZT8_3")])
#sd.zt12 <- rowMeans(sd.normalized.proteomic.data[, c("ZT12_1", "ZT12_2", "ZT12_3")])
#sd.zt16 <- rowMeans(sd.normalized.proteomic.data[, c("ZT16_1", "ZT16_2", "ZT16_3")])
#sd.zt20 <- rowMeans(sd.normalized.proteomic.data[, c("ZT20_1", "ZT20_2", "ZT20_3")])


#sd.protein.abundance <- cbind(sd.zt0, sd.zt4, sd.zt8, sd.zt12, sd.zt16, sd.zt20)
#colnames(sd.protein.abundance) <- c("ZT0", "ZT4", "ZT8", "ZT12", "ZT16", "ZT20")
#rownames(sd.protein.abundance) <- sd.protein.ids

#Esto se haría si las replicas fuese tecnicas, pero solo son biológicas por lo que nos quedamos con la matriz inicial

sd.protein.abundance <- sd.normalized.proteomic.data
rownames(sd.protein.abundance) <- sd.protein.ids
colnames(sd.protein.abundance) <- c(paste("zt0",1:3,sep="_"),
                                    paste("zt4",1:3,sep="_"),
                                    paste("zt8",1:3,sep="_"),
                                    paste("zt12",1:3,sep="_"),
                                    paste("zt16",1:3,sep="_"),
                                    paste("zt20",1:3,sep="_"))

write.table(x = sd.protein.abundance,file = "SD_proteomic_data_abundance.tsv",quote = F,sep = "\t")

# Ritmicidad con RAIN
library("rain")
results.sd.proteomics <- rain(t(sd.protein.abundance), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.sd.proteomics$pVal < 0.05)/nrow(sd.protein.abundance)
rhythmic.proteins.sd <- rownames(subset(results.sd.proteomics, pVal < 0.05))

results.sd.proteomics.12 <- rain(t(sd.protein.abundance), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.sd.proteomics.12$pVal < 0.05)/nrow(sd.protein.abundance)
rhythmic.proteins.sd.12 <- rownames(subset(results.sd.proteomics.12, pVal < 0.05))

results.sd.proteomics.8 <- rain(t(sd.protein.abundance), deltat=4, period=8, verbose=FALSE, nr.series=3)
sum(results.sd.proteomics.8$pVal < 0.05)/nrow(sd.protein.abundance)

complete.rhythmic.proteins.sd <- unique(c(rhythmic.proteins.sd, rhythmic.proteins.sd.12))
length(complete.rhythmic.proteins.sd)
non.rhythmic.proteins.sd <- setdiff(sd.protein.ids, complete.rhythmic.proteins.sd)
length(non.rhythmic.proteins.sd)


data.proteins.sd <- matrix(c(length(complete.rhythmic.proteins.sd),
                             length(non.rhythmic.proteins.sd)),nrow=2)

sum(data.proteins.sd)

# Dibujar barplot apilado
bp <- barplot(data.proteins.sd, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
              col=c("red", "white"), 
              border="black", 
              font.axis=2,
              ylab = "Number of Proteins", 
              cex.lab = 1.6, 
              xlab="",ylim=c(0,3500)) #178 445

legend("topright",
       inset = c(-2.2, 0),  # mueve hacia fuera del margen derecho
       legend = c("Non rhythmic", "Rhythmic SD"),
       fill = c("white", "red"),
       border = "black",
       xpd = TRUE,           # permite dibujar fuera del área de la gráfica
       bty = "n",            # sin caja
       cex = 1)

write(x = complete.rhythmic.proteins.sd,file = "rhythmic_proteins_sd.txt")


#PCA
pca.protein.abundance.sd <- data.frame(colnames(sd.protein.abundance[complete.rhythmic.proteins.sd,]),t(sd.protein.abundance[complete.rhythmic.proteins.sd,]))
colnames(pca.protein.abundance.sd)[1] <- "Time point"

res.pca.protein.sd <- PCA(pca.protein.abundance.sd, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.protein.sd <- HCPC(res.pca.protein.sd, graph=FALSE)   
fviz_dend(res.hcpc.protein.sd,k=3,
          cex = 1,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 75      # Augment the room for labels
)

fviz_pca_ind(res.pca.protein.sd, col.ind = c("ZT00", "ZT00", "ZT00", "ZT04", "ZT04", "ZT04",
                                             "ZT08", "ZT08", "ZT08", "ZT12", "ZT12", "ZT12",
                                             "ZT16", "ZT16", "ZT16", "ZT20", "ZT20", "ZT20"),  
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE) 


##Proteomica Día Largo

# Exploración de Datos
library(readr)
ld.swath.raw.data <- read_tsv("report.pg_matrix_LD.tsv")

nrow(ld.swath.raw.data)
ncol(ld.swath.raw.data)
head(ld.swath.raw.data)

# Renombrar las columnas de la tabla
colnames(ld.swath.raw.data) <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description","N.Sequences","N.Sequences",
                                    "ZT0_1", "ZT4_1", "ZT8_1", "ZT12_1", "ZT16_1", "ZT20_1",
                                    "ZT0_2", "ZT4_2", "ZT8_2", "ZT12_2", "ZT16_2", "ZT20_2",
                                    "ZT0_3", "ZT4_3", "ZT8_3", "ZT12_3", "ZT16_3", "ZT20_3")

# Verificar los nombres de las columnas
colnames(ld.swath.raw.data)

# Especificar el orden deseado de las columnas
ordered_columns_LD <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description","N.Sequences","N.Sequences",
                     "ZT0_1", "ZT0_2", "ZT0_3", "ZT4_1", "ZT4_2", "ZT4_3",
                     "ZT8_1", "ZT8_2", "ZT8_3", "ZT12_1", "ZT12_2", "ZT12_3",
                     "ZT16_1", "ZT16_2", "ZT16_3", "ZT20_1", "ZT20_2", "ZT20_3")

# Reorganizar el data frame según el orden especificado
ld.swath.raw.data.ordered <- ld.swath.raw.data[, ordered_columns_LD]

# Verificar el resultado
head(ld.swath.raw.data.ordered)

boxplot(ld.swath.raw.data.ordered[,7:24],outline=F,las=2,col=rep(rainbow(6),each=3),
        main="LD Proteomics Data Before Normalization",cex.main=2)


design <- data.frame(sample=colnames(ld.swath.raw.data.ordered[,7:24]),
                     group=c(rep("zt0",3),rep("zt4",3),rep("zt8",3),
                             rep("zt12",3),rep("zt16",3),rep("zt20",3)))

write.table(x = design,file = "LD_3days_design.tsv",quote = F,row.names = F,
            sep = "\t")


library("NormalyzerDE")

write.table(ld.swath.raw.data.ordered, file = "LD_data.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


normalyzer(jobName = "LD_normalization",designPath = "LD_3days_design.tsv",
           dataPath = "LD_data.tsv",outputDir = "./")

ld.normalized.proteomic.data <- read.table(
  file = "LD_normalization/Quantile-normalized.txt",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  fill = TRUE,
  stringsAsFactors = FALSE
)
head(ld.normalized.proteomic.data)

boxplot(ld.normalized.proteomic.data[,7:24],las=2,col=rep(rainbow(6),each=3),
        outline = F,main="LD Proteomics Data After Normalization",cex.main=2)

# Imputacion de Datos
ld.normalized.proteomic.data <- ld.normalized.proteomic.data[, -c(2:6)]
sum(is.na(ld.normalized.proteomic.data))
for(i in 1:nrow(ld.normalized.proteomic.data))
{
  if(sum(is.na(ld.normalized.proteomic.data[i,])) != 0)
  { 
    na.points <- colnames(ld.normalized.proteomic.data)[which(is.na(ld.normalized.proteomic.data[i,]))]
    for(j in 1:length(na.points))
    {
      zt <- strsplit(na.points[j],split="_")[[1]][1]  
      
      imputed.value <- mean(as.numeric(ld.normalized.proteomic.data[i,paste(zt,1:3,sep="_")]),na.rm = T)
      if(!is.nan(imputed.value))
      {
        ld.normalized.proteomic.data[i,na.points[j]] <- imputed.value
      }
    }
  }
}

sum(is.na(ld.normalized.proteomic.data))

for(i in 1:nrow(ld.normalized.proteomic.data))
{
  if(sum(is.na(ld.normalized.proteomic.data[i,])) != 0)
  {
    na.points <- colnames(ld.normalized.proteomic.data)[which(is.na(ld.normalized.proteomic.data[i,]))]
    
    zts <- unique(sapply(X = strsplit(na.points,split="_"), FUN = function(x){ return(x[1])}))
    
    for(j in 1:length(zts))
    {
      current.time.point <- as.numeric(strsplit(zts[j],split="T")[[1]][2])
      
      if(current.time.point > 0)
      {
        previous.time.point <- paste("ZT",current.time.point-4,sep="")
      } else
      {
        previous.time.point <- "ZT20"
      }
      
      if(current.time.point < 20)
      {
        next.time.point <- paste("ZT",current.time.point+4,sep="")
      } else
      {
        next.time.point <- "ZT0"
      }
      
      ld.normalized.proteomic.data[i,paste(zts[j],1:3,sep="_")] <-  mean(as.numeric(c(ld.normalized.proteomic.data[i,paste(previous.time.point,1:3,sep="_")],
                                                                                      ld.normalized.proteomic.data[i,paste(next.time.point,1:3,sep="_")])))
    }
  }
}

sum(is.na(ld.normalized.proteomic.data))

ld.normalized.proteomic.data[is.na(ld.normalized.proteomic.data)] <- min(ld.normalized.proteomic.data[,2:19],na.rm = T)
head(ld.normalized.proteomic.data)

sum(is.na(ld.normalized.proteomic.data))

# Matriz de abundancia
ld.protein.ids <- ld.normalized.proteomic.data$Protein.Group

ld.normalized.proteomic.data <- as.matrix(ld.normalized.proteomic.data[,2:19])
head(ld.normalized.proteomic.data)
rownames(ld.normalized.proteomic.data) <- ld.protein.ids

#ld.zt0 <- rowMeans(ld.normalized.proteomic.data[, c("ZT0_1", "ZT0_2", "ZT0_3")])
#ld.zt4 <- rowMeans(ld.normalized.proteomic.data[, c("ZT4_1", "ZT4_2", "ZT4_3")])
#ld.zt8 <- rowMeans(ld.normalized.proteomic.data[, c("ZT8_1", "ZT8_2", "ZT8_3")])
#ld.zt12 <- rowMeans(ld.normalized.proteomic.data[, c("ZT12_1", "ZT12_2", "ZT12_3")])
#ld.zt16 <- rowMeans(ld.normalized.proteomic.data[, c("ZT16_1", "ZT16_2", "ZT16_3")])
#ld.zt20 <- rowMeans(ld.normalized.proteomic.data[, c("ZT20_1", "ZT20_2", "ZT20_3")])


#ld.protein.abundance <- cbind(ld.zt0, ld.zt4, ld.zt8, ld.zt12, ld.zt16, ld.zt20)
#colnames(ld.protein.abundance) <- c("ZT0", "ZT4", "ZT8", "ZT12", "ZT16", "ZT20")
#rownames(ld.protein.abundance) <- ld.protein.ids

#Esto se haría si las replicas fuese tecnicas, pero solo son biológicas por lo que nos quedamos con la matriz inicial

ld.protein.abundance <- ld.normalized.proteomic.data
rownames(ld.protein.abundance) <- ld.protein.ids
colnames(ld.protein.abundance) <- c(paste("zt0",1:3,sep="_"),
                                    paste("zt4",1:3,sep="_"),
                                    paste("zt8",1:3,sep="_"),
                                    paste("zt12",1:3,sep="_"),
                                    paste("zt16",1:3,sep="_"),
                                    paste("zt20",1:3,sep="_"))

write.table(x = ld.protein.abundance,file = "LD_proteomic_data_abundance.tsv",quote = F,sep = "\t")

# Ritmicidad con RAIN
library("rain")
results.ld.proteomics <- rain(t(ld.protein.abundance), deltat=4, period=24, verbose=FALSE, nr.series=3)
sum(results.ld.proteomics$pVal < 0.05)/nrow(ld.protein.abundance)
rhythmic.proteins.ld <- rownames(subset(results.ld.proteomics, pVal < 0.05))

results.ld.proteomics.12 <- rain(t(ld.protein.abundance), deltat=4, period=12, verbose=FALSE, nr.series=3)
sum(results.ld.proteomics.12$pVal < 0.05)/nrow(ld.protein.abundance)
rhythmic.proteins.ld.12 <- rownames(subset(results.ld.proteomics.12, pVal < 0.05))

results.ld.proteomics.8 <- rain(t(ld.protein.abundance), deltat=4, period=8, verbose=FALSE, nr.series=3)
sum(results.ld.proteomics.8$pVal < 0.05)/nrow(ld.protein.abundance)

complete.rhythmic.proteins.ld <- unique(c(rhythmic.proteins.ld, rhythmic.proteins.ld.12))
length(complete.rhythmic.proteins.ld)
non.rhythmic.proteins.ld <- setdiff(ld.protein.ids, complete.rhythmic.proteins.ld)
length(non.rhythmic.proteins.ld)


data.proteins.ld <- matrix(c(length(complete.rhythmic.proteins.ld),
                             length(non.rhythmic.proteins.ld)),nrow=2)

sum(data.proteins.ld)


# Configurar márgenes para permitir espacio a la derecha
par(mar = c(5, 6, 4, 7))  # abajo, izq, arriba, der

# Dibujar barplot apilado
bp <- barplot(data.proteins.ld, lwd=3, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue", "white"), 
        border="black", 
        font.axis=2,
        ylab = "Number of Proteins", 
        cex.lab = 1.6, 
        xlab="",ylim=c(0,3500)) #178 445

legend("topright",
       inset = c(-2.2, 0),  # mueve hacia fuera del margen derecho
       legend = c("Non rhythmic", "Rhythmic LD"),
       fill = c("white", "blue"),
       border = "black",
       xpd = TRUE,           # permite dibujar fuera del área de la gráfica
       bty = "n",            # sin caja
       cex = 1)


write(x = complete.rhythmic.proteins.ld,file = "rhythmic_proteins_ld.txt")


#PCA
pca.protein.abundance.ld <- data.frame(colnames(ld.protein.abundance[complete.rhythmic.proteins.ld,]),t(ld.protein.abundance[complete.rhythmic.proteins.ld,]))
colnames(pca.protein.abundance.ld)[1] <- "Time point"

res.pca.protein.ld <- PCA(pca.protein.abundance.ld, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc.protein.ld <- HCPC(res.pca.protein.ld, graph=FALSE)   
fviz_dend(res.hcpc.protein.ld,k=3,
          cex = 1,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 75      # Augment the room for labels
)

fviz_pca_ind(res.pca.protein.ld, col.ind =  c("ZT00", "ZT00", "ZT00", "ZT04", "ZT04", "ZT04",
                                              "ZT08", "ZT08", "ZT08", "ZT12", "ZT12", "ZT12",
                                              "ZT16", "ZT16", "ZT16", "ZT20", "ZT20", "ZT20"), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE) 


# Diagrama de Venn: SD vs LD
library(grid)
library(VennDiagram)

grid.newpage()
draw.pairwise.venn(area1 = length(complete.rhythmic.proteins.ld),area2 = length(complete.rhythmic.proteins.sd),cross.area =  length(intersect(complete.rhythmic.proteins.ld,complete.rhythmic.proteins.sd)),lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c(" blue","red"),alpha = 0.3,cex = 2,cat.cex = 2)

