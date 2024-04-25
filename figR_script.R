ROOT_DIR<-"/fh/fast/furlan_s/grp/experiments"
stem<-"pilot_cd27"
DATA_DIR <- file.path(ROOT_DIR,  stem, "data")      # SPECIFY HERE
RES_DIR  <- file.path(ROOT_DIR,  stem, "res")     # SPECIFY HERE
RMD_DIR  <- file.path(ROOT_DIR,  stem, "rmd")     # SPECIFY HERE
CDS_DIR <- file.path(ROOT_DIR,  stem, "cds")
FIG_DIR <- file.path(ROOT_DIR,  stem, "figs")

suppressPackageStartupMessages({
  library(monocle3)
  library(reticulate)
  library(openxlsx)  
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(xfun)
  library(pals)
  library(RColorBrewer)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(future)
  library(GenomicRanges)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  library(stringr)
  library(ggsignif)
  library(ggpubr)
  library(ggplot2)
  library(ggrastr)
  library(BuenColors)
  library(FigR)
  library(ggrepel)
  library(fgsea)
  library(ArchR)
  library(parallel)
  library(scCustomize)
  library(pbmcapply)
  library(FigR)
})

set.seed(1234)

MO<-readRDS( file.path(CDS_DIR, "MO_pilot_cd27.rds"))

ATAC.se <- SummarizedExperiment(assays = list(counts = MO[["ATAC"]]@data) , rowRanges =MO[["ATAC"]]@ranges, colData = MO@meta.data)

RNAmat<- MO[["SCT"]]@data %>% as.sparse()

cisCor <- runGenePeakcorr(ATAC.se = ATAC.se,
                           RNAmat = RNAmat,
                           genome = "hg38", # Also supports mm10 and hg38
                           nCores = detectCores(), 
                           p.cut=NULL)

write.csv(cisCor, file.path(RES_DIR, "cisCor.csv"))