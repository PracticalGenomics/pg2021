###############################################################################
####### THIS CODE IS FOR INFORMATIONAL PURPOSES , PLEASE START FROM  ##########
####### Day4_PM_Exercise.Rmd file                                    ##########
###############################################################################

library(scRNAseq)
library(org.Hs.eg.db)
library(AnnotationHub)

sce.grun <- GrunPancreasData()
gene.ids <- mapIds(org.Hs.eg.db, keys=rowData(sce.grun)$symbol,
                   keytype="SYMBOL", column="ENSEMBL")

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.grun <- sce.grun[keep,]
rownames(sce.grun) <- gene.ids[keep]
sce.grun <- sce.grun[,sce.grun$donor %in% c("D17", "D7", "D2")]
saveRDS(sce.grun, file="sce.grun.rds")

sce.muraro <- MuraroPancreasData()

#--- gene-annotation ---#
edb <- AnnotationHub()[["AH73881"]]
gene.symb <- sub("__chr.*$", "", rownames(sce.muraro))
gene.ids <- mapIds(edb, keys=gene.symb, 
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.muraro <- sce.muraro[keep,]
rownames(sce.muraro) <- gene.ids[keep]
sce.muraro <- sce.muraro[,sce.muraro$donor != "D28"]
saveRDS(sce.muraro, file="sce.muraro.rds")

sce.seger <- SegerstolpePancreasData()

#--- gene-annotation ---#
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

# Removing duplicated rows.
keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]

#--- sample-annotation ---#
emtab.meta <- colData(sce.seger)[,c("cell type", "disease",
                                    "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Disease", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
  toupper(substr(sce.seger$CellType, 1, 1)),
  substring(sce.seger$CellType, 2))
sce.seger <-  sce.seger[,!sce.seger$Donor %in% c("HP1504901", "HP1509101")]
sce.seger <- sce.seger[, sce.seger$Quality != "low quality cell"]
sce.seger <- sce.seger[,librarySizeFactors(altExp(sce.seger)) > 0 & sce.seger$Donor!="AZ"]
saveRDS(sce.seger, file="sce.seger.rds")