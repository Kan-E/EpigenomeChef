##with RNAseq homer startボタン
##Motif plot order反映されない
##rmdの修正
library(rtracklayer) 
library(Rsubread)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(TxDb.Btaurus.UCSC.bosTau8.refGene)
library(TxDb.Cfamiliaris.UCSC.canFam3.refGene)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(TxDb.Ggallus.UCSC.galGal4.refGene)
library(TxDb.Mmulatta.UCSC.rheMac8.refGene)
library(TxDb.Ptroglodytes.UCSC.panTro4.refGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)
library(org.Bt.eg.db)
library(org.Cf.eg.db)
library(org.Dr.eg.db)
library(org.Gg.eg.db)
library(org.Mmu.eg.db)
library(org.Pt.eg.db)
library(biomaRt)
library(shiny)
library(DT)
library(rstatix)
library(multcomp)
library(tidyverse)
library(tools)
library(ggpubr)
library(ggrepel)
library(ggdendro)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(ggnewscale)
library(AnnotationDbi)
library(msigdbr)
library(genefilter)
library(ComplexHeatmap)
library(shinyBS, verbose = FALSE)
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library(BiocManager)
library(clusterProfiler.dplyr)
library(GenomicRanges)
library(rGREAT)
library(FindIT2)
library(limma)
library(ggseqlogo)
library(colorspace)
library(ggcorrplot)
library(RColorBrewer)

library(venn)
library(reshape2)
library(ggsci)
library(ggrastr) ##devtools::install_github('VPetukhov/ggrastr')
library(EnrichedHeatmap)
library(pdftools)
library(magick)
library(webshot)
library(clue)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(stringr)
library(readxl)
file.copy("pair_report.Rmd",file.path(tempdir(),"pair_report.Rmd"), overwrite = TRUE)
home_dir <- paste0(getwd(),"/")
Sys.setenv(OPENSSL_CONF="/dev/null")
jscode <- "shinyjs.closeWindow = function() { window.close(); }"
options(rsconnect.max.bundle.size=31457280000000000000)
species_list <- c("not selected", "Homo sapiens (hg19)","Homo sapiens (hg38)","Mus musculus (mm10)","Mus musculus (mm39)",
                  "Rattus norvegicus (rn6)","Drosophila melanogaster (dm6)","Caenorhabditis elegans (ce11)",
                  "Bos taurus (bosTau8)","Canis lupus familiaris (canFam3)","Danio rerio (danRer10)",
                  "Gallus gallus (galGal4)","Macaca mulatta (rheMac8)","Pan troglodytes (panTro4)")
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target","Position")
gene_set_list_genome <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "Transcription factor targets", "miRNA target")
org <- function(Species){
  if(Species != "not selected"){
    switch (Species,
            "Homo sapiens (hg38)" = org <- org.Hs.eg.db,
            "Homo sapiens (hg19)" = org <- org.Hs.eg.db,
            "Mus musculus (mm10)" = org <- org.Mm.eg.db,
            "Mus musculus (mm39)" = org <- org.Mm.eg.db,
            "Drosophila melanogaster (dm6)" = org <- org.Dm.eg.db,
            "Rattus norvegicus (rn6)" = org <- org.Rn.eg.db,
            "Caenorhabditis elegans (ce11)" = org <- org.Ce.eg.db,
            "Bos taurus (bosTau8)" = org <- org.Bt.eg.db,
            "Canis lupus familiaris (canFam3)" = org <- org.Cf.eg.db,
            "Danio rerio (danRer10)" = org <- org.Dr.eg.db,
            "Gallus gallus (galGal4)" = org <- org.Gg.eg.db,
            "Macaca mulatta (rheMac8)" = org <- org.Mmu.eg.db,
            "Pan troglodytes (panTro4)" = org <- org.Pt.eg.db,
            "Saccharomyces cerevisiae (sacCer3)" = org <- org.Sc.sgd.db,
            "Xenopus laevis (xenLae2)" = org <- org.Xl.eg.db,
            "Arabidopsis thaliana (tair10)" = org <- org.At.tair.db
            )
    return(org)
  }
}
bed_name <- function(name){
  name <- gsub(".+\\/","",name)
  name <- gsub("\\.bed$", "", name)
  name <- gsub("\\.narrowPeak$", "", name)
  name <- gsub("\\.narrowpeak$", "", name)
  name <- gsub("\\.sumit$", "", name)
  name <- gsub("\\.\\.\\.\\."," < ",name)
  return(name)
}
bigwig_name <- function(name){
  name <- gsub(".+\\/","",name)
  name <- gsub("\\.bigwig$", "", name)
  name <- gsub("\\.BigWig$", "", name)
  name <- gsub("\\.bw$", "", name)
  name <- gsub("\\.Bigwig$", "", name)
  name <- gsub("\\.\\.\\.\\."," < ",name)
  return(name)
}
bam_name <- function(name){
  name <- gsub(".+\\/","",name)
  name <- gsub("\\.bam$", "", name)
  name <- gsub("\\.Bam$", "", name)
  name <- gsub("\\.\\.\\.\\."," < ",name)
  return(name)
}
  
bws_ordering <- function(bws,sample_order,additional=TRUE){
  if(!is.null(bws) && !is.null(sample_order)){
    order <- sample_order
    bws <- bws
    order <- gsub("\\.","-",order)
    bws <- bws[order]
    return(bws)
  }else if(additional == FALSE)  validate("Bigwig files are required.")
}
read_df <- function(tmp, Species=NULL){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(tools::file_ext(tmp) == "xlsx") {
      df2 <- readxl::read_xlsx(tmp) 
      df2 <- as.data.frame(df2)
      df <- try(data.frame(row.names = df2[,1]),silent = T)
      if(class(df) != "try-error") {
        if(dim(df2)[2] == 2){
          df <- data.frame(row.names = df2[,1],a = df2[,2])
          colnames(df)[1] <- colnames(df2)[2]
        }else{
          rownames(df2) <- df2[,1]
          df <- df2[,-1]
          colnames(df) <- gsub("-",".",colnames(df))
        }
      }
    }
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
    if(tools::file_ext(tmp) == "txt" || tools::file_ext(tmp) == "tsv") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
    rownames(df) = gsub("\"", "", rownames(df))
    rownames(df) = gsub(":", ".", rownames(df))
    rownames(df) = gsub("\\\\", ".", rownames(df))
    if(length(grep("SYMBOL", colnames(df))) != 0){
      df <- df[, - which(colnames(df) == "SYMBOL")]
    }
    if(length(colnames(df)) != 0){
      if(str_detect(colnames(df)[1], "^X\\.")){
        colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
      }
    }
    return(df)
  }
}
read_dfs <- function(tmp, Volume=NULL,Species=NULL){
  name = c()
  upload = list()
  if(!is.null(Volume)){
    for(nr in 1:length(tmp)){
      df <- read_df(tmp[[nr]])
      file_name <- gsub(".+\\/","",tmp[[nr]])
      name <- c(name, gsub("\\..+$", "", file_name))
      upload[[nr]] <- df
    }
    names(upload) <- name
  }else{
  for(nr in 1:length(tmp[, 1])){
    df <- read_df(tmp[[nr, 'datapath']])
    file_name <- gsub(paste0("\\.",tools::file_ext(tmp[[nr, 'datapath']]),"$"), "", tmp[nr,]$name)
    name <- c(name, file_name)
    upload[[nr]] <- df
  }
  names(upload) <- name
  }
  return(upload)
}
txdb_function <- function(Species){
  if(Species == "Xenopus laevis (xenLae2)"){
    xenLae2<-makeTxDbFromUCSC(genome="xenLae2", tablename="ncbiRefSeq",
                              transcript_ids=NULL,
                              circ_seqs=NULL,
                              url="http://genome.ucsc.edu/cgi-bin/",
                              goldenPath.url=getOption("UCSC.goldenPath.url"),
                              taxonomyId=NA,
                              miRBaseBuild=NA)
  }
  switch (Species,
          "Mus musculus (mm10)" = txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene,
          "Homo sapiens (hg19)" = txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene,
          "Homo sapiens (hg38)" = txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene,
          "Mus musculus (mm39)" = txdb <- TxDb.Mmusculus.UCSC.mm39.refGene,
          "Drosophila melanogaster (dm6)" = txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene,
          "Rattus norvegicus (rn6)" = txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene,
          "Caenorhabditis elegans (ce11)" = txdb <- TxDb.Celegans.UCSC.ce11.refGene,
          "Bos taurus (bosTau8)" = txdb <- TxDb.Btaurus.UCSC.bosTau8.refGene,
          "Canis lupus familiaris (canFam3)" = txdb <- TxDb.Cfamiliaris.UCSC.canFam3.refGene,
          "Danio rerio (danRer10)" = txdb <- TxDb.Drerio.UCSC.danRer10.refGene,
          "Gallus gallus (galGal4)" = txdb <- TxDb.Ggallus.UCSC.galGal4.refGene,
          "Macaca mulatta (rheMac8)" = txdb <- TxDb.Mmulatta.UCSC.rheMac8.refGene,
          "Pan troglodytes (panTro4)" = txdb <- TxDb.Ptroglodytes.UCSC.panTro4.refGene,
          "Saccharomyces cerevisiae (sacCer3)" = txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,
          "Xenopus laevis (xenLae2)" = txdb <- xenLae2,
          "Arabidopsis thaliana (tair10)" = txdb <- TxDb.Athaliana.BioMart.plantsmart51)
  return(txdb)
}
files_name2ENTREZID <- function(files,Species,gene_type){
  df2 <- list()
  if(!is.null(gene_type)){
  for (name in names(files)) {
    count <- files[[name]]
    if(!is.na(suppressWarnings(as.numeric(rownames(count)[1]))) == TRUE){
      df2[[name]] <- count
    }else{
      if(gene_type[name] != "SYMBOL"){
        if(str_detect(rownames(count)[1], "FBgn")){
          count <- as.data.frame(count)
          count$ENTREZID <- rownames(count)
          count2 <- count %>% dplyr::select(ENTREZID, everything())
        }else{
        my.symbols <- gsub("\\..*","", rownames(count))
        gene_id <-id_convert(my.symbols,Species = Species,type = "ENSEMBL2ENTREZID")
        gene_id <- gene_id %>% distinct(ENSEMBL, .keep_all = T) %>% na.omit()
        gene_id <- data.frame(ENTREZID = gene_id$ENTREZID, row.names = gene_id$ENSEMBL)
        count2 <- merge(gene_id,count,by=0)
        count2 <- count2[,-1]
        count2 <- count2 %>% distinct(ENTREZID, .keep_all = T)
        }
      }else{
        #symbol
        count2 <- symbol2gene_id(data = count,org=org(Species))
        colnames(count2)[1]<-"ENTREZID"
      }
      df2[[name]] <- count2
    }
  }
  return(df2)
  }
}
filter_function <- function(filter, files){
  if(filter == "Reproducible_peaks"){
    if(length(names(files)) == 1) {
      consensus <- list(files[[1]],files[[1]])
    }else{
    name_list <- names(files) %>% sort()
    files <- files[name_list]
    unique_col <- unique(gsub("\\_.+$", "", name_list))
    total <- 0
    consensus <- list()
    for(i in 1:length(unique_col)){
      cond <- length(which(gsub("\\_.+$", "", name_list) == unique_col[i]))
      files2 <- files[(1 + total):(total + cond)]
      print(cond)
      if(cond == 1) files2 <- list(files2[[1]],files2[[1]])
        consensusToCount <- soGGi:::runConsensusRegions(GRangesList(files2), "none")
        occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
          rowSums
        consensusToCount <- consensusToCount[occurrences >= 2, ]
        consensus[[i]] <- consensusToCount
      total <- total + cond
      print(consensus)
    }
    if(length(unique(gsub("\\_.+$", "", names(files)))) == 1) consensus <- list(consensus[[1]],consensus[[1]])
    }
    consensusToCount <- soGGi:::runConsensusRegions(GRangesList(consensus), "none")
  }else{
    if(length(names(files)) == 1) files <- list(files[[1]],files[[1]])
      consensusToCount <- soGGi:::runConsensusRegions(GRangesList(files), "none")
  }
  return(consensusToCount)
}
promoter <- function(txdb, upstream, downstream,input_type = "Promoter",files = NULL, bam=F,RPM,filter="Reproducible_peaks"){
  if(input_type == "Promoter"){
  return(promoters(genes(txdb),upstream = upstream, downstream = downstream))
  }else{
   return(filter_function(files = files,filter=filter))
  }
}
promoter_clustering <- function(txdb, upstream, downstream,input_type = "Promoter",files = NULL, bam=F,RPM, filter="Reproducible_peaks"){
  if(input_type == "Promoter"){
    return(promoters(genes(txdb),upstream = upstream, downstream = downstream))
  }else{
    return(filter_function(files = files,filter=filter))
  }
}

Bigwig2count <- function(bw, promoter, Species, input_type = "Promoter"){
bed1<-promoter
write.table(bed1,file = paste0(tempdir(),"bed.bed"),sep = "\t")
bed1 <- read.table(paste0(tempdir(),"bed.bed"),header = T)
bed1 <- dplyr::arrange(bed1, seqnames)
data <- as.data.frame(bed1)[,1:3]
colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start,end)))
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))
perc <- 0
withProgress(message = "converting BigWig to gene count data",{
for(i in seq_len(length(bw))) {
    perc <- perc + 1
  last=1 	 
  coverage <- try(import(bw[[i]], as = 'RleList'))
  if(class(coverage) == "try-error") {
    validate(paste0("Error: the uploaded bigwig files are in an unexpected format. The original error message is as follows:\n",print(coverage)))
  }
  for(j in chromnames){
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage)){
      counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }else{
      newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
      counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
  }
  incProgress(1/length(colnames(counts)), message = paste("Finish '", colnames(counts)[i], "', ", 
                                            perc, "/", length(colnames(counts)),sep = ""))
}
})
if(input_type == "Promoter"){
  or <- org(Species)
  rownames(counts) <- bed1$gene_id
  if(!str_detect(rownames(counts)[1], "FBgn")){
my.symbols <- rownames(counts)
gene_IDs<-AnnotationDbi::select(or,keys = my.symbols,
                                keytype = "ENTREZID",
                                columns = c("ENTREZID","SYMBOL"))
colnames(gene_IDs) <- c("Row.names","SYMBOL")
gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
gene_IDs <- data.frame(SYMBOL = gene_IDs$SYMBOL, row.names = gene_IDs$Row.names)
data2 <- merge(gene_IDs,counts, by=0)
rownames(data2) <- data2$SYMBOL
counts <- data2[,-1:-2]
}
}else{
  a <- as.data.frame(bed)
  Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
  rownames(counts) <- Row.name
}
return(counts)
}



PCAplot <- function(data, plot,legend=NULL){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
  }
  pca <- prcomp(data, scale. = T)
  label<- colnames(data)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$rotation <- as.data.frame(pca$rotation)
  if(plot==TRUE){
    if(!is.null(legend)) {
      if(legend == "Legend"){
        legend_position <- "top" 
        label2<- NULL
      }else{
        legend_position <- "none"
        label2<- label
      }
    }
    g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                  y=pca$rotation[,2],
                                  col=gsub("\\_.+$", "", label), label = label2)) +
      geom_point()+
      theme(panel.background =element_rect(fill=NA,color=NA),
            panel.border = element_rect(fill = NA)) +
      xlab(lab_x) + ylab(lab_y) +
      theme(legend.position=legend_position, aspect.ratio=1)+ 
      guides(color=guide_legend(title=""))
    if(!is.null(legend)){
      if(legend == "Label") g1 <- g1 + geom_text_repel(show.legend = NULL)
    }
  rho <- cor(data,method="spearman")
  d <- dist(1-rho)
  mds <- as.data.frame(cmdscale(d))
  label<-colnames(data)
  label<-gsub("\\_.+$", "", label)
  g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                        col = gsub("\\_.+$", "", label), label = label2)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab("dim 1") + ylab("dim 2") +
    theme(legend.position=legend_position, aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  if(!is.null(legend)){
    if(legend == "Label") g2 <- g2 + geom_text_repel(show.legend = NULL)
  }
  x <- NULL
  y <- NULL
  xend <- NULL
  yend <- NULL
  data.t <- t(data)
  hc <- hclust(dist(data.t), "ward.D2")
  dendr <- dendro_data(hc, type="rectangle")
  g3 <- ggplot() +
    geom_segment(data=segment(dendr),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr),
              aes(x, y, label=label, hjust=0), size=3) +
    theme(legend.position = "none",
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),axis.text.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title.y=element_blank(),panel.background=element_rect(fill="white"))+
    coord_flip()+ scale_y_reverse(expand=c(0.6, 0))
  p2 <- plot_grid(g1, g2, g3, nrow = 1)
  return(p2)
  }else return(pca$rotation)
}
umap_plot <- function(data, n_neighbors){
  umap <- umap::umap(t(data),n_neighbors = n_neighbors, random_state = 123)
  data2 <- umap$layout %>% as.data.frame()
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  p<- ggplot(data2, mapping = aes(V1,V2, color = label, label = colnames(data)))+
    geom_point()+geom_text_repel()+ xlab("UMAP_1") + ylab("UMAP_2")+
    theme(panel.background =element_rect(fill=NA,color=NA),panel.border = element_rect(fill = NA),
          aspect.ratio=1)
  return(p)
}

GOIheatmap <- function(data.z, show_row_names = TRUE){
  ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                clustering_method_columns = 'ward.D2',
                show_row_names = show_row_names, show_row_dend = F,column_names_side = "top",
                row_names_gp = gpar(fontface = "italic"),use_raster = TRUE)
  return(ht)
}
pdf_h <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_hsize <- 44
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_hsize <- 40
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_hsize <- 36
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_hsize <- 25
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_hsize <- 25
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_hsize <- 20
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_hsize <- 15
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_hsize <- 15
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_hsize <- 15
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_hsize <- 10
  if (length(rowlist) == 4) pdf_hsize <- 10
  if (length(rowlist) == 3) pdf_hsize <- 5
  if (length(rowlist) == 2) pdf_hsize <- 5
  if (length(rowlist) == 1) pdf_hsize <- 5
  if (length(rowlist) > 200) pdf_hsize <- 30
  return(pdf_hsize)
}
pdf_w <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_wsize <- 44
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_wsize <- 40
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_wsize <- 36
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_wsize <- 30
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_wsize <- 25
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_wsize <- 20
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_wsize <- 17
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_wsize <- 17
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_wsize <- 14
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_wsize <- 17
  if (length(rowlist) == 4) pdf_wsize <- 14
  if (length(rowlist) == 3) pdf_wsize <- 17
  if (length(rowlist) == 2) pdf_wsize <- 14
  if (length(rowlist) == 1) pdf_wsize <- 8
  if (length(rowlist) > 200) pdf_wsize <- 30
  return(pdf_wsize)
}
pdfSize_for_GOI <- paste(strong("Heatmap:"),"height = 10, width = 7 <br>", 
                         strong("Boxplot:"),"<br>",
                         "Gene number = 1,","height = 3, width = 3 <br>",
                         "Gene number = 2,","height = 3, width = 6 <br>",
                         "Gene number = 3,","height = 3, width = 9 <br>",
                         "Gene number = 4,","height = 6, width = 6 <br>",
                         "Gene number = 5 ~ 6,","height = 6, width = 9 <br>",
                         "Gene number = 7 ~ 9,","height = 7.5, width = 6.75 <br>",
                         "Gene number = 10 ~ 12,","height = 7.5, width = 9 <br>",
                         "Gene number = 13 ~ 16,","height = 9, width = 9 <br>",
                         "Gene number = 17 ~ 25,","height = 11.5, width = 11.5 <br>",
                         "Gene number = 26 ~ 36,","height = 13.5, width = 13.5 <br>",
                         "Gene number = 37 ~ 49,","height = 15.75, width = 15.75 <br>",
                         "Gene number = 50 ~ 64,","height = 18, width = 18 <br>",
                         "Gene number = 65 ~ 81,","height = 20.5, width = 20.5 <br>",
                         "Gene number = 82 ~ 200,","height = 22.5, width = 22.5 <br>",
                         "Gene number > 200,", "height = 30, width = 30 <br>")

dotplot_for_output <- function(data, plot_genelist, Gene_set, Species){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      withProgress(message = "Plot results",{
        if(Species != "not selected"){
          print(plot_genelist)
        }
        incProgress(1)
      })
    }
  }
}
GeneList_for_enrichment <- function(Species, Gene_set, org, Custom_gene_list){
  if(Species != "not selected" && !is.null(Gene_set) && !is.null(org)){
    species <- substr(gsub("\\(.+$","",Species),1,nchar(gsub("\\(.+$","",Species))-1)
    if(Gene_set == "MSigDB Hallmark"){
      H_t2g <- msigdbr(species = species, category = "H") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description,ensembl_gene) 
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
    }
    if(Gene_set == "Position"){
      H_t2g <- msigdbr(species = species, category = "C1") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
    }
    if(Gene_set == "KEGG"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="KEGG_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Transcription factor targets"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "DoRothEA regulon (activator)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (activator)")) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "DoRothEA regulon (repressor)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (repressor)")) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "Reactome"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="REACTOME_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "miRNA target"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "MIR:MIRDB" | gs_subcat == "MIR:MIR_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "GO biological process"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOBP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO cellular component"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:CC") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOCC_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO molecular function"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:MF") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOMF_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Human phenotype ontology"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "HPO") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "WikiPathways"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="WP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "PID (Pathway Interaction Database)"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:PID") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PID_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "BioCarta"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:BIOCARTA") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="BIOCARTA_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Custom gene set"){
      if(!is.null(Custom_gene_list)){
        H_t2g <- gene_list_convert_for_enrichment(data= Custom_gene_list, Species = Species)
        H_t2g <- data.frame(gs_name = H_t2g$Group, entrez_gene = H_t2g$ENTREZID)
        H_t2g$gs_name <- gsub(":", "_", H_t2g$gs_name)
      }else H_t2g <- NULL
    }
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tnf", replacement = "TNF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tgf", replacement = "TGF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Pi3k_akt_mtor", replacement = "PI3K_Akt_mTOR")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il6_", replacement = "IL6_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il2_", replacement = "IL2_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Kras", replacement = "KRas")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Uv_r", replacement = "UV_r")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Dna_", replacement = "DNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Rna_", replacement = "RNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mrna_", replacement = "mRNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="E2f", replacement = "E2F")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="G2m", replacement = "G2M")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mtorc", replacement = "mTORC")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Ecm_", replacement = "ECM_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Abc_", replacement = "ABC_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="No1_", replacement = "NO1_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_mirna", replacement = "_miRNA")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="gtpase_", replacement = "GTPase_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="rho_", replacement = "Rho_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="jnk_", replacement = "JNK_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_jak", replacement = "_JAK")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_stat", replacement = "_STAT")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_nfkb", replacement = "_NFkB")
    print(head(H_t2g))
    print(species)
    return(H_t2g)
  }else return(NULL)
}
gene_list_for_enrichment_genome <- function(H_t2g, Species=NULL){
  df <- list()
  set <- unique(H_t2g$gs_name)
  for(name in set){
    data <- dplyr::filter(H_t2g, gs_name == name)
    if(Species == "Drosophila melanogaster (dm6)"){
      column <- c("ENTREZID", "ENSEMBL")
      key <- "ENTREZID"
    gene_IDs<-AnnotationDbi::select(org(Species),keys = as.character(data$entrez_gene),
                                    keytype = key,
                                    columns = column)
    colnames(gene_IDs) <- c("entrez_gene","ENSEMBL")
    data <- gene_IDs %>% distinct(ENSEMBL, .keep_all = T)
    df[[name]] <- data$ENSEMBL
    }
    df[[name]] <- data$entrez_gene
  }
  return(df)
}
dorothea <- function(species, confidence = "recommend",type){
  if(species == "Mus musculus (mm10)" || species == "Mus musculus (mm39)"){
    net <- dorothea::dorothea_mm
  }else{
    net <- dorothea::dorothea_hs
  }
  if(confidence == "recommend"){
    net2 <- net %>% filter(confidence != "D") %>% filter(confidence != "E")
  }else net2 <- net
  if(type == "DoRothEA regulon (activator)") net2 <- net2%>% filter(mor == 1)
  if(type == "DoRothEA regulon (repressor)") net2 <- net2%>% filter(mor == -1)
  my.symbols <- gsub("\\..*","", net2$target)
  gene_IDs<-AnnotationDbi::select(org(species),keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL", "ENTREZID"))
  colnames(gene_IDs) <- c("target", "ENTREZID")
  gene_IDs <- gene_IDs %>% distinct(target, .keep_all = T)
  gene_IDs <- na.omit(gene_IDs)
  net2 <- merge(net2, gene_IDs, by="target")
  net3 <- data.frame(gs_name = net2$tf, entrez_gene = net2$ENTREZID, target = net2$target, confidence = net2$confidence)
  net3 <- dplyr::arrange(net3, gs_name)
  if(species != "Mus musculus (mm10)"  && species != "Mus musculus (mm39)" && species != "Homo sapiens (hg19)" && 
     species != "Homo sapiens (hg38)"){
    withProgress(message = paste0("Gene ID conversion from human to ", species, "for the regulon gene set. It takes a few minutes."),{
      genes <- net3$entrez_gene
      switch (species,
              "Rattus norvegicus (rn6)" = set <- "rnorvegicus_gene_ensembl",
              "Xenopus tropicalis" = set <- "xtropicalis_gene_ensembl",
              "Drosophila melanogaster (dm6)" = set <- "dmelanogaster_gene_ensembl",
              "Caenorhabditis elegans (ce11)" = set <- "celegans_gene_ensembl",
              "Bos taurus (bosTau8)" = set <- "btaurus_gene_ensembl",
              "Canis lupus familiaris" = set <- "clfamiliaris_gene_ensembl",
              "Danio rerio (danRer10)" = set <- "drerio_gene_ensembl",
              "Gallus gallus (galGal4)" = set <- "ggallus_gene_ensembl",
              "Macaca mulatta (rheMac8)" = set <- "mmulatta_gene_ensembl",
              "Pan troglodytes (panTro4)" = set <- "ptroglodytes_gene_ensembl",
              "Saccharomyces cerevisiae (sacCer3)" = set <-"scerevisiae_gene_ensembl",
              "Arabidopsis thaliana (tair10)" = set <- "athaliana_eg_gene")
      convert = useMart("ensembl", dataset = set, host="https://dec2021.archive.ensembl.org")
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
      genes2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
                      values = genes ,mart = human,
                      attributesL = c("entrezgene_id"),
                      martL = convert, uniqueRows=T)
      colnames(genes2) <- c("entrez_gene", "converted_entrez_gene")
      genes2 <- genes2 %>% distinct(converted_entrez_gene, .keep_all = T)
      merge <- merge(net3, genes2, by = "entrez_gene") 
      net3 <- data.frame(gs_name = merge$gs_name, entrez_gene = merge$converted_entrez_gene, confidence = merge$confidence)
      net3 <- dplyr::arrange(net3, gs_name)
    })
  }
  return(net3)
}
cnet_for_output <- function(data, plot_data, Gene_set, Species){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      if(Species != "not selected"){
        withProgress(message = "cnet plot",{
          p <- plot_data
          print(p)
          incProgress(1)
        })
      }else return(NULL)
    }
  }
}
range_changer <- function(data){
  chr<-gsub("\\:.+$", "", rownames(data))
  start <- gsub(".+\\:","",gsub("\\-.+$", "", rownames(data)) )
  end <- gsub(".+\\-","",rownames(data))
  data$chr <- chr
  data$start <- as.numeric(start)
  data$end <- as.numeric(end)
  return(data)
}
enrich_for_table <- function(data, H_t2g, Gene_set){
  if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
    return(NULL)
  }else{
    colnames(data)[1] <- "gs_name"
    H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
    data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
    if(Gene_set == "DoRothEA regulon (activator)" || Gene_set == "DoRothEA regulon (repressor)"){
      data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                          Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                          p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
    }else{
      if(Gene_set == "Custom gene set"){
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
      }else{
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, ID = data2$gs_id, Description = data2$gs_description,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
        
      }
      data3$Group <- gsub("\n"," ",data3$Group)
      return(data3) 
    }
  }
}


files2GRangelist <- function(files){
  Glist <- GRangesList()
  for(name in names(files)){
    extension <- gsub(".+\\.","",files$name)
    file_name <- gsub(".+\\/","",files$name)
    file_name<- gsub("\\..+$", "", file_name)
    print(file_name)
    print(files$name)
    input <- toGRanges(files$name, format=extension, header=FALSE) 
    Glist[[file_name]] <- input 
  }
  return(Glist)
}

integrate_ChIP_RNA <- function (result_geneRP, result_geneDiff, lfc_threshold = 1, 
                                padj_threshold = 0.05,name=NULL) {
  if ("GRanges" %in% class(result_geneRP)) {
    stop("sorry, please use the the simplify result or metadata(fullRP_hit)$peakRP_gene", 
         call. = FALSE)
  }
  merge_result <- dplyr::left_join(result_geneRP, result_geneDiff, 
                                   by = "gene_id")
  allGenes_N <- as.double(nrow(merge_result))
  if(length(name) <= 2 && length(merge_result$padj) != 0){
    up_name=name[2]
    down_name=name[1]
  merge_result <- merge_result %>% 
    dplyr::mutate(diff_rank = rank(padj, na.last = "keep"), 
                  diff_rank = dplyr::case_when(is.na(diff_rank) ~ allGenes_N, TRUE ~ diff_rank), 
                  rankProduct = RP_rank * diff_rank, rankOf_rankProduct = rank(rankProduct)) %>% 
    dplyr::arrange(rankOf_rankProduct) %>% 
    dplyr::mutate(gene_category = dplyr::case_when(log2FoldChange > lfc_threshold & padj < padj_threshold ~ up_name, 
                                                   log2FoldChange < -lfc_threshold & padj < padj_threshold ~ down_name, 
                                                   TRUE ~ "NS"), gene_category = factor(gene_category, levels = c(up_name,down_name, "NS")))
  upGenes_rank <- filter(merge_result, gene_category == up_name)$RP_rank
  downGenes_rank <- filter(merge_result, gene_category == down_name)$RP_rank
  staticGenes_rank <- filter(merge_result, gene_category == 
                               "NS")$RP_rank
  if (length(upGenes_rank) == 0 & length(downGenes_rank) == 
      0) {
    warning("no significant genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  else if (length(upGenes_rank) == 0) {
    warning("no significant up genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  else if (length(downGenes_rank) == 0) {
    warning("no significant down genes, just returing rank product result", 
            call. = FALSE)
    return(merge_result)
  }
  up_static_pvalue <- suppressWarnings(ks.test(upGenes_rank, 
                                               staticGenes_rank)$p.value)
  if(up_static_pvalue > 1) up_static_pvalue <- 1
  down_static_pvalue <- suppressWarnings(ks.test(downGenes_rank, 
                                                 staticGenes_rank)$p.value)
  if(down_static_pvalue > 1) down_static_pvalue <- 1
  ks_test <- paste0("\n Kolmogorov-Smirnov Tests ", "\n pvalue of ",up_name," vs NS: ", 
                    format(up_static_pvalue, digits = 3, scientific = TRUE), 
                    "\n pvalue of ", down_name," vs NS: ", format(down_static_pvalue, 
                                                                  digits = 3, scientific = TRUE))
  df <- data.frame(Group=c(paste0(up_name, " vs NS"), paste0(down_name, " vs NS")),
                   statistics = c(ks.test(upGenes_rank, staticGenes_rank)$statistic,ks.test(downGenes_rank, staticGenes_rank)$statistic),
                   pvalue=c(up_static_pvalue, down_static_pvalue))
  padj <- p.adjust(df$pvalue,method="BH")
  df$padj <- padj
  annotate_df <- data.frame(xpos = -Inf, ypos = Inf, annotateText = ks_test, 
                            hjustvar = 0, vjustvar = 1)
  p <- merge_result %>% ggplot2::ggplot(aes(x = RP_rank)) + 
    ggplot2::stat_ecdf(aes(color = gene_category), geom = "line") + 
    ggplot2::geom_text(data = annotate_df, aes(x = xpos, 
                                               y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText)) + 
    ggplot2::xlab("Regulatory potential rank\n(RP-high gene -> RP-low gene)") + ggplot2::ylab("Cumulative Probability")+
    ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::theme_bw(base_size = 15) + guides(color=guide_legend(title="Expression\nstatus"))+ 
    ggplot2::scale_color_manual(breaks = c(up_name,down_name,"NS"),values = c("#F8766D","#00BFC4","grey"))+theme(aspect.ratio = 1)
  }else{
    merge_result$Group[is.na(merge_result$Group)] <- "Other"
    NS_rank <- filter(merge_result, Group == "Other")$RP_rank
    Genes_rank <- list()
    static_pvalue <- list()
    static_D <- list()
    df <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
    for(name in unique(merge_result$Group)){
      if(name != "Other"){
      Genes_rank[name] <- list(dplyr::filter(merge_result, Group == name)$RP_rank)
      static_pvalue[name] <- list(ks.test(Genes_rank[[name]],NS_rank)$p.value)
      static_D[name] <- list(ks.test(Genes_rank[[name]],NS_rank)$statistic)
      df2 <- data.frame(Group = paste0(name," vs Other"), statistics = static_D[[name]],
                        pvalue = static_pvalue[[name]])
      df <- rbind(df, df2)
      }
    }
    padj <- p.adjust(df$pvalue,method="BH")
    df$padj <- padj
    annotate_df <- data.frame(xpos = -Inf, ypos = Inf, annotateText = NA, 
                              hjustvar = 0, vjustvar = 1)
    p <- merge_result %>% 
      ggplot2::ggplot(aes(x = RP_rank)) + 
      ggplot2::stat_ecdf(aes(color = Group), geom = "line") + 
      ggplot2::geom_text(data = annotate_df, aes(x = xpos, 
                                                 y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText)) + 
      ggplot2::xlab("Regulatory potential rank\n(RP-high gene -> RP-low gene)") + ggplot2::ylab("Cumulative Probability")+
      ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::theme_bw(base_size = 15) + guides(color=guide_legend(title="Expression\nstatus"))+ 
      theme(aspect.ratio = 1)
  }
  p$statistics <- df
  return(p)
}

enrich_viewer_forMulti2 <- function(data3, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set)){
    if(is.null(data3)){
      return(NULL)
    }else{
      withProgress(message = "enrichment analysis",{
        if(is.null(H_t2g)){
          df <- NULL
        }else{
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            sum <- length(data3$ENTREZID[data3$Group == name])
            em <- clusterProfiler::enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(clusterProfiler::setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                df <- rbind(df, cnet1)
              }
            }
          }
        }
        if(length(df$ID) !=0){
          df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
          return(df)
        }else return(NULL)
      })
    }
  } 
}

enrich_gene_list <- function(data, Gene_set, H_t2g, org){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      if(is.null(H_t2g)){
        df <- NULL
      }else{
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        df <- list()
        for (name in unique(data$Group)) {
          sum <- length(data$ENTREZID[data$Group == name])
          em <- clusterProfiler::enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
              df[[name]] <- cnet1
            }
          }
        }
      }
      return(df)
    }
  }
}

enrich_genelist <- function(data, enrich_gene_list, showCategory=5,type=NULL,group_order=NULL){
  if(is.null(data) || is.null(enrich_gene_list)){
    return(NULL)
  }else{
    df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
    cluster_list <- c()
    for (name in names(enrich_gene_list)) {
      sum <- length(data$ENTREZID[data$Group == name])
      if(!is.null(type)){
        if(type=="pairwise") sum <- length(dplyr::filter(data, group == name)$ENTREZID)
        }else sum <- length(data$ENTREZID[data$Group == name])
      em <- enrich_gene_list[[name]]
      if (length(as.data.frame(em)$ID) != 0) {
        if(length(colnames(as.data.frame(em))) == 9){
          cnet1 <- as.data.frame(em)
          cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
          cluster_list <- c(cluster_list, paste(name, "\n","(",sum, ")",sep = ""))
          if(!is.null(group_order)) group_order[which(group_order == name)] <- paste(name, "\n","(",sum, ")",sep = "")
          cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
          if (length(cnet1$pvalue) > showCategory){
            cnet1 <- cnet1[1:showCategory,]
          }
          df <- rbind(df, cnet1)
        }
      }
    }
    if(!is.null(group_order)) group_order <- group_order[group_order %in% cluster_list]
    if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
      p1 <- NULL
    } else{
      if(!is.null(type)){
        if(type == "withRNAseq"){
          df$Group <- gsub(":", ":\n", df$Group)
          if(!is.null(group_order)) group_order <- gsub(":", ":\n", group_order)
        }
      }
      if(!is.null(group_order)) {
        df$Group <- factor(df$Group, levels=group_order)
        df <- df %>% dplyr::arrange(Group) %>% dplyr::filter(!is.na(Group))
      }
      df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
      df <- dplyr::filter(df, !is.na(qvalue))
      df$Description <- gsub("_", " ", df$Description)
      df <- dplyr::mutate(df, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
      df$x <- gsub(":","", df$x)
      df <- dplyr::arrange(df, Group, x)
      idx <- order(df[["Group"]], df[["x"]], decreasing = FALSE)
      df$Description <- factor(df$Description,
                               levels=rev(unique(df$Description[idx])))
      df$log10padj <- -log10(df$qvalue)
      p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=log10padj,size=GeneRatio))+
                      geom_point() +
                      scale_color_continuous(low="blue", high="red",name="-log10\n(q value)") +
                      scale_size(range=c(1, 6),name="Gene\nratio")+ DOSE::theme_dose(font.size=12)+ylab(NULL)+xlab(NULL)+
                      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      p <- plot_grid(p1)
      return(p)
    }
  }
}

symbol2gene_id <- function(data,org){
  if(str_detect(rownames(data)[1], "FBgn")) {
    data$gene_id <- rownames(data)
  }else{
    my.symbols <- rownames(data)
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL","ENTREZID"))
  colnames(gene_IDs) <- c("SYMBOL","gene_id")
  gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
  gene_IDs <- data.frame(gene_id = gene_IDs$gene_id, row.names = gene_IDs$SYMBOL)
  data <- merge(gene_IDs,data,by=0)
  rownames(data)<-data$Row.names
  data <- data[,-1]
  }
  return(data)
}

data_trac <- function(y,gene_position,gen,txdb,org,filetype=NULL,bw_files,
                      bam_files,track_additional_files){
  chr <- as.character(gene_position$seqnames)
  print(chr)
  print(gen)
  grtrack <- GeneRegionTrack(txdb,
                             chromosome = chr, name = "UCSC known genes",geneSymbol = TRUE,
                             transcriptAnnotation = "symbol",genome = gen,
                             background.title = "grey",cex = 1.25)
  if(gen != "dm6"){
    ID <- "ENTREZID"
  }else ID <- "ENSEMBL"
  symbols <- unlist(mapIds(org, gene(grtrack), "SYMBOL", ID, multiVals = "first"))
  symbol(grtrack) <- symbols[gene(grtrack)]
  displayPars(grtrack) <- list(fontsize = 15)
  if(!is.null(filetype)){
  switch(filetype,
         "Row1" = bw_files <- bw_files,
         "Row2" = bw_files <- bam_files)
  }
  if(!is.null(track_additional_files)) bw_files <- c(bw_files, track_additional_files)
  df <- list()
  c <- 0
  unique <- unique(gsub("\\_.+$", "", names(bw_files)))
  num <- length(unique)
  col <- list()
  col_count = num + 1
  for(name in unique){
    col_count <- col_count - 1
    col[[name]] <- rainbow_hcl(num,c = 100)[col_count]
  }
  for(name in names(bw_files)){
    c <- c + 1
    if(!is.null(filetype)){
    if(filetype == "Row2"){ 
      bai <- file.exists(paste0(bw_files[[name]],".bai"))
      if(bai == FALSE) a <- indexBam(bw_files[[name]])
    }}
    name2 <- gsub("\\_.+$", "", name)
    df[[name]] <- DataTrack(range = bw_files[[name]], type = "l",genome = gen,
                            name = gsub("\\..+$", "", name), window = -1,
                            chromosome = chr, background.title = col[[name2]],cex = 0.8,
                            col.histogram = col[[name2]], cex.axis=0.8,cex.main=0.8, cex.title = 0.8,
                            fill.histogram = col[[name2]])
  }
  df[["grtrack"]] <- grtrack
  return(df)
}

RNAseqDEGimport <- function(tmp,exampleButton){
  withProgress(message = "Importing an RNA-seq DEG regult file, please wait",{
    if(is.null(tmp) && exampleButton > 0 )  tmp = "data/RNAseq.txt"
    if(is.null(tmp)) {
      return(NULL)
    }else{
      if(tools::file_ext(tmp) == "xlsx") {
        df2 <- readxl::read_xlsx(tmp) 
        df2 <- as.data.frame(df2)
        df <- try(data.frame(row.names = df2[,1]),silent = T)
        if(class(df) != "try-error") {
          if(dim(df2)[2] == 2){
            df <- data.frame(row.names = df2[,1],a = df2[,2])
            colnames(df)[1] <- colnames(df2)[2]
          }else{
            rownames(df2) <- df2[,1]
            df <- df2[,-1]
            colnames(df) <- gsub("-",".",colnames(df))
          }
        }
      }
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
      rownames(df) = gsub("\"", "", rownames(df))
      if(length(colnames(df)) != 0){
        if(str_detect(colnames(df)[1], "^X\\.")){
          colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
        }
      }
      return(df)
    }
  })
}
RNAseqDEG_ann <- function(RNAdata,Species,gene_type,input_type=NULL){
  if(input_type != "List") RNAdata$log2FoldChange <- -RNAdata$log2FoldChange
  if(str_detect(rownames(RNAdata)[1], "FBgn")){
    RNAdata$gene_id <- rownames(RNAdata)
    data <- RNAdata
  }else{
  if(gene_type != "SYMBOL"){
    my.symbols <- gsub("\\..*","", rownames(RNAdata))
    gene_IDs<-id_convert(my.symbols,Species,type="ENSEMBL")
    colnames(gene_IDs) <- c("EnsemblID","Symbol","gene_id")
    RNAdata$EnsemblID <- gsub("\\..*","", rownames(RNAdata))
    gene_IDs <- gene_IDs %>% distinct(EnsemblID, .keep_all = T)
    data <- merge(RNAdata, gene_IDs, by="EnsemblID")
  }else{
    my.symbols <- rownames(RNAdata)
    gene_IDs<-id_convert(my.symbols, Species,type="SYMBOL_single")
    colnames(gene_IDs) <- c("Symbol", "gene_id")
    gene_IDs <- gene_IDs %>% distinct(Symbol, .keep_all = T)
    RNAdata$Symbol <- rownames(RNAdata) 
    data <- merge(RNAdata, gene_IDs, by="Symbol")
  }
  }
  return(data)
}

GetGRanges <- function (LoadFile, simple = FALSE,sepr = "\t", simplify = FALSE) {
    if (sum(class(LoadFile) == "character")) {
      RangesTable <- read.delim(LoadFile, sep = sepr, header = FALSE, 
                                comment.char = "#")
      if(str_detect(RangesTable[1,2], "tart") == TRUE || is.na(RangesTable[1,2])) {
        RangesTable <- RangesTable[-1,]
      }else RangesTable <- RangesTable
    }
    Chromosomes <- as.vector(RangesTable[, 1])
    Start <- as.numeric(as.vector(RangesTable[, 2]))
    End <- as.numeric(as.vector(RangesTable[, 3]))
    RegionRanges <- GRanges(seqnames = Chromosomes, ranges = IRanges(start = Start, 
                                                                     end = End))
    if (simple == FALSE) {
      if (ncol(RangesTable) > 4) {
        ID <- as.vector(RangesTable[, 4])
        Score <- as.vector(RangesTable[, 5])
        if (ncol(RangesTable) > 6) {
          Strand <- rep("*", nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[, 
                                                       -c(1:6)])
          elementMetadata(RegionRanges) <- cbind(ID, 
                                                 Score, Strand, RemainderColumn)
        }
        else {
          elementMetadata(RegionRanges) <- cbind(ID, 
                                                 Score)
        }
      }
    }
  return(RegionRanges)
}

id_convert <- function(my.symbols,Species,type){
  if(Species != "Drosophila melanogaster (dm6)") {
    if(type == "ENTREZID"){
      column <- c("SYMBOL", "ENTREZID")
      key <- "ENTREZID"
    }
    if(type == "SYMBOL_single"){
      column <- "ENTREZID"
      key <- "SYMBOL"
    }
    if(type == "SYMBOL_double"){
      column <- c("SYMBOL", "ENTREZID")
      key <- "SYMBOL"
    }
    if(type == "ENSEMBL"){
      column <- c("ENSEMBL","SYMBOL","ENTREZID")
      key <- "ENSEMBL"
    }
    if(type == "ENSEMBL2ENTREZID"){
      column <- c("ENSEMBL","ENTREZID")
      key <- "ENSEMBL"
    }
  }else {
    if(type == "ENTREZID"){
    column <- c("SYMBOL", "ENSEMBL")
    key <- "ENSEMBL"
    }
    if(type == "SYMBOL_single"){
      column <- "ENSEMBL"
      key <- "SYMBOL"
    }
    if(type == "SYMBOL_double"){
      column <- c("SYMBOL", "ENSEMBL")
      key <- "ENSEMBL"
    }
    if(type == "ENSEMBL"){
      column <- c("ENSEMBL","SYMBOL")
      key <- "ENSEMBL"
    }
    if(type == "ENSEMBL2ENTREZID"){
      column <- c("ENSEMBL","ENTREZID")
      key <- "ENSEMBL"
    }
  }
  gene_IDs<-AnnotationDbi::select(org(Species),keys = my.symbols,
                                  keytype = key,
                                  columns = column)
  if(Species == "Drosophila melanogaster (dm6)" && type == "ENSEMBL") gene_IDs$ENTREZID <- gene_IDs$ENSEMBL
  return(gene_IDs)
}
ref_for_GREAT <- function(Species){
  switch (Species,
          "Mus musculus (mm10)" = source <- "TxDb.Mmusculus.UCSC.mm10.knownGene",
          "Homo sapiens (hg19)" = source <- "TxDb.Hsapiens.UCSC.hg19.knownGene",
          "Homo sapiens (hg38)" = source <- "TxDb.Hsapiens.UCSC.hg38.knownGene",
          "Mus musculus (mm39)" = source <- "TxDb.Mmusculus.UCSC.mm39.refGene",
          "Drosophila melanogaster (dm6)" = source <- "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
          "Rattus norvegicus (rn6)" = source <- "TxDb.Rnorvegicus.UCSC.rn6.refGene",
          "Caenorhabditis elegans (ce11)" = source <- "TxDb.Celegans.UCSC.ce11.refGene",
          "Bos taurus (bosTau8)" = source <- "TxDb.Btaurus.UCSC.bosTau8.refGene",
          "Canis lupus familiaris (canFam3)" = source <- "TxDb.Cfamiliaris.UCSC.canFam3.refGene",
          "Danio rerio (danRer10)" = source <- "TxDb.Drerio.UCSC.danRer10.refGene",
          "Gallus gallus (galGal4)" = source <- "TxDb.Ggallus.UCSC.galGal4.refGene",
          "Macaca mulatta (rheMac8)" = source <- "TxDb.Mmulatta.UCSC.rheMac8.refGene",
          "Pan troglodytes (panTro4)" = source <- "TxDb.Ptroglodytes.UCSC.panTro4.refGene",
          "Saccharomyces cerevisiae (sacCer3)" = source <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
          "Xenopus laevis (xenLae2)" = source <- "xenLae2",
          "Arabidopsis thaliana (tair10)" = source <- "TxDb.Athaliana.BioMart.plantsmart51")
  return(source)
}


peak_pattern_function <- function(grange, files,rg = NULL,additional=NULL,plot=TRUE){
  feature.recentered <- ChIPpeakAnno::reCenterPeaks(grange, width=4000)
  if(!is.null(additional)) files <- c(files, additional)
  cvglists <- sapply(files, import,which=feature.recentered,as="RleList")
  names(cvglists) <- gsub(".+\\/","",gsub("\\..+$", "", names(files)))
  feature.center <- ChIPpeakAnno::reCenterPeaks(grange, width=1)
  sig <- ChIPpeakAnno::featureAlignedSignal(cvglists, feature.center,
                              upstream=2000, downstream=2000)
  if(plot == TRUE){
  if(is.null(rg)){
  rg <- c()
  for(name in names(sig)){
    rg <- c(rg, mean(sig[[name]][,50]))
  }
  rg <- max(rg) + max(rg)*0.1
  }
  range <- c()
  for(n in length(names(files))) {range <- c(range, rg)}
  heat <- ChIPpeakAnno::featureAlignedHeatmap(sig, feature.center,
                                upstream=2000, downstream=2000,
                                upper.extreme=range,color = brewer.pal(n = 9, name = 'OrRd'))
  line <- ChIPpeakAnno::featureAlignedDistribution(sig, feature.center,
                                     upstream=2000, downstream=2000,
                                     type="l")
  df <- list()
  df[["heat"]] <- heat
  df[["line"]] <- line
  return(df)
  }else return(sig)
}
bigwig_breakline <- function(bigwig){
  names(bigwig) <- gsub("-", "- ", names(bigwig))
  for(i in 1:length(names(bigwig))){
    names(bigwig)[i] <- paste(strwrap(names(bigwig)[i], width = 15),collapse = "\n")
  }
  names(bigwig) <- gsub(" ", "", names(bigwig))
  return(bigwig)
}

batch_lineplot <- function(files2,files_bw){
  line_list <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
  for(k in 1:length(names(files_bw))){
    cvg_com <- data.frame(matrix(rep(NA, 100), nrow=1))[numeric(0), ]
    line_com <- data.frame(matrix(rep(NA, 100), nrow=1))[numeric(0), ]
    for(i in 1:length(names(files2))){
      cvg_list <- peak_pattern_function(grange=files2[[i]],files = files_bw[[k]],plot=F)
      line_com <- rbind(line_com,apply(cvg_list[[1]],2,mean))
    }
    line_com_t <- t(line_com)
    rownames(line_com_t) <- seq(from=-2000,to=2000,length=100)
    colnames(line_com_t) <- names(files2)
    y <- melt(line_com_t)
    colnames(y) <- c("distance","bed","density")
    y$bigwig <- names(files_bw)[[k]]
    line_list <- rbind(line_list, y)
  }
  plot_list <- list()
  g <- ggplot(line_list, aes(x = distance, y = density, color = bed))+
    facet_wrap(~ bigwig, scales="free_y")+ labs(color = "intersection")
  g <- g + scale_color_nejm()+ geom_line() + 
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA))
  g <- g +  theme(legend.position = "top",strip.text.x = element_text(size = 15),
                  title = element_text(size = 15),text = element_text(size = 12),
                  axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),legend.text = element_text(size=15),
                  legend.background = element_rect(fill=NA,color=NA))
  plot_list[["bed"]] <- g
  g <- ggplot(line_list, aes(x = distance, y = density, color = bigwig))+
    facet_wrap(~ bed, scales="free_y")
  g <- g + scale_color_nejm(name=)+ geom_line() +
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA))
  g <- g +  theme(legend.position = "top",strip.text.x = element_text(size = 15),
                  title = element_text(size = 15),text = element_text(size = 12),
                  axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),legend.text = element_text(size=15),
                  legend.background = element_rect(fill=NA,color=NA))
  plot_list[["bigwig"]] <- g
  return(plot_list)
}
pdf_h <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_hsize <- 24.5
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_hsize <- 22.25
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_hsize <- 19
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_hsize <- 17.75
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_hsize <- 15.5
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_hsize <- 13.5
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_hsize <- 11
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_hsize <- 9
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_hsize <- 9
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_hsize <- 8
  if (length(rowlist) == 4) pdf_hsize <- 8
  if (length(rowlist) == 3) pdf_hsize <- 4
  if (length(rowlist) == 2) pdf_hsize <- 4
  if (length(rowlist) == 1) pdf_hsize <- 4
  if (length(rowlist) > 200) pdf_hsize <- 30
  return(pdf_hsize)
}
pdf_w <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_wsize <- 24.5
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_wsize <- 22.25
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_wsize <- 19
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_wsize <- 17.75
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_wsize <- 15.5
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_wsize <- 15.5
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_wsize <- 14
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_wsize <- 14
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_wsize <- 10
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_wsize <- 12
  if (length(rowlist) == 4) pdf_wsize <- 10
  if (length(rowlist) == 3) pdf_wsize <- 12
  if (length(rowlist) == 2) pdf_wsize <- 8
  if (length(rowlist) == 1) pdf_wsize <- 4
  if (length(rowlist) > 200) pdf_wsize <- 30
  return(pdf_wsize)
}

batch_heatmap <- function(files2,files_bw,maxrange=NULL,type=NULL,
                          color=c("white", "red"),signal="signal",withRNAseq=FALSE,ylim=NULL){
  ht_list <- NULL
  perc <- 0
  if(withRNAseq == F){
    names(files2) <- gsub(".bw$", "",names(files2))
    names(files2) <- gsub(".bigwig$", "",names(files2))
    names(files2) <- gsub(".BigWig$", "",names(files2))
    names(files2) <- gsub("_", " ",names(files2))
    for(i in 1:length(names(files2))){
      names(files2)[i] <- paste(strwrap(names(files2)[i], width = 15),collapse = "\n")
    }
  }
  withProgress(message = "Preparing heatmap",{
  for(k in 1:length(names(files_bw))){
    perc<-perc+1
    num_list <- c()
    glist <- list()
    for(i in 1:length(names(files2))){
      feature.recentered <- ChIPpeakAnno::reCenterPeaks(files2[[i]], width=1)
      num <- dim(as.data.frame(feature.recentered))[1]
      num_list <- c(num_list, rep(names(files2)[[i]],num))
      glist[[i]] <- feature.recentered
    }
    target <- unlist(as(glist,"GRangesList"))
    mat1 = normalizeToMatrix(import(files_bw[[k]],as="GRanges"), target = target, value_column = "score",
                             extend = 2000, mean_mode = "absolute", w = 20, keep = c(0, 0.99))
    axis_name <- c("-2kb","0","2kb")
    if(!is.null(type)){
    if(type == "Promoter"){
      axis_name <- c("-2kb","TSS","2kb")
    }}
    name <- gsub("\\..+$", "", names(files_bw)[[k]])
    name <- gsub("_", " ", name)
    name <- gsub("-", " ", name)
    name <- paste(strwrap(name, width = 8),collapse = "\n")
    name <- gsub(" ", "\\_", name)
    print("heat")
    print(unique(num_list))
    if(is.null(ylim)) {
      ht <- EnrichedHeatmap(mat1, split = num_list, col = color,name=name,
                            axis_name = axis_name,pos_line = F,
                            column_title =name,use_raster=TRUE, row_title_rot = 0,
                            top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 1:8, lty = 1),
                                                                                        axis_param = list(side = "right", facing = "inside"))))
    }else{

      ht <- EnrichedHeatmap(mat1, split = num_list, col = circlize::colorRamp2(c(0, ylim), color),name=name,
                            axis_name = axis_name,pos_line = F,
                            column_title =name,use_raster=TRUE, row_title_rot = 0,
                            top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = 1:8, lty = 1),ylim = c(0,ylim),
                                                                                        axis_param = list(side = "right", facing = "inside")))) 
    }
    ht_list <- ht_list + ht
    list <- list()
    list[["heatmap"]] <- ht_list
    list[["mat"]] <- mat1
    incProgress(1/length(names(files_bw)), message = paste("Finish '", names(files_bw)[k], "', ", 
                                                            perc, "/", length(names(files_bw)),sep = ""))
  }
  })
  return(list)  
}
lgd <- function(files2,withRNAseq=FALSE){
  if(withRNAseq == F){
  leg_name <- gsub("\\..+$", "",names(files2))
  leg_name <- gsub("_", " ",leg_name)
  for(i in 1:length(leg_name)){
    leg_name[i] <- paste(strwrap(leg_name[i], width = 15),collapse = "\n")
  }
  }else leg_name <- names(files2)
  print(leg_name)
  leg_name <- sort(leg_name)
  lgd = Legend(at = leg_name, title = "Group", 
               type = "lines", legend_gp = gpar(col = 1:8))
  return(lgd)
}

gene_type <- function(my.symbols,org,Species){
  if(Species != "not selected"){
      ENSEMBL<-try(AnnotationDbi::select(org,keys = my.symbols,
                                         keytype = "ENSEMBL",
                                         columns = c("ENSEMBL", "ENTREZID")))
      SYMBOL <-try(AnnotationDbi::select(org,keys = my.symbols,
                                         keytype = "SYMBOL",
                                         columns = c("SYMBOL", "ENTREZID")))
      if(class(ENSEMBL) == "try-error" && class(SYMBOL) != "try-error") {type <- "SYMBOL"
      }else if(class(ENSEMBL) != "try-error" && class(SYMBOL) == "try-error") {type <- "ENSEMBL"
      }else if(class(ENSEMBL) == "try-error" && class(SYMBOL) == "try-error") {validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
      }else{
        if(dim(ENSEMBL)[1] > dim(SYMBOL)[1]) type <- "ENSEMBL" else type <- "SYMBOL"
      }
  }else type <- "not selected"
  return(type)
}
consensus_kmeans = function(mat, centers, km_repeats) {
  partition_list = lapply(seq_len(km_repeats), function(i) {
    as.cl_hard_partition(kmeans(mat, centers))
  })
  partition_list = cl_ensemble(list = partition_list)
  partition_consensus = cl_consensus(partition_list)
  as.vector(cl_class_ids(partition_consensus)) 
}

gene_type_for_integrated_heatmap <- function(files, Species, org) {
  gene_types <- c()
  name_list <- c()
  if(!is.null(files)){
    for(name in names(files)){
      type <- gene_type(my.symbols=rownames(files[[name]]),org=org,Species=Species)
      gene_types <- c(gene_types, type)
      name_list <- c(name_list,name) 
    }
    names(gene_types) <- name_list
    return(gene_types)
  }
}

rnaseqDEGs_for_integrated_heatmap <- function(files,Species,gene_type){
  if(!is.null(files)){
    df <- files_name2ENTREZID(files = files,Species=Species,gene_type=gene_type)
    if(length(names(df)) != 1){
      matrix_list <- list()
      for (name in names(df)) {
        matrix <- as.data.frame(df[name])
        if(str_detect(colnames(matrix)[1], "ENTREZID")) {
          rownames(matrix) <- matrix[,1]
          matrix <- matrix[,-1]
        }
        matrix[is.na(matrix)] <- 0
        matrix <- merge(matrix, matrix, by = 0)[,-2:-(1 + length(colnames(matrix)))]
        matrix_list[[name]] <- matrix
      }
      base <- matrix_list[[1]]
      int_matrix <- lapply(matrix_list[-1], function(i) base <<- merge(base, i, by = "Row.names"))
      rownames(base) <- base$Row.names
      colnames(base) <- gsub("\\.y$", "", colnames(base))
      rna <- as.data.frame(base[,-1])
      print(head(rna))
    }else{
      rna <- df[[names(df)]]
      if(str_detect(colnames(rna)[1], "ENTREZID")) {
        rownames(rna) <- rna$ENTREZID
        rna <- rna[,-1]
      }
      rna[is.na(rna)] <- 0
      rna <- as.data.frame(rna)
    }
    rna <- dplyr::select(rna, contains("log2FoldChange"))
    return(rna)
  }
}
rnaseqCounts_for_integrated_heatmap <- function(files,Species,gene_type, pre_zscoring =NULL){
  if(!is.null(files)){
    df <- files_name2ENTREZID(files = files,Species=Species,gene_type = gene_type)
    if(length(names(df)) != 1){
      if(is.null(pre_zscoring)) validate("")
      matrix_list <- list()
      matrix_z_list <- list()
      for (name in names(df)) {
        matrix <- as.data.frame(df[name])
        if(str_detect(colnames(matrix)[1], "ENTREZID")) {
          rownames(matrix) <- matrix[,1]
          matrix <- matrix[,-1]
        }
        if(pre_zscoring == "TRUE"){
        matrix_z <- genescale(matrix, axis = 1, method = "Z")
        matrix_z[is.na(matrix_z)] <- 0
        matrix_z <- merge(matrix, matrix_z, by = 0)[,-2:-(1 + length(colnames(matrix)))]
        matrix_z_list[[name]] <- matrix_z
        }else {
        matrix_2 <- matrix
        matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
        matrix_list[name] <- list(matrix_3)
        }
      }
      if(pre_zscoring == "TRUE"){
        base_z <- matrix_z_list[[1]]
        int_matrix <- lapply(matrix_z_list[-1], function(i) base_z <<- merge(base_z, i, by = "Row.names"))
      }else{
        base_z <- matrix_list[[1]]
        int_matrix <- lapply(matrix_list[-1], function(i) base_z <<- merge(base_z, i, by = "Row.names"))
      }
      rownames(base_z) <- base_z$Row.names
      colnames(base_z) <- gsub("\\.y$", "", colnames(base_z))
      rna <- base_z[,-1]
      if(pre_zscoring != "TRUE"){
        rna <- genefilter::genescale(rna, axis = 1, method = "Z")
        rna[is.na(rna)] <- 0
      }
    }else{
      rna <- df[[names(df)]]
      if(str_detect(colnames(rna)[1], "ENTREZID")) {
        rownames(rna) <- rna$ENTREZID
        rna <- rna[,-1]
      }
      print(head(rna))
      rna <- genescale(rna, axis = 1, method = "Z")
      rna[is.na(rna)] <- 0
    }
    rna <- as.data.frame(rna)
    return(rna)
  }
}

pre_integrated_additional <- function(integrated_bw){
if(!is.null(integrated_bw)){
  if(length(list.files("./Volume/")) > 0){
    files <- integrated_bw
    names(files) <- bigwig_name(integrated_bw)
  }else{
    files<-c()
    name<-c()
    for(nr in 1:length(integrated_bw[, 1])){
      file <- integrated_bw[[nr, 'datapath']]
      name <- c(name, bigwig_name(integrated_bw[nr,]$name))
      files <- c(files,file)
    }
    names(files)<-name
  }
  return(bigwig_breakline(files))
}
}


donut_replot <- function(dat){
  if(dim(dat)[1]>1){
  labs <- list(geneLevel=c(promoter="Promoter",
                           geneDownstream="Downstream",
                           geneBody="Gene body",
                           distalIntergenic="Distal Intergenic"),
               ExonIntron=c(exon="Exon",
                            intron="Intron",
                            intergenic="Intergenic"),
               Exons=c(utr5="5' UTR",
                       utr3="3' UTR",
                       CDS="CDS",
                       otherExon="Other exon"))
  labelCols = c(promoter="#D55E00",
                geneDownstream="#E69F00",
                geneBody="#51C6E6",
                distalIntergenic="#AAAAAA",
                exon="#009DDA",
                intron="#666666",
                intergenic="#DDDDDD",
                utr5="#0072B2",
                utr3="#56B4E9",
                CDS="#0033BF",
                otherExon="#009E73",
                undefined="#FFFFFF")
  l <- unlist(unname(labs))
  l1 <- paste0(l, " (",
               round(dat$percentage[match(names(l), dat$type)]*100,
                     digits = 1), "%)")
  names(l1) <- names(l)
  l1 <- c(l1, undefined="")
  p <- ggplot(dat, 
              aes_string(x="category", y="percentage", 
                         fill="type"))+
    geom_col() +
    coord_polar("y") + 
    geom_text(data = subset(dat, !duplicated(dat$category)),
              aes_string(x="category", label="category"),
              y=1) +
    theme_void() + 
    scale_fill_manual(values = labelCols, labels=l1, name=NULL,
                      guide = guide_legend(reverse=TRUE))
  return(p)
  }else return(NULL)
}
chisquare_for_annotation <- function(annotation){
  if(length(names(annotation$peaks))!=1){
  df <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
  for(name in names(annotation$peaks)){
    df2 <- data.frame(annotation$peaks[[name]])
    if(dim(df2)[1] >0){
    df2$Group <- name
    df <- rbind(df, df2)
    }
  }
  res <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
  for(level in c("geneLevel","ExonIntron","Exons")){
    table <- df %>% dplyr::group_by_(level) %>% 
      dplyr::count(Group) %>% tidyr::spread(key=level, value=n) %>%
      as.data.frame()
    rownames(table) <- table$Group
    table <- table[,-1]
    table[is.na(table)] <- 0
    if(level == "geneLevel") level_rename <-"Gene level" else 
      if(level == "ExonIntron") level_rename <- "Exon/Intron/Intergenic" else
        if(level == "Exons") level_rename <- "Exon level"
    combination <- combn(x=rownames(table),m=2)
    print(table)
    for(i in 1:dim(combination)[2]){
      table2 <- table[combination[,i],]
      chisq <- chisq.test(table2)
      res2 <- data.frame(group = paste0(paste(combination[,i],collapse = " vs ")," (",level_rename,")"), 
                         statistic=chisq$statistic,
                         method="Chi-square test",
                         pvalue=chisq$p.value,
                         row.names = NULL)
      res <- rbind(res, res2)
    } 
  }
  res$padj<- p.adjust(res$pvalue,method = "BH")
  return(res)
  }
}
annotation_barplot <- function(annotation){
  df <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
  for(name in names(annotation$peaks)){
    df2 <- data.frame(Group = name, peakN = log10(length(annotation$peaks[[name]])))
    df <- rbind(df, df2)
  }
  p <- ggbarplot(df,x="Group",y="peakN",fill = "Group") +
    coord_flip()+ theme_bw(base_size = 15) + xlab(NULL) + ylab("log10(Number of peaks)")
  return(p)
}

plot_annoDistance <- function (mmAnno, title=NULL){
  if(!is.null(mmAnno)){
  mmAnno$abs_dist <- abs(mmAnno$distanceToTSS) + 1
  summary_value <- summary(mmAnno$abs_dist) - 1
  df <- data.frame(row.names = names(summary_value),summary=round(as.numeric(summary_value),digits = 2))
  p1 <- ggplot2::ggplot(data.frame(mmAnno), aes(x = abs_dist)) + 
    ggplot2::geom_histogram(bins = 50) + ggplot2::scale_x_log10() +
    ggplot2::theme_bw() + ggplot2::xlab("abs(distanceToTSS) + 1")+ 
    ggplot2::ylab("Number of peaks") + ggtitle(title)
  p <- plot_grid(p1,tableGrob(unname(df),theme = ttheme_minimal()),nrow = 1,rel_widths = c(3,1))
  return(p)
  }
}

GREAT_dotplot <- function(data,data3,df,type,group_order=NULL){
  cluster_list <- c()
  for(name in names(df)){
    if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
      if(type !="pairwise") group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]])) else group1 <- as.data.frame(df[[name]])
      if(length(group1$id)==0){
        group1 <- NULL
      }else{
        if(type=="pairwise"){
          sum <- length(dplyr::filter(data3, group == name)$ENTREZID)
        }else if(type=="pair") {
          sum <- length(dplyr::filter(data3, group == name)$group)
        }else{
          sum <- length(as.data.frame(data3[[name]])$start)
        }
        group1$Group <- paste(name, "\n(",sum,")",sep = "")
        if(!is.null(group_order)) group_order[which(group_order == name)] <- paste(name, "\n","(",sum, ")",sep = "")
        cluster_list <- c(cluster_list, paste(name, "\n","(",sum, ")",sep = ""))
        if (length(group1$p_adjust_hyper) > 5){
          group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
          group1 <- group1[1:5,]
        }
      }
    }else group1 <- NULL
    data <- rbind(data, group1)
  }
  if(!is.null(group_order)) group_order <- group_order[group_order %in% cluster_list]
  colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                         "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
  if(length(data$Description) != 0){
    data$Group <- gsub("_", " ", data$Group)
    for(i in 1:length(data$Group)){
      data$Group[i] <- paste(strwrap(data$Group[i], width = 15),collapse = "\n")
    }
    data$Group <- gsub(" \\(", "\n\\(", data$Group)
    if(!is.null(group_order)) {
      group_order <- gsub("_", " ", group_order)
      for(i in 1:length(group_order)){
        group_order[i] <- paste(strwrap(group_order[i], width = 15),collapse = "\n")
      }
      group_order <- gsub(" \\(", "\n\\(", group_order)
      data$Group <- factor(data$Group, levels=group_order)
      data <- data %>% dplyr::arrange(Group) %>% dplyr::filter(!is.na(Group))
    }
    data$Description <- gsub("_", " ", data$Description)
    data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
    data$x <- gsub(":","", data$x)
    data <- dplyr::arrange(data, x)
    idx <- order(data[["x"]], decreasing = FALSE)
    data$Description <- factor(data$Description,
                               levels=rev(unique(data$Description[idx])))
    data$padj <- -log10(data$p_adjust_hyper+(1e-20))
    p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="padj",size="fold_enrichment_hyper"))+
                    geom_point() +
                    scale_color_continuous(low="blue", high="red",name="-log10\n(padj)") +
                    scale_size(range=c(1, 6),name = "Enrichment\nratio")+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                    scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
  }else p1 <- NULL
  return(p1)
}



#withRNA-seq-------
modified_calcRP_TFHit <- function (mmAnno, Txdb, decay_dist = 1000, report_fullInfo = FALSE, 
                                   verbose = TRUE,mode="fcw_combined") 
{
  peak_scaned_center <- GenomicRanges::resize(mmAnno, width = 1, 
                                              fix = "center")
  all_gene_location <- GenomicFeatures::genes(Txdb)
  gene_scaned <- all_gene_location[mmAnno$gene_id]
  gene_scaned_TSS <- GenomicRanges::resize(gene_scaned, width = 1, 
                                           fix = "start")
  mmAnno$centerToTSS <- GenomicRanges::distance(peak_scaned_center, 
                                                gene_scaned_TSS)
  scan_result <- as.data.frame(mmAnno)
  if (verbose) {
    message(">> calculating RP using centerToTSS and TF hit\t\t", 
            format(Sys.time(), "%Y-%m-%d %X"))
  }
  if(length(mmAnno$Log2FoldChange) != 0) {
    if(mode == "fcw_combined" || mode == "fcw_separate"){
      print("RP definition")
      mmAnno$RP <- 2^(-mmAnno$centerToTSS/decay_dist/2^abs(mmAnno$Log2FoldChange)) 
    }else mmAnno$RP <- 2^(-mmAnno$centerToTSS/decay_dist)
  }else mmAnno$RP <- 2^(-mmAnno$centerToTSS/decay_dist)
  if (verbose) {
    message(">> merging all info together\t\t", format(Sys.time(), 
                                                       "%Y-%m-%d %X"))
  }
  peakRP_gene <- mcols(mmAnno) %>% data.frame(stringsAsFactors = FALSE) %>% 
    dplyr::group_by(gene_id) %>% dplyr::summarise(withPeakN = dplyr::n(), 
                                                  sumRP = sum(RP)) %>% dplyr::arrange(-sumRP) %>% dplyr::mutate(RP_rank = rank(-sumRP))
  if (verbose) {
    message(">> done\t\t", format(Sys.time(), "%Y-%m-%d %X"))
  }
  if (report_fullInfo) {
    metadata(mmAnno)$peakRP_gene <- peakRP_gene
    return(mmAnno)
  }
  else {
    return(peakRP_gene)
  }
}
mmAnno <- function(peak,genomic_region=NULL,txdb,peak_distance,mode=NULL,group_name,distribution,DAR=NULL){
  if(!is.null(peak)){
    if(genomic_region == "Promoter"){ 
      peak <- peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      peak <- with(peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
      if(length(as.data.frame(peak)$seqnames) != 0){
        mcols(peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(peak))))
        mmAnno_up <- mm_nearestGene(peak,txdb)
        mmAnno_up$Group <- group_name
      }
    }else{
      range <- peak_distance * 1000
      if(mode == "Gene_scan"){
        if(is.null(peak_distance)) validate("")
        if(genomic_region != "Promoter"){ 
          withProgress(message = paste0("Preparing annotation of ", group_name," peaks"),{
            anno <- distribution
            incProgress(1)
          })
          promoter <- subset(anno, geneLevel %in% "promoter")
          other <- subset(anno, ! geneLevel %in% "promoter")
          if(length(as.data.frame(promoter)$seqnames) != 0){
            mcols(promoter) <- DataFrame(feature_id = paste0("peak_", seq_len(length(promoter))))
            promoter_scan <- mm_geneScan(promoter, txdb,upstream = 2000,downstream = 100)
            promoter_scan <- promoter_scan %>% plyranges::filter(distanceToTSS >= -2000 & distanceToTSS <= 100)
          }else promoter_scan <- NULL
          if(length(as.data.frame(other)$seqnames) != 0){
            mcols(other) <- DataFrame(feature_id = paste0("peak_", seq_len(length(other))))
            other_scan <- mm_geneScan(other, txdb,upstream = range,downstream = range)
            other_scan <- other_scan %>% plyranges::filter(abs(distanceToTSS) <= range)
          }else other_scan <- NULL
          mmAnno_up <- plyranges::bind_ranges(promoter_scan, other_scan)
        }
        mmAnno_up$Group <- group_name
      }else{
        if(length(as.data.frame(peak)$seqnames) != 0){
          mcols(peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(peak))))
          mmAnno_up <- mm_nearestGene(peak,txdb)
          mmAnno_up <- mmAnno_up %>% plyranges::filter(abs(distanceToTSS) <= range)
          mmAnno_up$Group <- group_name
        }
      }
    }
    if(!is.null(DAR)){
      if(genomic_region == "Promoter") {
        DAR <- data.frame(locus = DAR$locus,log2FoldChange = DAR[,2])
      }else DAR <- data.frame(locus = rownames(DAR),log2FoldChange = DAR[,2])
      mmano <- as.data.frame(mmAnno_up)
      mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
      print("mmano")
      print(head(mmano))
      print(head(DAR))
      merge <- merge(mmano,DAR, by = "locus")
      merge2<-with(merge,GRanges(seqnames,IRanges(start,end),
                                 feature_id = feature_id,
                                 gene_id = gene_id, distanceToTSS = distanceToTSS,
                                 Group = Group, Log2FoldChange = log2FoldChange))
      print("merge2")
      print(head(merge2))
      return(merge2)
    }else return(mmAnno_up)
  }else return(NULL)
}
DAR_withRNAseq <- function(DAR,genomic_region=NULL,Species){
  if(genomic_region == "Promoter"){ 
    if(str_detect(rownames(DAR)[1], "FBgn")){
      DAR <- DAR
    }else{
      my.symbols <- rownames(DAR)
      gene_IDs<-id_convert(my.symbols, Species,type="SYMBOL_single")
      colnames(gene_IDs) <- c("Symbol", "gene_id")
      gene_IDs <- gene_IDs %>% distinct(Symbol, .keep_all = T) %>% distinct(gene_id, .keep_all = T)
      gene_IDs <- data.frame(row.names = gene_IDs$Symbol, gene_id = gene_IDs$gene_id)
      data <- merge(DAR, gene_IDs, by=0)
      locus <- as.data.frame(promoter_region())
      locus$locus <- paste0(locus$seqnames,":",locus$start,"-",locus$end)
      DAR <- merge(data,locus,by="gene_id")
      DAR <- DAR[,-2]
    }
  }  
  return(DAR)
} 

RP_f <- function(mmAnno,txdb, mode="combined_fcw"){
  if(!is.null(mmAnno)) {
    result_geneRP_up <- modified_calcRP_TFHit(mmAnno = mmAnno,Txdb = txdb,mode = mode)
  }else result_geneRP_up <- NULL
  if(!is.null(result_geneRP_up)) {
    result_geneRP <-result_geneRP_up
    result_geneRP <- result_geneRP %>% dplyr::arrange(-sumRP)
    result_geneRP$RP_rank <- rownames(result_geneRP) %>% as.numeric()
    return(result_geneRP)
  }else return(NULL)
}
regulatory_potential_f <- function(species,data,result_geneRP,DEG_fc,DEG_fdr,name){
  if(species != "not selected"){
    merge_data <- integrate_ChIP_RNA(
      result_geneRP = result_geneRP,
      result_geneDiff = data,lfc_threshold = log(DEG_fc,2),padj_threshold = DEG_fdr,
      name=name
    )
    return(merge_data)
  }
}
RNAseq_box <- function(RNA,rna_name,select,data_type,RP,mmAnno,statistics){
  cond1 <- gsub("\\_.+$", "", rna_name[1])
  cond2 <- gsub("\\_.+$", "", rna_name[2])
  RNA <- dplyr::filter(RNA, !is.na(gene_id))
  df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
  genelist <- list()
  bed_list <- list()
  for(name in select){
    data <- merge(RNA,RP[[name]], by="gene_id",all=T)
    if(data_type != "List") data <- dplyr::filter(data, baseMean != 0)
    data$group <- "Others"
    data$group[data$sumRP > 2] <- "2 < RP"
    data$group[data$sumRP <= 2 & data$sumRP > 1] <- "1 < RP < 2"
    data$group[data$sumRP <= 1 & data$sumRP > 0.1] <- "0.1 < RP < 1"
    data$group[data$sumRP <= 0.1 & data$sumRP > 0.01] <- "0.01 < RP < 0.1"
    data$group[data$sumRP <= 0.01 & data$sumRP > 0] <- "0 < RP < 0.01"
    level <- NULL
    col <- NULL
    if(dim(dplyr::filter(data, group == "Others"))[1] != 0) {
      level <- c(level,"Others")
      col <- c(col,"gray")
    }
    if(dim(dplyr::filter(data, group == "0 < RP < 0.01"))[1] != 0) {
      level <- c(level,"0 < RP < 0.01")
      col <- c(col,"blue")
    }
    if(dim(dplyr::filter(data, group == "0.01 < RP < 0.1"))[1] != 0) {
      level <- c(level,"0.01 < RP < 0.1")
      col <- c(col,"#00BFC4")
    }
    if(dim(dplyr::filter(data, group == "0.1 < RP < 1"))[1] != 0) {
      level <- c(level,"0.1 < RP < 1")
      col <- c(col,"lightgreen")
    }
    if(dim(dplyr::filter(data, group == "1 < RP < 2"))[1] != 0) {
      level <- c(level,"1 < RP < 2")
      col <- c(col,"#F8766D")
    }
    if(dim(dplyr::filter(data, group == "2 < RP"))[1] != 0) {
      level <- c(level,"2 < RP")
      col <- c(col,"red")
    }
    data$group <- factor(data$group,levels=level,ordered=TRUE)
    data$intersection <- as.factor(name)
    df <- rbind(df,data)
    if(data_type != "List") {
      genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,RNAlog2FC = -1*data$log2FoldChange,
                                     sumRP = data$sumRP, gene_id = data$gene_id)
    }else{
      genelist[[name]] <- data.frame(Symbol = data$Symbol, group = data$group,
                                     sumRP = data$sumRP, gene_id = data$gene_id)
    }
    bed_list1 <- list()
    for(name2 in unique(data$group)){
      if(name2 != "Others"){
        table <- genelist[[name]] %>% dplyr::filter(group == name2)
        print(head(table))
        gene <- table$gene_id
        up_peak3 <- NULL
        if(length(gene) != 0){
          up_peak <- subset(mmAnno[[name]], gene_id %in% gene)
          up_peak2 <- as.data.frame(up_peak)
          up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
          up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
          up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
          mcols(up_peak3) <- DataFrame(Group = paste0(name,":",name2))
          print(up_peak3)
          bed_list1[[name2]] <- up_peak3
        }
      }
    }
    bed_list[[name]] <- bed_list1
  }
  data <- df 
  for(name in select){
    check <- data %>% dplyr::filter(intersection == name) %>% 
      dplyr::filter(group != "Others") %>% summarise(n())
    if(check <= 1) data <- data %>% dplyr::filter(intersection != name)
  }
  data$intersection <- gsub("-","-\n",data$intersection)
  collist <- unique(data$group)
  
  col <-c("gray","blue","#00BFC4","lightgreen","#F8766D","red")
  if(data_type != "List"){
    data$log10FoldChange <- log10(2^data$log2FoldChange)
    data$log10FoldChange <- as.numeric(data$log10FoldChange)
    if(dim(data)[1] == 0) validate("boxplot: There are few genes with |RP| > 1")
    if (length(collist) >= 3){
      stat.test <- data %>% dplyr::group_by(intersection) %>% 
        tukey_hsd(log10FoldChange ~ group)
      stat.test <- stat.test %>% add_significance("p.adj")
      stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
    }else{
      group1 <- dplyr::filter(data, group == collist[1])
      group2 <- dplyr::filter(data, group == collist[2])
      if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
        stat.test <- data %>% dplyr::group_by(intersection) %>% 
          t_test(log10FoldChange ~ group)
        stat.test <- stat.test %>% add_significance()
        stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      }else stat.test <- NULL
    }
    p <- try(ggpubr::ggboxplot(data, x = "group", y = "log10FoldChange",
                               fill = "group", scales = "free", 
                               xlab = FALSE, ylab = paste0("RNAseq log10(",cond2,"/",cond1,")"))+theme_bw(base_size = 15)+
               xlab(NULL)+scale_fill_manual(values = col) + scale_x_discrete(labels = label_wrap_gen(8)))
    if(statistics !="not selected") p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)         
    p <- facet(p, facet.by = "intersection",
               panel.labs.background = list(fill = "transparent", color = "transparent"),
               scales = "free", short.panel.labs = T, panel.labs.font = list(size=15))
    stat.test <- stat.test[,-2]
    stat.test <- stat.test[,1:9]
    colnames(stat.test)[1] <- "Peak"
    df <- list()
    df[["plot"]] <- p
    df[["statistical_test"]] <- stat.test
  }else {
    df <- list()
    table2 <- data %>% group_by(group, intersection) %>%
      summarise(count = n(), .groups = "drop")
    p <- ggplot(table2,aes(x = group,y= count, fill= group)) +
      geom_col(position=position_dodge2(preserve = "single"))+scale_fill_manual(values = col) +
      theme_bw(base_size = 15)+facet_wrap(~intersection,scales = "free") +
      xlab(NULL) + ylab("Gene number") + scale_x_discrete(labels = label_wrap_gen(8))
    df[["plot"]] <- p
  }
  df[["genelist"]] <- genelist
  df[["bedlist"]] <- bed_list
  return(df)}

peak_population <- function(select,RNA,mmAnno,rna_name,data_type,DEG_fdr,DEG_fc,RNAseq_mode){
  df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
  for(name in select){
    mmano <- as.data.frame(mmAnno[[name]])
    mmano$locus <- paste0(mmano$seqnames,":",mmano$start,"-",mmano$end)
    up_name <- paste0(rna_name[2]," gene only")
    down_name <- paste0(rna_name[1]," gene only")
    merge <- merge(mmano,RNA, by = "gene_id")
    if(data_type != "List"){
      merge <- merge  %>% dplyr::mutate(Up_gene = if_else(padj < DEG_fdr & log2FoldChange > log(DEG_fc,2), 1, 0),
                                        Down_gene = if_else(padj < DEG_fdr & log2FoldChange < -log(DEG_fc,2), 1, 0),
                                        NS_gene = if_else(padj >= DEG_fdr | baseMean == 0 | is.na(padj) |
                                                            (padj <= DEG_fdr & abs(log2FoldChange) <= log(DEG_fc,2)), 1, 0))
      merge[is.na(merge)] <- 0
      merge2 <- merge %>% group_by(locus,Group) %>% summarise(Total_associated_gene = sum(Up_gene)+sum(Down_gene)+sum(NS_gene),Group = Group,
                                                              Up_gene = sum(Up_gene), Down_gene = sum(Down_gene), NS_gene = sum(NS_gene))
      table <- merge2 %>% dplyr::mutate(type = if_else(Up_gene > 0 & Down_gene == 0 & NS_gene == 0, up_name,
                                                       if_else(Up_gene == 0 & Down_gene > 0 & NS_gene == 0, down_name,
                                                               if_else(Up_gene == 0 & Down_gene == 0 & NS_gene > 0, "unaffected genes", "Either combination\nof the above"))))
      table$type <- factor(table$type,levels = c(down_name,up_name,"unaffected genes","Either combination\nof the above"))
    }else{
      merge2 <- merge %>% group_by(locus, Group.y) %>% summarise(Count = n()) 
      merge2 <- spread(merge2,Group.y,Count)
      merge2[is.na(merge2)] <- 0
      namelist <- unique(merge$Group.y)
      merge2$Total_associated_gene <- rowSums(merge2[,-1],na.rm = TRUE)
      vec <- c()
      merge3 <- merge2[,-1]
      merge3 <- merge3[,-length(colnames(merge3))]
      for(i in 1:dim(merge3)[1]){
        type <- merge3[i,]
        type <- paste(names(type[which(type != 0)]),collapse = "&")
        vec <- c(vec,type) 
      }
      merge2$type <- vec
      table <- merge2
    }
    table$Group <- paste0(table$Group," peak")
    if(max(table$Total_associated_gene) < 5) table$Total_associated_gene <- as.character(table$Total_associated_gene) else table$Total_associated_gene <- as.numeric(table$Total_associated_gene)
    table$intersection <- name
    df <- rbind(df, table)
  }
  table <- df
  table$intersection <- gsub("-","\n",table$intersection)
  
  p2 <- ggplot(table,aes(x = Total_associated_gene, fill = type))
  if(RNAseq_mode=="Nearest") p2 <- p2 + geom_bar(position = "stack",width=0.3) else p2 <- p2 + geom_bar(position = "stack")
  p2 <- p2 + theme_bw(base_size = 15)+facet_wrap(~intersection,scales = "free") +
    ylab("Number of peaks") +
    xlab("Number of associated genes")+guides(fill=guide_legend(title="associated_gene_type"))
  if(RNAseq_mode != "Nearest") p2 <- p2 +  scale_x_continuous(breaks= scales::pretty_breaks())
  if(data_type != "List"){
    col <- c("#00BFC4","#F8766D","grey","black")
    p2 <- p2 + scale_fill_manual(values = col,name=c(down_name,up_name,"unaffected genes","Either combination\nof the above"))
  }
  return(p2)
}




##Docker only----------
library(marge)
options('homer_path' = "/usr/local/homer")
check_homer()

read_known_results<-function (path, homer_dir = TRUE) {
  if (homer_dir == TRUE) {
    path <- paste0(path, "/knownResults.txt")
  }
  if (!file.exists(path)) {
    warning(paste("File", path, "does not exist"))
    return(NULL)
  }
  col_spec <- readr::cols("c", "c", "d", "-", "d", "d", "c", 
                          "d", "c")
  raw <- readr::read_tsv(path, col_types = col_spec)
  colnames(raw) <- c("motif_name", "consensus", "log_p_value", 
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct")
  tmp <- raw %>% tidyr::separate_("motif_name", c("motif_name", 
                                                  "experiment", "database"), "/", extra = "drop")
  parsed <- .parse_homer_subfields(tmp) %>% dplyr::mutate_at(vars(contains("pct")), 
                                                             .parse_pcts)
  known <- .append_known_pwm(parsed)
  names(known$motif_pwm) <- known$motif_name
  return(known)
}
.parse_homer_subfields <- function(motif_tbl) {
  cond <- stringr::str_detect(motif_tbl$motif_name, "/") %>%
    sum(., na.rm = TRUE) > 0
  if (cond == TRUE) {
    motif_tbl <- motif_tbl %>%
      tidyr::separate_('motif_name',
                       c('motif_name', 'experiment', 'database'),
                       '/', extra = "drop", fill = "right")
  }
  
  ## Detect if parentheses are present in motif_name
  ## to break apart into motif_name vs. motif_family
  cond <- stringr::str_detect(motif_tbl$motif_name, '\\(') %>%
    sum(., na.rm = TRUE) > 0
  if (cond == TRUE) {
    motif_tbl <- motif_tbl %>%
      tidyr::separate_('motif_name',
                       c('motif_name', 'motif_family'),
                       '\\(', extra = "drop", fill = "right")
    motif_tbl$motif_family <- stringr::str_replace(motif_tbl$motif_family, '\\)', '')
  }
  
  ## Detect If parentheses are present in experiment
  ## to break apart into experiment vs. accession
  if ("experiment" %in% colnames(motif_tbl)) {
    cond <- stringr::str_detect(motif_tbl$experiment, '\\(') %>%
      sum(., na.rm = TRUE) > 0
    if (cond == TRUE) {
      motif_tbl <- motif_tbl %>%
        tidyr::separate_('experiment',
                         c('experiment', 'accession'),
                         '\\(', extra = "drop", fill = "right")
      motif_tbl$accession <- stringr::str_replace(motif_tbl$accession, '\\)', '')
    }
  }
  
  return(motif_tbl)
}
.parse_pcts <- function(x) {
  stringr::str_replace(x, "%", "") %>%
    as.numeric() * 0.01
}
.append_known_pwm <- function(known_results) {
  data(HOMER_motifs, envir = environment())
  hm <- HOMER_motifs %>%
    ##        dplyr::filter(rlang::UQ(rlang::sym('log_odds_detection')) > 0) %>%
    dplyr::select("motif_name", "motif_family", "experiment", "accession",
                  "motif_pwm", "log_odds_detection")
  dplyr::inner_join(known_results, hm, 
                    by = c("motif_name", "motif_family", "experiment", "accession"))
}   

findMotif <- function(df,anno_data = NULL,Species,type = "Genome-wide",section=NULL,venn=NULL,
                      motif,size,back="random",bw_count=NULL,other_data=NULL,motif_length){
  ref <- gsub(".+\\(","",gsub(")", "", Species))
  setwd(tempdir())
  switch(motif,
         "known motif" = time <- "10 ~ 20",
         "known and de novo motifs" = time <- "20 ~ 30")
  withProgress(message = paste0("Motif analysis takes about ",time," min per group"),{
    if(type == "Genome-wide" || type == "Other") {
      group_name <- names(df)
      group_file <- length(names(df))
    }else{
      group_name <- unique(df$group)
      group_file <- length(unique(df$group))
    }
    print("start")
    print(head(df))
    perc <- 0
    df2 <- list()
    path <- paste0(format(Sys.time(), "%Y%m%d_%H%M_homer_"),back,"_size-",size)
    dir.create(path = path)
    for(name in group_name){
      group_dir <- paste0(path, "/",name)
      group_dir <- gsub(" < ","....",group_dir)
      dir.create(path = group_dir)
      perc <- perc + 1
      if(!is.null(anno_data)){
        if(type == "Genome-wide"){
          data <- df[[name]]
          data2 <- anno_data %>% dplyr::filter(locus %in% rownames(data))
          y <- with(data2, GRanges(seqnames = seqnames, 
                                   ranges = IRanges(start,end)))
        }else{
          data <- dplyr::filter(df, group == name)
          y <- subset(anno_data, gene_id %in% data$ENTREZID)
          y <- as.data.frame(y)
          y <- with(y, GRanges(seqnames = seqnames, 
                               ranges = IRanges(start,end)))
        }
      }else y<- df[[name]]
      y <- as.data.frame(y)
      print(head(y))
      if(dim(y)[1] != 0){
        if(type=="Genome-wide"){
          if(back == "random"){
            bg <-'automatic'
          }else if(back == "all_peaks"){
            data <- anno_data %>% dplyr::filter(! locus %in% rownames(data))
            data2 <- range_changer(data)
            bg <- data.frame(seqnames = data2$chr,start=data2$start,end=data2$end)
            write.table(bg,paste0(tempdir(),"bed.bed"),row.names = F,col.names = F,quote = F, sep = "\t")
            bg <- paste0(tempdir(),"bed.bed")
          }
        }
        if(type== "Promoter") {
          bg <- subset(anno_data, ! gene_id %in% data$ENTREZID) 
          bg <- as.data.frame(bg)
          bg <- with(bg, GRanges(seqnames = seqnames, 
                                 ranges = IRanges(start,end)))
          bg <- as.data.frame(bg)
          write.table(bg,paste0(tempdir(),"bed.bed"),row.names = F,col.names = F,quote = F, sep = "\t")
          bg <- paste0(tempdir(),"bed.bed")
        }
        if(type=="Other"){
          if(back == "random"){
            bg <-'automatic'
          }else if(back == "custom"){
            bg <- as.data.frame(other_data)
            write.table(bg,paste0(tempdir(),"bed.bed"),row.names = F,col.names = F,quote = F, sep = "\t")
            bg <- paste0(tempdir(),"bed.bed")
          }else{
            if(!is.null(venn)) bg <- venn else bg <- df
            if(length(names(bg)) > 1){
              bg <- soGGi:::runConsensusRegions(GRangesList(bg), "none")
              write.table(bg,paste0(tempdir(),"bed.bed"),row.names = F,col.names = F,quote = F, sep = "\t")
              bg <- paste0(tempdir(),"bed.bed")
            }
          }
        }
        print("background sequence")
        print(head(bg))
        switch(motif,
               "known motif" = motif_type <- TRUE,
               "known and de novo motifs" = motif_type <- FALSE)
        print("tested sequence")
        print(head(y))
        find_motifs_genome(
          y,
          path = group_dir,
          genome = ref, 
          motif_length = motif_length,
          scan_size = size,
          optimize_count = 25,
          background = bg,
          local_background = FALSE,
          only_known = motif_type, only_denovo = FALSE,
          cores = 2, cache = 500,
          overwrite = TRUE, keep_minimal = FALSE
        )
        df2[[name]] <- group_dir
      }
      incProgress(1/group_file, message = paste("Finish motif analysis of Group '", name, "', ", perc, "/", group_file,sep = ""))
    }
    return(df2)
  })
}

known_motif <- function(df){
  df2 <- data.frame(matrix(rep(NA, 14), nrow=1))[numeric(0), ]
  colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct", "log_odds_detection","Group")
  for(name in names(df)){
    if(file.exists(paste0(df[[name]],"/knownResults.txt")) == TRUE){
      known <-  as.data.frame(read_known_results(path = df[[name]]))[,-13]
      known$Group <- name
      df2 <- rbind(df2,known)
    }
  }
  colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                     "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct", "log_odds_detection","Group")
  return(df2)
}

denovo_motif <- function(df){
  df2 <- data.frame(matrix(rep(NA, 18), nrow=1))[numeric(0), ]
  colnames(df2) <- c("consensus","motif_name", "log_odds_detection","motif_id",  "log_p_value_detection", 
                     "tgt_num", "tgt_pct", "bgd_num", "bgd_pct"," log_p_value","fdr","tgt_pos","tgt_std","bgd_pos","bgd_std","strand_bias","multiplicity", "Group")
  for(name in names(df)){
    if(file.exists(paste0(df[[name]],"/homerResults.html")) == TRUE){
      known <-  as.data.frame(read_denovo_results(path = df[[name]]))[,-5]
      known$Group <- name
      df2 <- rbind(df2,known)
    }
  }
  return(df2)
}


homer_Motifplot <- function(df, showCategory=5,section=NULL,group_order=NULL){
  df2 <- data.frame(matrix(rep(NA, 15), nrow=1))[numeric(0), ]
  cluster_list <- c()
  for(name in names(df)){
    print(file.exists(paste0(df[[name]],"/knownResults.txt")))
    
    if(file.exists(paste0(df[[name]],"/knownResults.txt")) == TRUE){
      res <- as.data.frame(as.data.frame(read_known_results(path = df[[name]])))
      if(length(res$motif_name) != 0){
        name <- gsub("\\.\\.\\.\\."," < ",name)
        res$Group <- name
        cluster_list <- c(cluster_list, name)
        if(!is.null(group_order)) group_order[which(group_order == name)] <- name
        if(length(rownames(res)) > showCategory){
          res <- res[1:showCategory,]
        }
        df2 <- rbind(df2, res)
      }
    }
  }
  if(!is.null(group_order)) group_order <- group_order[group_order %in% cluster_list]
  if(length(df2$motif_name) == 0){
    return(NULL)
  }else{
    colnames(df2) <- c("motif_name", "motif_family", "experiment", "accession", "database", "consensus", "p_value",
                       "fdr", "tgt_num", "tgt_pct", "bgd_num", "bgd_pct","motif_pwm","log_odds_detection","Group")
    df2$Group <- order_reform(group_order = df2$Group,section = section)
    if(!is.null(group_order)) {
      group_order <- order_reform(group_order = group_order,section = section)
      print(head(df2$Group))
      print(head(group_order))
      df2$Group <- factor(df2$Group, levels=group_order)
      df2 <- df2 %>% dplyr::arrange(Group) %>% dplyr::filter(!is.na(Group))
    }
    df2 <- dplyr::mutate(df2, x = paste0(Group, 1/-log10(eval(parse(text = "p_value")))))
    df2$x <- gsub(":","", df2$x)
    df2 <- dplyr::arrange(df2, Group, x)
    idx <- order(df2[["Group"]], df2[["x"]], decreasing = FALSE)
    df2$motif_name <- factor(df2$motif_name,
                             levels=rev(unique(df2$motif_name[idx])))
    df2$log10pval <- -log10(df2$p_value)
    d <- ggplot(df2, aes(x = Group,y= motif_name,color=log10pval,size=log_odds_detection))+
      geom_point() +
      scale_color_continuous(low="blue", high="red",name="-log10\n(p value)") +
      scale_size(range=c(1, 6),name="log\nodds\nratio")+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) +
      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top")+
      theme(plot.margin=margin(l=-0.75,unit="cm"))
    df2 <- df2 %>% distinct(motif_name, .keep_all = T)
    pfm <- df2$motif_pwm
    pfm1 <- list()
    for(name in names(pfm)){
      pfm1[[name]] <- t(pfm[[name]])
    }
    Seqlogo <- try(as.grob(ggseqlogo(pfm1,ncol = 1)+ theme(axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           title = element_blank(),
                                                           text = element_blank()
                                                            ))) 
    p <- try(plot_grid(plot_grid(NULL, Seqlogo, ncol = 1, rel_heights = c(0.05:10)),as.grob(d)))
    if(length(p) == 1) if(p == "try-error") p <- d
    return(p)
  }
}

order_reform <- function(group_order=NULL,section =NULL){
  if(!is.null(section)){
    if(section == "withRNAseq"){
      if(!is.null(group_order)) {
        group_order <- gsub("--", ":\n", group_order)
        group_order <- gsub("_", " ", group_order)
        for(i in 1:length(group_order)){
          group_order[i] <- paste(strwrap(group_order[i], width = 8),collapse = "\n")
        }
      }
    }else if(section == "venn"){
      if(!is.null(group_order)){
        group_order <- gsub("- ", "- ", group_order)
        for(i in 1:length(group_order)){
          group_order[i] <- paste(strwrap(group_order[i], width = 8),collapse = "\n")
        }
      }
    }else{
      if(!is.null(group_order)){
        group_order <- gsub("_", " ", group_order)
        for(i in 1:length(group_order)){
          group_order[i] <- paste(strwrap(group_order[i], width = 8),collapse = "\n")
        }
      } 
    }
    if(!is.null(group_order)) group_order <- gsub(" \\(", "\n\\(", group_order)
  }
    return(group_order)
}

robust.system <- function (cmd) {
  stderrFile = tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
  stdoutFile = tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
  
  retval = list()
  retval$exitStatus = system(paste0(cmd, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
  retval$stdout = readLines(stdoutFile)
  retval$stderr = readLines(stderrFile)
  
  unlink(c(stdoutFile, stderrFile))
  return(retval)
}

homer_Motifplot_denovo <- function(df, name,showCategory=5,section=NULL,group_order=NULL){
  df2 <- data.frame(matrix(rep(NA, 15), nrow=1))[numeric(0), ]
  print(file.exists(paste0(df[[name]],"/homerResults.html")))
  
  if(file.exists(paste0(df[[name]],"/homerResults.html")) == TRUE){
    res <- read_denovo_html_results(path = df[[name]])
    motif_list <- read_denovo_results(path = df[[name]]) %>% 
      as.data.frame()
    motif_list <- motif_list[1:showCategory,]
    motif_list_for_pvalue <- data.frame(matrix(rep(NA, 15), nrow=1))[numeric(0), ]
    match_list <- data.frame(matrix(rep(NA, 15), nrow=1))[numeric(0), ]
    order <- c()
    for(i in 1:showCategory){
      ii <- as.character(i)
      order <- c(order, ii)
      motif_list_for_pvalue <- rbind(motif_list_for_pvalue,res[ii][[ii]]$Motif_information)
      match <- res[ii][[ii]]$Matches_to_known_motifs[1:5,]
      match$Group <- ii
      match_list <- rbind(match_list, match)
    }
  }
  if(length(motif_list$motif_name) == 0){
    return(NULL)
  }else{
    motif_list$motif_name <- factor(motif_list$motif_name,
                                    levels=rev(motif_list$motif_name))
    motif_list$log_p_value <- as.numeric(gsub("1e-","",motif_list_for_pvalue$p_value))
    motif_list$Group <- "1"
    d <- ggplot(motif_list, aes(x=Group,y= motif_name,color=log_p_value,size=log_odds_detection))+
      geom_point() +
      scale_color_continuous(low="blue", high="red",name="-log10\n(p value)") +
      scale_size(range=c(1, 12),name="log\nodds\nratio")+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) +
      theme(plot.margin=margin(l=-0.75,unit="cm"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            text = ggplot2::element_text(size = 15))
    pfm <- motif_list$motif_pwm
    pfm1 <- list()
    for(name in names(pfm)){
      pfm1[[name]] <- t(pfm[[name]])
    }
    Seqlogo <- try(as.grob(ggseqlogo(pfm1,ncol = 1)+ theme(axis.text.x = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           axis.title.y = element_blank(),
                                                           title = element_blank(),
                                                           text = element_blank()
    ))) 
    p <- try(plot_grid(Seqlogo,plot_grid(NULL, as.grob(d),nrow = 1, rel_widths = c(0.1:10000))))
    if(length(p) == 1) if(p == "try-error") p <- d
    
    
    bar_list <- list()
    for(ii in unique(match_list$Group)){
      match_list2 <- match_list %>% dplyr::filter(Group==ii)
      match_list2$xlab <- paste0(match_list2$rank,". ",match_list2$motif_name)
      match_list2$xlab <- factor(match_list2$xlab,levels = match_list2$xlab)
      match_list2$score <- as.numeric(match_list2$score)
      bar_list[[ii]] <- ggbarplot(match_list2,x="xlab",y="score",
                                  fill = "motif_name") + 
        scale_x_discrete(limit=rev(match_list2$xlab))+
        coord_flip()+ xlab(NULL) + ylab(NULL) + ylim(c(0,1))+
        theme(legend.position="none")
    }
    if(length(unique(match_list$Group))==10) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],bar_list[["6"]],bar_list[["7"]],bar_list[["8"]],bar_list[["9"]],bar_list[["10"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==9) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],bar_list[["6"]],bar_list[["7"]],bar_list[["8"]],bar_list[["9"]],bar_list[["10"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==8) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],bar_list[["6"]],bar_list[["7"]],bar_list[["8"]],bar_list[["9"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==7) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],bar_list[["6"]],bar_list[["7"]],bar_list[["8"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==6) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],bar_list[["6"]],bar_list[["7"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==5) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],bar_list[["5"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==4) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],bar_list[["4"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==3) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],bar_list[["3"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==2) bar <- try(plot_grid(bar_list[["1"]],bar_list[["2"]],ncol = 1,align = "v"))
    if(length(unique(match_list$Group))==1) bar <- try(plot_grid(bar_list[["1"]],ncol = 1,align = "v"))
    c <- plot_grid(bar,p,nrow = 1)
    return(c)
  }
}

