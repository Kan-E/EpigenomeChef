library(rtracklayer) 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(shiny)
library(DT)
library(gdata)
library(rstatix)
library(multcomp)
library(tidyverse)
library(tools)
library(ggpubr)
library(venn)
library(ggrepel)
library(ggdendro)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(ggnewscale)
library(IHW)
library(qvalue)
library(DEGreport)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(genefilter)
library(ComplexHeatmap)
library(shinyBS, verbose = FALSE)
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library(BiocManager)
library(clusterProfiler.dplyr)
library(dorothea)
library(umap)
library(biomaRt)
library(monaLisa)
library(GenomicRanges)
library(BiocParallel)
library(SummarizedExperiment)
library(JASPAR2020)
library(soGGi) ##devtools::install_github("ColeWunderlich/soGGi")
library(ChIPQC)
library(ChIPseeker)
library(ChIPpeakAnno)
library(rGREAT)
library(FindIT2)
options(repos = BiocManager::repositories())
options(rsconnect.max.bundle.size=3145728000)
species_list <- c("not selected", "Mus musculus (mm10)","Homo sapiens (hg19)","Homo sapiens (hg38)")
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target")
gene_set_list_genome <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "Transcription factor targets", "miRNA target")
org <- function(Species){
  if(Species != "not selected"){
    switch (Species,
            "Homo sapiens (hg38)" = org <- org.Hs.eg.db,
            "Homo sapiens (hg19)" = org <- org.Hs.eg.db,
            "Mus musculus (mm10)" = org <- org.Mm.eg.db
            )
    return(org)
  }
}

txdb_function <- function(Species){
  switch (Species,
          "Mus musculus (mm10)" = txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene,
          "Homo sapiens (hg19)" = txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene,
          "Homo sapiens (hg38)" = txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene)
  return(txdb)
}

promoter <- function(txdb, upstream, downstream,input_type = "Promoter",files = NULL){
  if(input_type == "Promoter"){
  return(promoters(genes(txdb),upstream = upstream, downstream = downstream))
  }else{
    consensusToCount <- soGGi:::runConsensusRegions(GRangesList(files), "none")
    occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
      rowSums
    consensusToCount <- consensusToCount[occurrences >= 2, ]
    return(consensusToCount)
  }
  }

Bigwig2count <- function(bw, promoter, Species, input_type = "Promoter"){

  switch (Species,
          "Mus musculus (mm10)" = org <- org.Mm.eg.db,
          "Homo sapiens (hg19)" = org <- org.Hs.eg.db,
          "Homo sapiens (hg38)" = org <- org.Hs.eg.db)
bed1<-promoter
write.table(bed1,file = "bed.bed",sep = "\t")
bed1 <- read.table("bed.bed",header = T)
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
  coverage <- import(bw[[i]], as = 'RleList')
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
  rownames(counts) <- bed1$gene_id
my.symbols <- rownames(counts)
gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                keytype = "ENTREZID",
                                columns = c("ENTREZID","SYMBOL"))
colnames(gene_IDs) <- c("Row.names","SYMBOL")
gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
gene_IDs <- data.frame(SYMBOL = gene_IDs$SYMBOL, row.names = gene_IDs$Row.names)
data2 <- merge(gene_IDs,counts, by=0)
rownames(data2) <- data2$SYMBOL
counts <- data2[,-1:-2]
}else{
  a <- as.data.frame(bed)
  Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
  rownames(counts) <- Row.name
}
return(counts)

}


PCAplot <- function(data, plot){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
  }
  pca <- prcomp(data, scale. = T)
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$rotation <- as.data.frame(pca$rotation)
  if(plot==TRUE){
  g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                y=pca$rotation[,2],
                                col=label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab(lab_x) + ylab(lab_y) + geom_text_repel()  +
    theme(legend.position="none", aspect.ratio=1)
  rho <- cor(data,method="spearman")
  d <- dist(1-rho)
  mds <- as.data.frame(cmdscale(d))
  label<-colnames(data)
  label<-gsub("\\_.+$", "", label)
  g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                        col = label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab("dim 1") + ylab("dim 2") +
    geom_text_repel() + theme(legend.position="none", aspect.ratio=1)
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

GOIheatmap <- function(data.z, show_row_names = TRUE){
  ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                clustering_method_columns = 'ward.D2',
                show_row_names = show_row_names, show_row_dend = F,column_names_side = "top",
                row_names_gp = gpar(fontface = "italic"))
  return(ht)
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
  if(Species != "not selected" || is.null(Gene_set) || is.null(org)){
    switch (Species, 
            "Mus musculus (mm10)" = species <- "Mus musculus",
            "Homo sapiens (hg19)" = species <- "Homo sapiens",
            "Homo sapiens (hg38)" = species <- "Homo sapiens")
    if(Gene_set == "MSigDB Hallmark"){
      H_t2g <- msigdbr(species = species, category = "H") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
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
    return(H_t2g)
  }else return(NULL)
}
gene_list_for_enrichment_genome <- function(H_t2g){
  df <- list()
  set <- unique(H_t2g$gs_name)
  for(name in set){
    data <- dplyr::filter(H_t2g, gs_name == name)
    df[[name]] <- data$entrez_gene
  }
  return(df)
}
dorothea <- function(species, confidence = "recommend",type){
  if(species == "Mus musculus (mm10)"){
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
      return(data3) 
    }
  }
}

pwms <- function(Species){
  library(TFBSTools)
  if(Species == "Mus musculus (mm10)"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome = BSgenome.Mmusculus.UCSC.mm10
    tax <- 10090
  }
  if(Species == "Homo sapiens (hg19)"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    tax <- 9606
  }
  if(Species == "Homo sapiens (hg38)"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome = BSgenome.Hsapiens.UCSC.hg38
    tax <- 9606
  }
  pwms_res <- getMatrixSet(JASPAR2020,
                       opts = list(matrixtype = "PWM",
                                   tax_group = "vertebrates",
                                   species = tax
                       ))
  return(pwms_res)
}

MotifAnalysis <- function(df, anno_data =NULL,Species,pwms, type ="Genome-wide",consensus = NULL){
  withProgress(message = "Motif analysis takes about 2 min per group",{
    library(TFBSTools)
    if(Species == "Mus musculus (mm10)"){
      library(BSgenome.Mmusculus.UCSC.mm10)
      genome = BSgenome.Mmusculus.UCSC.mm10
      tax <- 10090
    }
    if(Species == "Homo sapiens (hg19)"){
      library(BSgenome.Hsapiens.UCSC.hg19)
      genome = BSgenome.Hsapiens.UCSC.hg19
      tax <- 9606
    }
    if(Species == "Homo sapiens (hg38)"){
      library(BSgenome.Hsapiens.UCSC.hg38)
      genome = BSgenome.Hsapiens.UCSC.hg38
      tax <- 9606
    }
    
    perc <- 0
    df2 <- list()
    if(type == "Genome-wide") {
      group_name <- names(df)
      group_file <- length(names(df))
    }else{
    group_name <- unique(df$group)
    group_file <- length(unique(df$group))
    }
    for(name in group_name){
      perc <- perc + 1
      if(!is.null(anno_data)){
      if(type == "Genome-wide"){
        data <- df[[name]]
      data2 <- anno_data %>% dplyr::filter(locus %in% rownames(data))
      y <- with(data2, GRanges(seqnames = seqnames, 
                           ranges = IRanges(start,end)))
      }else{
        data <- dplyr::filter(df, group == name)
        print(data)
        print(anno_data)
        y <- subset(anno_data, gene_id %in% data$ENTREZID)
      }
      }else y<- df[[name]]
      y <- reCenterPeaks(y, width=mean(as.data.frame(y)$width))
      seq <- getSeq(genome, y)
      print(seq)
      print(pwms)
      if(type == "Promoter") {
        genome.region = subset(anno_data, ! gene_id %in% data$ENTREZID)
      }else{
        if(!is.null(consensus)){
          consensus_nega <- setdiff(consensus,y)
          genome.region = reCenterPeaks(consensus_nega,width=mean(as.data.frame(y)$width))
          print(genome.region)
        }else genome.region = NULL
      }
      se <- calcBinnedMotifEnrR(seqs = seq,
                                pwmL = pwms,
                                background = "genome",
                                genome = genome,
                                genome.oversample = 2,
                                genome.regions = genome.region,
                                BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                                verbose = TRUE)
      res <- data.frame(motif.id = elementMetadata(se)$motif.id, motif.name = elementMetadata(se)$motif.name,
                        motif.percentGC = elementMetadata(se)$motif.percentGC,
                        negLog10P = assay(se,"negLog10P"),negLog10Padj = assay(se,"negLog10Padj"), 
                        log2enr = assay(se,"log2enr"),pearsonResid = assay(se,"pearsonResid"),
                        expForegroundWgtWithHits = assay(se,"expForegroundWgtWithHits"),
                        sumForegroundWgtWithHits = assay(se,"sumForegroundWgtWithHits"),
                        sumBackgroundWgtWithHits = assay(se,"sumBackgroundWgtWithHits"),
                        Group = name)
      df2[[name]] <- res
      incProgress(1/group_file, message = paste("Finish motif analysis of Group '", name, "', ", perc, "/", group_file,sep = ""))
    }
    return(df2)
  })
}

MotifRegion <- function(df, anno_data =NULL, target_motif, Species, type="Genome-wide"){
  if(Species == "Mus musculus (mm10)"){
    genome = BSgenome.Mmusculus.UCSC.mm10
  }
  if(Species == "Homo sapiens (hg19)"){
    genome = BSgenome.Hsapiens.UCSC.hg19
  }
  if(Species == "Homo sapiens (hg38)"){
    genome = BSgenome.Hsapiens.UCSC.hg38
  }
  name <- gsub("\\\n.+$", "", target_motif$Group)
  if(!is.null(anno_data)){
  if(type == "Genome-wide"){
    data <- df[[name]]
  data2 <- anno_data %>% dplyr::filter(locus %in% rownames(data))
  y <- with(data2, GRanges(seqnames = seqnames, 
                           ranges = IRanges(start,end),
                           locus = locus))
  }else{
    data <- dplyr::filter(df, group == name)
    y <- subset(anno_data, gene_id %in% data$ENTREZID)
  }
  }else y <- df[[name]]
  if(length(rownames(as.data.frame(y))) == 0) stop("Incorrect species")
  seq <- getSeq(genome, y)
  if(!is.null(anno_data)){
  if(type == "Genome-wide") names(seq)<- y$locus
  }else{
    anno <- as.data.frame(y)
    Row.name <- paste0(anno$seqnames,":",anno$start,"-",anno$end)
    rownames(seq) <- Row.name
    rownames(anno) <- Row.name
  }
  print(seq)
  pfm <- getMatrixByID(JASPAR2020,target_motif$motif.id)
  pwm <- toPWM(pfm)
  res <- findMotifHits(query = pwm,
                       subject = seq,
                       min.score = 6.0,
                       method = "matchPWM",
                       BPPARAM = BiocParallel::SerialParam()) %>% as.data.frame()
  print(res)
  if(!is.null(anno_data)){
  if(type == "Genome-wide"){
  anno <- data.frame(seqnames = data2$locus, NearestGene = data2$NearestGene)
  }else anno <- data.frame(Gene = data$Row.names, seqnames = data$ENTREZID)
  res2<-merge(anno,res,by="seqnames")
  if(type == "Genome-wide"){
  colnames(res2)[1] <-"locus"
  start <- data2$start + res2$start
  end <- data2$start + res2$end
  res2$start<- start
  res2$end <- end
  }
  }else res2 <- res
  return(res2)
}

Motifplot <- function(df2, padj = 0.05){
  df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
  for(name in names(df2)){
    res <- df2[[name]]
    res <- dplyr::filter(res, X1.1 > -log10(padj))
    res <- res %>% dplyr::arrange(-X1.1)
    if(length(rownames(res)) > 5){
      res <- res[1:5,]
    }
    df <- rbind(df, res)
  }
  colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                    "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                    "Group")
  if(length(df$motif.id) == 0){
    return(NULL)
  }else{
    df$padj <- 10^(-df$negLog10Padj)
    df <- dplyr::mutate(df, x = paste0(Group, 1/-log10(eval(parse(text = "padj")))))
    df$x <- gsub(":","", df$x)
    df <- dplyr::arrange(df, x)
    idx <- order(df[["x"]], decreasing = FALSE)
    df$motif.name <- factor(df$motif.name,
                            levels=rev(unique(df$motif.name[idx])))
    d <- ggplot(df, aes(x = Group,y= motif.name,color=padj,size=log2enr))+
      geom_point() +
      scale_color_continuous(low="red", high="blue",
                             guide=guide_colorbar(reverse=TRUE)) +
      scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) +
      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top")+
      theme(plot.margin=margin(l=-0.75,unit="cm"))
    
    df <- df %>% distinct(motif.id, .keep_all = T)
    width.seqlogo = 2
    highlight <- NULL
    clres <- FALSE
    optsL <- list(ID = df$motif.id)
    pfm1 <- TFBSTools::getMatrixSet(JASPAR2020, opts = optsL)
    maxwidth <- max(vapply(TFBSTools::Matrix(pfm1), ncol, 0L))
    grobL <- lapply(pfm1, seqLogoGrob, xmax = maxwidth, xjust = "center")
    hmSeqlogo <- HeatmapAnnotation(logo = annoSeqlogo(grobL = grobL, 
                                                      which = "row", space = unit(1, "mm"),
                                                      width = unit(width.seqlogo, "inch")), 
                                   show_legend = FALSE, show_annotation_name = FALSE, 
                                   which = "row")
    tmp <- matrix(rep(NA, length(df$motif.id)),ncol = 1, 
                  dimnames = list(df$motif.name, NULL))
    
    hmMotifs <- Heatmap(matrix = tmp, name = "names", width = unit(0, "inch"), 
                        na_col = NA, col = c(`TRUE` = "green3",`FALSE` = "white"), 
                        cluster_rows = clres, show_row_dend = show_dendrogram, 
                        cluster_columns = FALSE, show_row_names = TRUE, row_names_side = "left", 
                        show_column_names = FALSE, show_heatmap_legend = FALSE,
                        left_annotation = hmSeqlogo)
    h <- grid.grabExpr(print(hmMotifs),wrap.grobs=TRUE)
    p <- plot_grid(plot_grid(NULL, h, ncol = 1, rel_heights = c(0.05:10)),as.grob(d))
    
    return(p)
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
                                padj_threshold = 0.05) {
  if ("GRanges" %in% class(result_geneRP)) {
    stop("sorry, please use the the simplify result or metadata(fullRP_hit)$peakRP_gene", 
         call. = FALSE)
  }
  merge_result <- dplyr::left_join(result_geneRP, result_geneDiff, 
                                   by = "gene_id")
  allGenes_N <- as.double(nrow(merge_result))
  merge_result <- merge_result %>% dplyr::mutate(diff_rank = rank(padj, 
                                                                  na.last = "keep"), diff_rank = dplyr::case_when(is.na(diff_rank) ~ 
                                                                                                                    allGenes_N, TRUE ~ diff_rank), rankProduct = RP_rank * 
                                                   diff_rank, rankOf_rankProduct = rank(rankProduct)) %>% 
    dplyr::arrange(rankOf_rankProduct) %>% dplyr::mutate(gene_category = dplyr::case_when(log2FoldChange > 
                                                                                            lfc_threshold & padj < padj_threshold ~ "up", log2FoldChange < 
                                                                                            -lfc_threshold & padj < padj_threshold ~ "down", TRUE ~ 
                                                                                            "static"), gene_category = factor(gene_category, levels = c("up", 
                                                                                                                                                        "down", "static")))
  upGenes_rank <- filter(merge_result, gene_category == "up")$RP_rank
  downGenes_rank <- filter(merge_result, gene_category == "down")$RP_rank
  staticGenes_rank <- filter(merge_result, gene_category == 
                               "static")$RP_rank
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
  down_static_pvalue <- suppressWarnings(ks.test(downGenes_rank, 
                                                 staticGenes_rank)$p.value)
  ks_test <- paste0("\n Kolmogorov-Smirnov Tests ", "\n pvalue of up vs static: ", 
                    format(up_static_pvalue, digits = 3, scientific = TRUE), 
                    "\n pvalue of down vs static: ", format(down_static_pvalue, 
                                                            digits = 3, scientific = TRUE))
  annotate_df <- data.frame(xpos = -Inf, ypos = Inf, annotateText = ks_test, 
                            hjustvar = 0, vjustvar = 1)
  p <- merge_result %>% ggplot2::ggplot(aes(x = RP_rank)) + 
    ggplot2::stat_ecdf(aes(color = gene_category), geom = "line") + 
    ggplot2::geom_text(data = annotate_df, aes(x = xpos, 
                                               y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText)) + 
    ggplot2::xlab("Regulatory potential rank") + ggplot2::ylab("Cumulative Probability")+
    ggplot2::scale_x_continuous(expand = c(0,0)) + ggplot2::theme_bw(base_size = 12)
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
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                df <- rbind(df, cnet1)
              }
            }
          }
        }
        if(length(df$ID) !=0){
          df$GeneRatio <- parse_ratio(df$GeneRatio)
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
          em <- enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- setReadable(em, org, 'ENTREZID')
              df[[name]] <- cnet1
            }
          }
        }
      }
      return(df)
    }
  }
}

enrich_genelist <- function(data, enrich_gene_list, showCategory=5){
  if(is.null(data) || is.null(enrich_gene_list)){
    return(NULL)
  }else{
    df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
    for (name in names(enrich_gene_list)) {
      sum <- length(data$ENTREZID[data$Group == name])
      em <- enrich_gene_list[[name]]
      if (length(as.data.frame(em)$ID) != 0) {
        if(length(colnames(as.data.frame(em))) == 9){
          cnet1 <- as.data.frame(em)
          cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
          cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
          if (length(cnet1$pvalue) > showCategory){
            cnet1 <- cnet1[1:showCategory,]
          }
          df <- rbind(df, cnet1)
        }
      }
    }
    if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
      p1 <- NULL
    } else{
      df$GeneRatio <- parse_ratio(df$GeneRatio)
      df <- dplyr::filter(df, !is.na(qvalue))
      df$Description <- gsub("_", " ", df$Description)
      df <- dplyr::mutate(df, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
      df$x <- gsub(":","", df$x)
      df <- dplyr::arrange(df, x)
      idx <- order(df[["x"]], decreasing = FALSE)
      df$Description <- factor(df$Description,
                               levels=rev(unique(df$Description[idx])))
      p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                      geom_point() +
                      scale_color_continuous(low="red", high="blue",
                                             guide=guide_colorbar(reverse=TRUE)) +
                      scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL)+
                      scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      p <- plot_grid(p1)
      return(p)
    }
  }
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
      return(data3) 
    }
  }
}
symbol2gene_id <- function(data,org){
  my.symbols <- rownames(data)
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL","ENTREZID"))
  colnames(gene_IDs) <- c("SYMBOL","gene_id")
  gene_IDs <- gene_IDs %>% distinct(SYMBOL, .keep_all = T) %>% na.omit()
  gene_IDs <- data.frame(gene_id = gene_IDs$gene_id, row.names = gene_IDs$SYMBOL)
  data <- merge(gene_IDs,data,by=0)
  data <- data[,-1]
  return(data)
}
