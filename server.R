popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=10000000*1024^2)
  output$Spe <- renderText({
    if(input$Species == "not selected" && input$Genomic_region == "Promoter") print("Please select 'Species'")
  })
  output$Spe_track <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_dist <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_dist_promoter <- renderText({
    if(input$Genomic_region == "Promoter") print("In the case of 'Promoter' mode, this function is not available. This function is for 'Genome-wide' mode. ")
  })
  output$Spe_motif <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_int <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  output$Spe_clustering <- renderText({
    if(input$Species_clustering == "not selected" && input$Genomic_region_clustering == "Promoter") print("Please select 'Species'")
  })
  output$Spe_clustering_track <- renderText({
    if(input$Species_clustering == "not selected") print("Please select 'Species'")
  })
  output$Spe_int_clustering <- renderText({
    if(input$Species_clustering == "not selected") print("Please select 'Species'")
  })
  output$Spe_venn_distribution <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  # pair-wise ------------------------------------------------------------------------------
  org1 <- reactive({
    return(org(Species = input$Species))
  })
  txdb <- reactive({
    if(input$Species != "not selected"){
      return(txdb_function(Species = input$Species))
    }
  })
  promoter_region <- reactive({
      if(input$Genomic_region == "Promoter"){
        if(input$Species != "not selected"){
        return(promoter(txdb(),upstream = input$upstream, downstream = input$downstream))
        }
      }else return(promoter(upstream = input$upstream, downstream = input$downstream,
                            input_type = "Genome-wide",files =peak_call_files()))
  })
  gene_position <- reactive({
    if(input$Species != "not selected"){
      return(genes(txdb()))
    }
  })
  
  bws <- reactive({
    if(input$data_file_type == "Row1"){
      if(is.null(input$file1)){
        if(input$goButton > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file1[, 1])){
          file <- input$file1[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$file1[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
        return(files)
      }
    }
  })
  bam <- reactive({
    if(input$data_file_type == "Row2"){
      if(is.null(input$file_bam)){
        if(input$goButton > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file_bam[, 1])){
          file <- input$file_bam[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$file_bam[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
        return(files)
      }
    }
  })
  
  peak_call_files <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      if(is.null(input$peak_call_file1)){
        if(input$goButton > 0 ){
          df <- list()
          df[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          df[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          df[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          df[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          name <- c("A_1","A_2","B_1","B_2")
          files2 <- lapply(df, GetGRanges, simple = TRUE)
          names(files2)<-name
          return(files2)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$peak_call_file1[, 1])){
          file <- input$peak_call_file1[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$peak_call_file1[nr,]$name))
          files <- c(files,file)
        }
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        names(files2)<-name
        print(files2)
        return(files2)
      }
    }
  })
  
  bw_files <- reactive({
    return(BigWigFileList(bws()))
  })
  
  output$input_bw_files <- DT::renderDataTable({
    if(input$data_file_type == "Row1"){
      uploaded_files = names(bw_files())
      as.data.frame(uploaded_files)
    }else{
      uploaded_files = names(bam())
      as.data.frame(uploaded_files) 
    }
  })
  
  output$input_peak_call_files <- DT::renderDataTable({
    if(input$Genomic_region == "Genome-wide"){
      uploaded_files = names(peak_call_files())
      as.data.frame(uploaded_files)
    }
  })
  
  bw_count <- reactive({
    if(input$data_file_type == "Row1" && !is.null(bw_files())){
      if(input$Genomic_region == "Promoter"){
        if(input$Species != "not selected"){
          return(Bigwig2count(bw = bw_files(),promoter_region(),
                              Species = input$Species,input_type =input$Genomic_region)) 
        }
      }else{
        if(!is.null(peak_call_files())){
          return(Bigwig2count(bw = bw_files(),promoter_region(),input_type =input$Genomic_region))
        }
      }
    }
    if(input$Species != "not selected" && input$data_file_type == "Row2" && 
       !is.null(bam())){
      withProgress(message = "converting bam to gene count data",{
        switch(input$Pair_or_single,
               "Paired-end" = seq <- TRUE,
               "Single-end" = seq <- FALSE)
        regionsToCount <- data.frame(GeneID = paste("ID", seqnames(promoter_region()), 
                                                    start(promoter_region()), end(promoter_region()), sep = "_"), Chr = seqnames(promoter_region()), 
                                     Start = start(promoter_region()), End = end(promoter_region()), Strand = strand(promoter_region()))
        a <- as.data.frame(promoter_region())
        Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
        if(input$Genomic_region == "Promoter"){
          count <- featureCounts(bam(),annot.ext = regionsToCount,isPairedEnd = seq,
                                 countMultiMappingReads =F,maxFragLength = 100)$counts
          rownames(count) <- Row.name
          colnames(count) <- names(bam())
          return(count)
        }else{
          if(!is.null(peak_call_files())){
            count <- featureCounts(bam(),annot.ext = regionsToCount,isPairedEnd = seq,
                                   countMultiMappingReads =F,maxFragLength = 100)$counts
            print(count)
            rownames(count) <- Row.name
            colnames(count) <- names(bam())
            return(count)
          }
        }
      })
    }
  })
  output$raw_count_table <- DT::renderDataTable({
    bw_count()
  })
  
  output$download_raw_count_table <- downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("bigwig2count-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("bigwig2count-genomeWide",".txt")
      }},
    content = function(file){write.table(bw_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  ##pair-wise deseq2--------
  dds <- reactive({
    count <- bw_count()
    collist <- gsub("\\_.+$", "", colnames(count))
    withProgress(message = "DESeq2",{
      group <- data.frame(con = factor(collist))
      dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
      dds$con <- factor(dds$con, levels = unique(collist))
      dds <- DESeq(dds)
      incProgress(1)
    })
    return(dds)
  })
  dds_limma <- reactive({
    withProgress(message = "Detecting differential accessible region (limma-trend)",{
    count <- log(bw_count() + 1,2)
    collist <- gsub("\\_.+$", "", colnames(count))
      group <- factor(collist)
      eset = new("ExpressionSet", exprs=as.matrix(count))
      design <- model.matrix(~0+collist)
      colnames(design) <- unique(collist)
      fit <- lmFit(eset, design)
      comparisons <-  paste(unique(collist)[1],"-",unique(collist)[2],sep="")
      cont.matrix <- makeContrasts(contrasts=comparisons, levels=design)
      fit <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit,trend = TRUE)
      result =topTable(fit2, number = 1e12)
      colnames(result) <- c("log2FoldChange","baseMean","t","pval","padj","B")
      incProgress(1)
    })
    return(result)
  })
  
  
  deg_result <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      if(input$data_file_type == "Row2"){
      count <- bw_count()
      collist <- gsub("\\_.+$", "", colnames(count))
      dds <- dds()
      contrast <- c("con", unique(collist))
      res <- results(dds,  contrast = contrast)
      res <- as.data.frame(res)
      }
      if(input$data_file_type == "Row1"){
        res <- dds_limma()
      }
      return(res)
    }
  })
  
  deg_norm_count <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      count <- bw_count()
      if(input$data_file_type == "Row2"){
      collist <- gsub("\\_.+$", "", colnames(count))
      group <- data.frame(con = factor(collist))
      dds <- dds()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
      }
      if(input$data_file_type == "Row1"){
        normalized_counts <- count
      }
      return(normalized_counts)
    }
  })
  
  observeEvent(bws(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="raw_count_panel")
  }))
  observeEvent(bam(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="raw_count_panel")
  }))
  observeEvent(peak_call_files(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="peak_call_files_panel")
  }))
  
  data_degcount <- reactive({
    data <- deg_result()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      Row.names <- NULL
      log2FoldChange <- NULL
      value <- NULL
      data <- merge(data,count, by=0)
      data <- dplyr::filter(data, apply(data[,8:(7 + Cond_1 + Cond_2)],1,mean) > input$basemean)
      data$log2FoldChange <- -1 * data$log2FoldChange
      if(input$Genomic_region == "Promoter"){
        if(input$Species != "not selected"){
          my.symbols <- data$Row.names
          gene_IDs<-id_convert(my.symbols,input$Species,type="SYMBOL_double")
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(data, gene_IDs, by="Row.names")
        }
      }
      genenames <- as.vector(data$Row.names)
      return(data)
    }
  })
  
  data_degcount2 <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- dplyr::filter(data, abs(data$log2FoldChange) > log(input$fc, 2))
      if(nrow(data2) != 0){
        data2$group <- "Up"
        data2$group[data2$log2FoldChange < 0] <- "Down"
        data3 <- dplyr::filter(data2, abs(data2$padj) < input$fdr)
        return(data3)
      }else{return(NULL)}
    }
  })
  
  deg_result_anno <- reactive({
    if(!is.null(deg_result())){
      withProgress(message = "preparing annotation",{
        data <- promoter_region()
        data2 <- annotatePeak(data, TxDb= txdb())
        return(data2)
      })
    }
  })
  deg_result_anno2 <- reactive({
    data <- as.data.frame(as.GRanges(deg_result_anno()))
    Row.name <- paste0(data$seqnames,":",data$start,"-",data$end)
    data$locus <- Row.name
    my.symbols <- data$geneId
    gene_IDs<-id_convert(my.symbols,input$Species,type="ENTREZID")
    colnames(gene_IDs) <- c("geneId","NearestGene")
    data <- merge(data, gene_IDs, by="geneId")
    data <- data %>% distinct(locus, .keep_all = T)
    rownames(data)<-data$locus
    return(data)
  })
  
  deg_result2 <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      anno <- deg_result_anno2()
      anno <- dplyr::arrange(anno, locus)
      data <- data.frame(NearestGene = anno$NearestGene,
                         locus = anno$locus)
      deg_result <- deg_result()
      deg_result$locus <- rownames(deg_result())
      deg_result <- dplyr::arrange(deg_result, locus)
      data2 <- merge(data,deg_result,by="locus")
      rownames(data2) <- data2$locus
      data2 <- data2[, - which(colnames(data2) == "locus")]
      return(data2)
    }
  })
  
  output$DEG_result <- DT::renderDT({
    if(!is.null(deg_result())){
      if(input$Genomic_region == "Promoter" || input$Species == "not selected"){
        deg_result() %>%
          datatable(
            selection = "single",
            filter = "top")
      }else{
        deg_result2()  %>%
          datatable(
            selection = "single",
            filter = "top")
      }
    }
  })
  output$Normalized_Count_matrix <- DT::renderDataTable({
    deg_norm_count()
  })
  output$download_pair_DEG_result = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("DESeq2_result-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("DESeq2_result-genomeWide",".txt")
      }},
    content = function(file){
      if(input$Genomic_region == "Promoter"){
        write.table(deg_result(), file, row.names = T, sep = "\t", quote = F)
      }else write.table(deg_result2(), file, row.names = T, sep = "\t", quote = F)
    }
  )
  output$download_pair_norm_count = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("DESeq2_normCount-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("DESeq2_normCount-genomeWide", ".txt")
      }},
    content = function(file){write.table(deg_norm_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  # pair-wise PCA ------------------------------------------------------------------------------
  output$download_pair_PCA = downloadHandler(
    filename = "PCA-MDS-dendrogram.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(PCAplot(data = deg_norm_count(),plot=TRUE))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA <- renderPlot({
    withProgress(message = "preparing result table",{
    if(is.null(deg_norm_count())){
      return(NULL)
    }else{
      print(PCAplot(data = deg_norm_count(),plot=TRUE))
    }
    })
  })
  
  output$pair_PCA_data <- DT::renderDataTable({
    PCAplot(data = deg_norm_count(),plot=FALSE)
  })
  
  output$download_pair_PCA_table = downloadHandler(
    filename ="PCA_table.txt",
    content = function(file){write.table( PCAplot(data = deg_norm_count(),plot=FALSE), 
                                          file, row.names = T, sep = "\t", quote = F)}
  )
  ###GOI----------------------
  output$volcano_x <- renderUI({
      if(!is.null(data_degcount())){
      withProgress(message = "Preparing volcano plot",{
        data <- as.data.frame(data_degcount())
        min <- floor(min(data$log2FoldChange))
        max <- ceiling(max(data$log2FoldChange))
    sliderInput("xrange","X_axis range:",min = min-1,
                max=max+1, step = 1,
                value = c(min, max))
    })
      }
  })
  output$volcano_y <- renderUI({
    if(!is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
    max <- ceiling(max(-log10(data$padj)))
    sliderInput("yrange","Y_axis range:",min = 0, max= max+1, step = 1,
                value = max)
    }
  })
  
  pair_volcano <- reactive({
    if(!is.null(input$xrange) && !is.null(input$yrange)){
      withProgress(message = "volcano plot",{
        data <- as.data.frame(data_degcount())
        count <- deg_norm_count()
        if(!is.null(input$DEG_result_rows_selected)){
          if(input$Genomic_region == "Promoter"){
            label_data <- rownames(deg_result()[input$DEG_result_rows_selected,])
          }else label_data <- rownames(deg_result2()[input$DEG_result_rows_selected,])
        }else label_data <- NULL
        data$color <- "NS"
        data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- paste("down:", length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]))
        data$color[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- paste("up:",length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))
        data$padj[data$padj == 0] <- 10^(-300)
        if(!is.null(label_data)) {
          for(name in label_data){
            data$color[data$Row.names == name] <- "GOI"
          }
          Color <- c("blue","green","darkgray","red")
          data$color <- factor(data$color, levels = c(paste("down:", length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr])),
                                                      "GOI","NS", paste("up:",length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))))
          if(length(data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]) == 0) Color <- c("green","darkgray","red")
        }else{
          Color <- c("blue","darkgray","red")
          data$color <- factor(data$color, levels = c(paste("down:", length(data$Row.names[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr])),
                                                      "NS", paste("up:",length(data$Row.names[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr]))))
          if(length(data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr]) == 0) Color <- c("darkgray","red")
        }
        
        v <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = color),size = 0.4)
        v <- v  + geom_vline(xintercept = c(-log2(input$fc), log2(input$fc)), linetype = c(2, 2), color = c("black", "black")) +
          geom_hline(yintercept = c(-log10(input$fdr)), linetype = 2, color = c("black"))
        v <- v +theme_bw()+ scale_color_manual(values = Color)+
          theme(legend.position = "top" , legend.title = element_blank(),
                axis.text.x= ggplot2::element_text(size = 12),
                axis.text.y= ggplot2::element_text(size = 12),
                text = ggplot2::element_text(size = 12),
                title = ggplot2::element_text(size = 12)) +
          xlab("log2 fold change") + ylab("-log10(padj)") +
          xlim(input$xrange)+
          ylim(c(0, input$yrange))
        if(!is.null(label_data)) {
          v <- v + geom_point(data=dplyr::filter(data, Row.names == label_data),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, Row.names == label_data), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), fontface = "bold.italic")
        }
        return(v)
        incProgress(1)
      })
    }else return(NULL)
  })
  
  output$volcano1 <- renderPlot({
    withProgress(message = "Preparing results",{
    if(!is.null(input$xrange)){
      if(is.null(bw_count())){
        return(NULL)
      }else{
        
        print(pair_volcano())
      }
    }
    })
  })
  
  pair_heatmap <- reactive({
    withProgress(message = "Heatmap",{
      count <- deg_norm_count()
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- data_degcount2()
      if(is.null(data2)){
        ht <- NULL
      }else{
        data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
        ht <- GOIheatmap(data.z,show_row_names = FALSE)
      }
      return(ht)
      incProgress(1)
    })
  })
  output$GOIheatmap <- renderPlot({
    pair_heatmap()
  })
  volcano_heatmap <- reactive({
    vol <- as.grob(pair_volcano())
    heat <- as.grob(pair_heatmap())
    return(plot_grid(vol,heat, rel_widths = c(2, 1)))
  })
  
  output$download_pair_volcano = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("Volcano-promoter(", -input$upstream,"-",input$downstream,").pdf")
      }else{
        paste0("Volcano-genomeWide",".pdf")
      }},
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 4
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_volcano())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_pair_heatmap = downloadHandler(
    filename = function() {
      if (input$Genomic_region=='Promoter'){
        paste0("Heatmap-promoter(", -input$upstream,"-",input$downstream,").pdf")
      }else{
        paste0("Heatmap-genomeWide",".pdf")
      }},
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 3
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_heatmap())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##trackplot----------
  ref <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species))
    return(ref)
  })
  
  
  
  observeEvent(input$DEG_result_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_pair_DEG", open="Trackplot_panel")
  }))
  observeEvent(input$DEG_result_rows_selected,({
    if(input$Species != "not selected"){
    if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
      updateSliderInput(session,"igv_uprange","Range:",
                        value = c(start_position,end_position),
                        step = 100, min = start_position - 10000, max = end_position + 10000)
    }
    }
  }))
  
  output$igv_uprange <- renderUI({
    if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("igv_uprange","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_ylim <- renderUI({
    numericInput("igv_ylim","peak range:", value = 1, min = 0)
  })
  
  gtrack <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=0.8)
      return(gtrack)
    })
  })
  
  goi_promoter_position<- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      library(Gviz)
      if(input$Genomic_region == "Promoter"){
        label_data <- rownames(deg_result()[input$DEG_result_rows_selected,])
        gene_IDs<-id_convert(label_data,input$Species,type="SYMBOL_single")
        y <- as.data.frame(subset(promoter_region(), gene_id %in% gene_IDs))
      }else{
        label_data <- rownames(deg_result2()[input$DEG_result_rows_selected,])
        y <- dplyr::filter(deg_result_anno2(),locus == label_data)
      }
      return(y)
    }
  })
  
  goi_gene_position <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      if(input$Genomic_region == "Promoter"){
        label_data <- rownames(deg_result()[input$DEG_result_rows_selected,])
      }else{
        label_data <- deg_result2()[input$DEG_result_rows_selected,]$NearestGene
      }
      gene_IDs<- id_convert(label_data,input$Species,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  output$trackplot_additional <- renderUI({
    fileInput("trackplot_additional1",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  track_additional_files <-reactive({
    if(!is.null(input$trackplot_additional1)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$trackplot_additional1[, 1])){
        file <- input$trackplot_additional1[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$trackplot_additional1[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  data_track <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      return(data_trac(y=goi_promoter_position(),gene_position=goi_gene_position(),
                       gen=ref(),txdb=txdb(),org=org1(),filetype=input$data_file_type,
                       bw_files=bws(),bam_files=bam(),
                       track_additional_files=track_additional_files()))
    }
  })
  
  highlight_trackplot <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      library(Gviz)
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      chr <- gene_position$seqnames
      df <- data_track()
      if(y$start < input$igv_uprange[2] && y$end > input$igv_uprange[1]){
        withProgress(message = "Highlighting the selected region",{
          ht <- HighlightTrack(trackList = df,
                               start = y$start, width = y$width,
                               chromosome = chr, alpha=0.5)
        })
      }else ht <- NULL
      return(ht)
    }
  })
  goi_trackplot <- reactive({
    if(!is.null(input$DEG_result_rows_selected) &&
       !is.null(goi_gene_position()) && 
       !is.null(input$igv_uprange)){
      library(Gviz)
      if(!is.null(highlight_trackplot())){
        plot<- plotTracks(list(gtrack(), highlight_trackplot()),
                          from = input$igv_uprange[1], 
                          to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                          type="hist")
      }else{
        df <- data_track()
        df[["gtrack"]] <- gtrack()
        plot<- plotTracks(df,
                          from = input$igv_uprange[1], 
                          to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                          type="hist")
      }
      return(plot)
    }
  })
  output$trackplot_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$DEG_result_rows_selected)  &&
         !input$Species == "not selected"){
        goi_trackplot()
      }
    })
  })
  
  output$download_pair_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(highlight_trackplot())){
          plotTracks(list(gtrack(), highlight_trackplot()),
                            from = input$igv_uprange[1], 
                            to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                            type="hist")
        }else{
          df <- data_track()
          df[["gtrack"]] <- gtrack()
          plotTracks(df,
                            from = input$igv_uprange[1], 
                            to = input$igv_uprange[2],ylim=c(0,input$igv_ylim),
                            type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  
  ##Functional enrichment analysis-------------
  output$Gene_set <- renderUI({
    selectInput('Gene_set', 'Gene Set', gene_set_list)
  })
  data_degcount_up <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- data_degcount2()
      up_all <- dplyr::filter(data2, log2FoldChange > 0)
      rownames(up_all) <- up_all$Row.names
      up_all <- up_all[,8:(7 + Cond_1 + Cond_2)]
      return(up_all)
    }
  })
  data_degcount_up_bed <-reactive({
    if(input$Genomic_region == "Genome-wide"){
      data <- range_changer(data_degcount_up())
      data2 <- data.frame(chr = data$chr, start = data$start, end = data$end)
    }else{
      up <- symbol2gene_id(data_degcount_up(),org1())
      up2 <- subset(promoter_region(), gene_id %in% up$gene_id) %>% as.data.frame()
      data2 <- data.frame(chr = up2$seqnames, start = up2$start, end = up2$end)
    }
    return(data2)
  })
  data_degcount_down_bed <-reactive({
    if(input$Genomic_region == "Genome-wide"){
      data <- range_changer(data_degcount_down())
      data2 <- data.frame(chr = data$chr, start = data$start, end = data$end)
    }else{
      down <- symbol2gene_id(data_degcount_down(),org1())
      down2 <- subset(promoter_region(), gene_id %in% down$gene_id) %>% as.data.frame()
      data2 <- data.frame(chr = down2$seqnames, start = down2$start, end = down2$end)
    }
    return(data2)
  })
  output$DEG_up_result <- DT::renderDT({
    data_degcount_up_bed()
  })
  output$DEG_down_result <- DT::renderDT({
    data_degcount_down_bed()
  })
  output$download_pair_DEG_up_result = downloadHandler(
    filename = function() {
      paste0("Up-DAR",".bed")
    },
    content = function(file){write.table(data_degcount_up_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  output$download_pair_DEG_down_result = downloadHandler(
    filename = function() {
      paste0("Down-DAR",".bed")
    },
    content = function(file){write.table(data_degcount_down_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  
  data_degcount_down <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- data_degcount2()
      down_all <- dplyr::filter(data2, log2FoldChange < 0)
      rownames(down_all) <- down_all$Row.names
      down_all <- down_all[,8:(7 + Cond_1 + Cond_2)]
      return(down_all)
    }
  })
  Hallmark_set <- reactive({
    if(!is.null(input$Gene_set)){
      return(GeneList_for_enrichment(Species = input$Species, Gene_set = input$Gene_set, org = org1()))
    }
  })
  
  data_degcount2_anno <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      data <- data_degcount2()
      rownames(data) <- data$Row.names
      data <- data[, -1]
      anno <- deg_result_anno2()
      colnames(anno)[1] <- "ENTREZID"
      data2 <- merge(data,anno,by=0)
      return(data2)
    }
  })
  
  enrichment_enricher <- reactive({
    data3 <- data_degcount2()
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(data3)){
      if(input$Genomic_region == "Promoter"){
        withProgress(message = "enrichment analysis",{
          H_t2g <- Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          em_up <- try(enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
          em_down <- try(enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
          df <- list()
          df[["Up"]] <- em_up
          df[["Down"]] <- em_down
          for(name in names(df)){
            if (length(as.data.frame(df[[name]])$ID) == 0) {
              df[[name]] <- NULL
            } else{
              df[[name]] <- setReadable(df[[name]], org1(), 'ENTREZID')
            }
          }
          incProgress(1)
          return(df)
        })
      }else{ 
        data3 <- deg_result_anno2()
        H_t2g <- Hallmark_set()
        source <- ref_for_GREAT(input$Species)
        data_degcount_up2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_up()))
        data_degcount_up3 <- with(data_degcount_up2, GRanges(seqnames = seqnames, 
                                                             ranges = IRanges(start,end),
                                                             ENTREZID = geneId))
        res_up = rGREAT::great(gr = data_degcount_up3,gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species),source)
        data_degcount_down2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_down()))
        data_degcount_down3 <- with(data_degcount_down2, GRanges(seqnames = seqnames, 
                                                                 ranges = IRanges(start,end),
                                                                 ENTREZID = geneId))
        res_down = rGREAT::great(data_degcount_down3,gene_list_for_enrichment_genome(H_t2g,input$Species), source)
        df <- list()
        df[["Up"]] <- res_up
        df[["Down"]] <- res_down
        incProgress(1)
        return(df)
      }
    }else{return(NULL)}
  })
  
  enrichment_1_1 <- reactive({
    df <- enrichment_enricher()
    if(input$Genomic_region == "Promoter"){
      if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
        data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
        for(name in names(df)){
          if(!is.null(df[[name]])) {
            group1 <- as.data.frame(df[[name]])
          }else group1 <- NULL
          group1$Group <- name
          data <- rbind(data, group1)
        }
        if(length(data$Description) != 0) {
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- parse_ratio(data$GeneRatio)
          return(data)
        }else return(NULL)
      }else{return(NULL)}
    }else{
      if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
        data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
        colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                            "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
        for(name in names(df)){
          if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
            group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
            group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
            group1$Group <- name
          }else group1 <- NULL
          data <- rbind(data, group1)
        }
        if(length(data$id) != 0) {
          return(data)
        }else return(NULL)
      }else{return(NULL)}
    }
  })
  
  pair_enrich1_H <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      count <- deg_norm_count()
      df <- enrichment_enricher()
      if(input$Genomic_region == "Promoter"){
        data3 <- data_degcount2()
        H_t2g <- Hallmark_set()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        if (is.null(df[["Up"]]) && is.null(df[["Down"]]))  {
          p1 <- NULL
        } else{
          data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          for(name in names(df)){
            if(!is.null(df[[name]])) {
              group1 <- as.data.frame(df[[name]])
              group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$ENTREZID),")",sep = "")
              if (length(group1$pvalue) > 5){
                group1 <- group1[sort(group1$pvalue, decreasing = F, index=T)$ix,]
                group1 <- group1[1:5,]
              }
            }else group1 <- NULL
            data <- rbind(data, group1)
          }
          if(length(data$Description) != 0){
            data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
            data$GeneRatio <- parse_ratio(data$GeneRatio)
            if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
              p1 <- NULL
            } else{
              data$Description <- gsub("_", " ", data$Description)
              data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
              data$x <- gsub(":","", data$x)
              data <- dplyr::arrange(data, x)
              idx <- order(data[["x"]], decreasing = FALSE)
              data$Description <- factor(data$Description,
                                         levels=rev(unique(data$Description[idx])))
              p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="qvalue",size="GeneRatio"))+
                              geom_point() +
                              scale_color_continuous(low="red", high="blue",
                                                     guide=guide_colorbar(reverse=TRUE)) +
                              scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                              scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
            }}else p1 <- NULL
        }
      }else{
        data3 <- data_degcount2_anno()
        if (is.null(df[["Up"]]) && is.null(df[["Down"]]))  {
          p1 <- NULL
        } else{
          data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
          for(name in names(df)){
            if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
              group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
              group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$group),")",sep = "")
              if (length(group1$p_adjust_hyper) > 5){
                group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
                group1 <- group1[1:5,]
              }
            }else group1 <- NULL
            data <- rbind(data, group1)
          }
          colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                                 "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
          if(length(data$Description) != 0){
            data$Description <- gsub("_", " ", data$Description)
            data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                            scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
          }else p1 <- NULL
        }
      }
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment1 <- renderPlot({
    dotplot_for_output(data = bw_count(),
                       plot_genelist = pair_enrich1_H(), Gene_set = input$Gene_set, 
                       Species = input$Species)
  })
  
  output$download_pair_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H()
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  pair_enrich2 <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      if(input$Genomic_region == "Promoter"){
        data3 <- data_degcount2()
      }else data3 <- data_degcount2_anno()
      count <- deg_norm_count()
      df <- enrichment_enricher()
      upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
      downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
      p <- list()
      for(name in names(df)){
        if(length(as.data.frame(df[[name]])$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- setReadable(df[[name]], org1(), 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          if(name == "Up") genes <- upgene
          if(name == "Down") genes <- downgene
          geneList <- genes$log2FoldChange
          names(geneList) = as.character(genes$ENTREZID)
          p2 <- try(as.grob(cnetplot(cnet1, foldChange=geneList,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }
        }
        p[[name]] <- p2
      }
      p <- plot_grid(p[["Up"]], p[["Down"]], nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  
  pair_enrich_table <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1()), H_t2g = Hallmark_set(), Gene_set = input$Gene_set))
  })
  
  output$pair_enrichment_result <- DT::renderDataTable({
    as.data.frame(enrichment_1_1())
  })
  
  output$download_pair_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$whichGeneSet <- renderUI({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Up_or_Down) && !is.null(as.data.frame(enrichment_1_1())$id)){
      group <- input$Up_or_Down
      set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1()), Group == group)$id)
      selectInput('Pathway_list', 'Pathway list', set_list)
    }
  })
  
  region_gene_associate <- reactive({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list)){
      set_list <- input$Pathway_list
      df <- enrichment_enricher()
      data3 <- data_degcount2()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0 &&
           length(as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))$seqnames) != 0) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$group),")",sep = "")
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_associations <- DT::renderDT({
    if(input$Genomic_region == "Genome-wide"){
      region_gene_associate()
    }
  })
  output$whichGenes <- renderUI({
    if(input$Genomic_region == "Genome-wide"){
      if(!is.null(enrichment_enricher())){
        selectInput('Up_or_Down', ' Up or Down', names(enrichment_enricher()))
      }
    }
  })
  region_gene_associate_plot <- reactive({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list) && 
       !is.null(input$Up_or_Down) && !is.null(enrichment_enricher())){
      set_list <- input$Pathway_list
      df <- enrichment_enricher()
      name <- input$Up_or_Down
      if(length(as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))$id) != 0) {
        res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
      }else res <- NULL
      return(res)
    }
  })
  output$region_gene_associations_plot <- renderPlot({
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Pathway_list) && 
       !is.null(input$Up_or_Down) && !is.null(enrichment_enricher())){
    if(class(try(region_gene_associate_plot())) !="try-error") region_gene_associate_plot()
    }
  })
  output$download_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$Pathway_list,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$pair_pdf_width
        set_list <- input$Pathway_list
        df <- enrichment_enricher()
        name <- input$Up_or_Down
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations",input$Pathway_list,".txt")
    },
    content = function(file){write.table(region_gene_associate(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ##Motif enrichment analysis----------
  output$homer_showCategory <- renderUI({
    sliderInput("homer_showCategory","Most significant motifs", value=5, min=1,max=20)
  })
  output$homer_size <- renderUI({
    radioButtons('homer_size','Type of the region for motif finding',
                 c('given (exact size)'="given",
                   'custom size'="custom"
                 ),selected = "custom")
  })
  output$homer_bg <- renderUI({
    if(input$Genomic_region == "Genome-wide"){
    radioButtons('homer_bg','Background sequence',
                 c('random'="random",
                   'peak call files'="peakcalling"
                 ),selected = "random")
    }
  })
  output$homer_size2 <- renderUI({
    if(!is.null(input$homer_size)){
    if(input$homer_size == "custom"){
    numericInput('homer_size2','Size of the region for motif finding',value=200, step=100)}}
  })
  output$homer_unknown <- renderUI({
    selectInput("homer_unknown","Type of enrichment analysis",c("known motif","known and de novo motifs"), selected = "known motif")
  })
  observeEvent(input$homer_unknown,({
    updateActionButton(session,"motifButton", "Start")
  }))
  
  preMotif_list <- reactive({
    df <- list()
    df[["Up"]] <- data_degcount_up()
    df[["Down"]] <- data_degcount_down()
    return(df)
  })
  
  enrich_motif <- reactive({
    if(input$homer_size == "given") size <- "given"
    if(input$homer_size == "custom") size <- input$homer_size2
    if(input$motifButton > 0 && !is.null(preMotif_list()) 
       && input$Species != "not selected" && !is.null(input$homer_unknown)){
      if(input$Genomic_region == "Genome-wide"){
        return(findMotif(df= preMotif_list(), anno_data = deg_result_anno2(),back = input$homer_bg,
                             Species = input$Species, motif=input$homer_unknown,size=size,bw_count=bw_count()))
      }else return(findMotif(df= data_degcount2(), anno_data = promoter_region(),
                                 Species = input$Species,
                                 type = "Promoter", motif=input$homer_unknown,size=size))
    }
  })
  
  output$motif_plot <- renderPlot({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && 
       !is.null(input$homer_unknown) && input$Species != "not selected"){
      homer_Motifplot(df = enrich_motif(),showCategory = input$homer_showCategory)
    }
  })
  output$download_motif_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = enrich_motif(),showCategory = input$homer_showCategory)
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_homer_report = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  
  output$download_motif_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  motif_table <- reactive({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      return(known_motif(enrich_motif()))
    }
  })
  denovo_motif_table <- reactive({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      return(denovo_motif(enrich_motif()))
    }
  })
  
  output$motif_result <- DT::renderDT({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && input$Species != "not selected"){
      motif_table()
    }
  })
  output$denovo_motif_result <- DT::renderDT({
    if(input$motifButton > 0 && !is.null(enrich_motif()) && input$Species != "not selected"){
      denovo_motif_table()
    }
  })
  
  
  output$download_denovo_motif_table = downloadHandler(
    filename = function() {
      paste0("denovo_motif_table",".txt")
    },
    content = function(file){write.table(denovo_motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ##peak distribution-------
  input_peak_list <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      files<-list()
      name <- c()
      if(!is.null(peak_call_files())){
        a <- peak_call_files()
        Glist <- GRangesList()
        for(name in names(a)){
          Glist[[name]] <- a[[name]]
        }
        Glist[["Consensus_region"]] <- promoter_region()
        return(Glist)
      }else{
        if(input$goButton > 0 ){
          files[["data/peakcall/A_1_peaks.narrowPeak"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          files[["data/peakcall/A_2_peaks.narrowPeak"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          files[["data/peakcall/B_1_peaks.narrowPeak"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          files[["data/peakcall/B_2_peaks.narrowPeak"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          Glist <- files2GRangelist(files)
          Glist[["Consensus_region"]] <- promoter_region()
          return(Glist)
        }else return(NULL)
      }
    }
  })
  output$input_peak_distribution <- renderPlot({
    withProgress(message = "Plot peak distribution for peak call files",{
      if(!is.null(input_peak_list()) && !is.null(txdb()) && 
         input$Genomic_region == "Genome-wide"){
        genomicElementDistribution(deg_peak_list()[["Up"]], 
                                   TxDb = txdb(),) 
      }
    })
  })
  deg_peak_list <- reactive({
    if(!is.null(data_degcount_up()) && !is.null(data_degcount_down())){
      up <- range_changer(data_degcount_up())
      down <- range_changer(data_degcount_down())
      Glist <- GRangesList()
      up2 <- with(up, GRanges(seqnames = chr,ranges = IRanges(start=start,
                                                              end=end)))
      down2 <- with(down, GRanges(seqnames = chr,ranges = IRanges(start=start,
                                                                  end=end)))
      Glist[["Up"]] <- up2
      Glist[["Down"]] <- down2
      return(Glist)
    }
  })  
  output$deg_peak_distribution <- renderPlot({
    withProgress(message = "Plot peak distribution for DAR",{
      if(input$Genomic_region == "Genome-wide"){
      if(!is.null(deg_peak_list()) && input$Species != "not selected"){
        genomicElementDistribution(deg_peak_list()[["Down"]], 
                                   TxDb = txdb())
      }
      }
    })
  })
  
  
  output$download_input_peak_distribution = downloadHandler(
    filename = function() {
      paste0("Up_DAR_peak_distribution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        genomicElementDistribution(deg_peak_list()[["Up"]], 
                                   TxDb = txdb()) 
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_deg_peak_distribution = downloadHandler(
    filename = function() {
      paste0("Down_DAR_peak_distribution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        genomicElementDistribution(deg_peak_list()[["Down"]], 
                                   TxDb = txdb()) 
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ###Peak pattern comparison up--------
  output$peak_pattern_up_additional <- renderUI({
    fileInput("peak_pattern_up_add",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  up_additional <-reactive({
    if(!is.null(input$peak_pattern_up_add)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$peak_pattern_up_add[, 1])){
        file <- input$peak_pattern_up_add[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$peak_pattern_up_add[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  output$peak_pattern_up_heat_range <- renderUI({
    numericInput("peak_up_range","Intensity range",value=10,min = 0.1,step=2)
  })
  peak_up_grange <- reactive({
    if(!is.null(data_degcount_down())){
    if(input$Genomic_region == "Genome-wide"){
      up <- range_changer(data_degcount_up())
      up <- with(up, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
    }else{
      up <- symbol2gene_id(data_degcount_up(),org1()) %>% distinct(gene_id, .keep_all = T)
      up <- subset(promoter_region(), gene_id %in% up$gene_id) 
    }
    return(up)
    }else return(NULL)
  })
  peak_up_cvglists <-reactive({
    up2 <- peak_up_grange()
    feature.recentered <- reCenterPeaks(up2, width=4000)
    files <- bws()
    if(!is.null(up_additional())) files <- c(files, up_additional())
    cvglists <- sapply(files, import,which=feature.recentered,as="RleList")
    names(cvglists) <- gsub(".+\\/","",gsub("\\..+$", "", names(files)))
    return(cvglists)
  })
  
  peak_up_sig <- reactive({
    up2 <- peak_up_grange()
    cvglists <- peak_up_cvglists()
    feature.center <- reCenterPeaks(up2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000)
    return(sig)
  })
  
  peak_up_feature_center <- reactive({
    up2 <- peak_up_grange()
    feature.center <- reCenterPeaks(up2, width=1)
    return(feature.center)
  })
  peak_up_alinedHeatmap <- reactive({
    if(!is.null(input$peak_up_range)){
      range <- c()
      for(n in length(names(bws()))) {range <- c(range, input$peak_up_range)}
      heatmap <- featureAlignedHeatmap(peak_up_sig(), peak_up_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = brewer.pal(n = 9, name = 'OrRd'))
      return(heatmap)
    }
  })
  output$peak_pattern_up_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(deg_result())){
        plot_grid(peak_up_alinedHeatmap())
      }
    })
  })
  output$peak_pattern_up_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(deg_result())){
        featureAlignedDistribution(peak_up_sig(), peak_up_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
      }
    })
  })
  
  output$download_peak_pattern_up_heatmap = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_up_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_up_alinedHeatmap()))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_up_line = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_up_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        featureAlignedDistribution(peak_up_sig(), peak_up_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
        dev.off()
        incProgress(1)
      })
    }
  )
  ###Peak pattern comparison down--------
  output$peak_pattern_down_additional <- renderUI({
    fileInput("peak_pattern_down_add",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  output$peak_pattern_down_heat_range <- renderUI({
    numericInput("peak_down_range","Intensity range",value=10,min = 0.1,step=2)
  })
  peak_down_grange <- reactive({
    if(!is.null(data_degcount_down())){
    if(input$Genomic_region == "Genome-wide"){
      down <- range_changer(data_degcount_down())
      down <- with(down, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
    }else{
      down <- symbol2gene_id(data_degcount_down(),org1()) %>% distinct(gene_id, .keep_all = T)
      down <- subset(promoter_region(), gene_id %in% down$gene_id)
    }
    return(down)
    }else return(NULL)
  })
  peak_down_cvglists <-reactive({
    down2 <- peak_down_grange()
    feature.recentered <- reCenterPeaks(down2, width=4000)
    files <- bws()
    if(!is.null(up_additional())) files <- c(files, up_additional())
    cvglists <- sapply(files, import,which=feature.recentered,as="RleList")
    names(cvglists) <- gsub(".+\\/","",gsub("\\..+$", "", names(files)))
    return(cvglists)
  })
  
  peak_down_sig <- reactive({
    down2 <- peak_down_grange()
    cvglists <- peak_down_cvglists()
    feature.center <- reCenterPeaks(down2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000)
    return(sig)
  })
  
  peak_down_feature_center <- reactive({
    down2 <- peak_down_grange()
    feature.center <- reCenterPeaks(down2, width=1)
    return(feature.center)
  })
  peak_down_alinedHeatmap <- reactive({
    if(!is.null(input$peak_down_range)){
      range <- c()
      for(n in length(names(bws()))) {range <- c(range, input$peak_down_range)}
      heatmap <- featureAlignedHeatmap(peak_down_sig(), peak_down_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = brewer.pal(n = 9, name = 'OrRd'))
      return(heatmap)
    }
  })
  output$peak_pattern_down_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(deg_result())){
        plot_grid(peak_down_alinedHeatmap())
      }
    })
  })
  output$peak_pattern_down_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(deg_result())){
        featureAlignedDistribution(peak_down_sig(), peak_down_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
      }
    })
  })
  output$download_peak_pattern_down_heatmap = downloadHandler(
    filename = function(){
      paste0("Peak_pattern_down_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_down_alinedHeatmap()))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_down_line = downloadHandler(
    filename = function(){
      paste0("Peak_pattern_down_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        featureAlignedDistribution(peak_down_sig(), peak_down_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##Regulatory potential------
  output$pairRNAseqresult <- renderUI({
    fileInput("pair_DEG_result",
              "Select RNA-seq DEG result file",
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  RNAseqDEG <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      tmp <- input$pair_DEG_result$datapath
      if(is.null(input$pair_DEG_result) && input$goButton > 0 )  tmp = "data/RNAseq.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
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
  })
  
  output$pair_DEG_result <- renderDataTable({
    if(!is.null(RNAseqDEG())){
      RNAseqDEG() 
    }
  })
  
  RNAseqDEG_anno <- reactive({
    RNAdata <- RNAseqDEG()
    RNAdata$log2FoldChange <- -RNAdata$log2FoldChange
    if(str_detect(rownames(RNAdata)[1], "ENS")){
      my.symbols <- gsub("\\..*","", rownames(RNAdata))
      gene_IDs<-id_convert(my.symbols,input$Species,type="ENSEMBL")
      colnames(gene_IDs) <- c("EnsemblID","Symbol","gene_id")
      RNAdata$EnsemblID <- gsub("\\..*","", rownames(RNAdata))
      gene_IDs <- gene_IDs %>% distinct(EnsemblID, .keep_all = T)
    }else{
      my.symbols <- rownames(RNAdata)
      gene_IDs<-id_convert(my.symbols, input$Species,type="SYMBOL_single")
      colnames(gene_IDs) <- c("Symbol", "gene_id")
      gene_IDs <- gene_IDs %>% distinct(Symbol, .keep_all = T)
      RNAdata$Symbol <- rownames(RNAdata) 
    }
    data <- merge(RNAdata, gene_IDs, by="Symbol")
    return(data)
  })
  
  
  
  mmAnno_up <- reactive({
    up_peak <- peak_up_grange()
    if(!is.null(up_peak) && !is.null(input$peak_distance)){
    if(input$Genomic_region == "Promoter"){ 
      up_peak <- up_peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      up_peak <- with(up_peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
    }
    range <- input$peak_distance * 1000
    mcols(up_peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(up_peak))))
    if(length(as.data.frame(up_peak)$seqnames) != 0){
    mmAnno_up <- mm_geneScan(up_peak, txdb(),upstream = range,downstream = range)
    }else mmAnno_up <- NULL
    return(mmAnno_up)
    }else return(NULL)
  })
  mmAnno_down <- reactive({
    down_peak <- peak_down_grange()
    if(!is.null(down_peak) && !is.null(input$peak_distance)){
    if(input$Genomic_region == "Promoter"){ 
      down_peak <- down_peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      down_peak <- with(down_peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
    }
    range <- input$peak_distance * 1000
    mcols(down_peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(down_peak))))
    if(length(as.data.frame(down_peak)$seqnames) != 0){
    mmAnno_down <- mm_geneScan(down_peak, txdb(),upstream = range,downstream = range)
    }else mmAnno_down <- NULL
    return(mmAnno_down)
    }else return(NULL)
  })
  
  RP <- reactive({
    mmAnno_up <- mmAnno_up()
    if(!is.null(mmAnno_up)) {
      result_geneRP_up <- calcRP_TFHit(mmAnno = mmAnno_up,Txdb = txdb())
      colnames(result_geneRP_up)[2:4] <- paste0(colnames(result_geneRP_up)[2:4], "_up")
    }else result_geneRP_up <- NULL
    mmAnno_down <- mmAnno_down()
    if(!is.null(mmAnno_down)) {
      result_geneRP_down <- calcRP_TFHit(mmAnno = mmAnno_down,Txdb = txdb())
      colnames(result_geneRP_down)[2:4] <- paste0(colnames(result_geneRP_down)[2:4], "_down")
    }else result_geneRP_down <- NULL
    if(!is.null(result_geneRP_up) && !is.null(result_geneRP_down)){
      result_geneRP <- merge(result_geneRP_up,result_geneRP_down,by="gene_id",all =T) 
      result_geneRP$sumRP <- apply(data.frame(result_geneRP$sumRP_up,
                                              -result_geneRP$sumRP_down),1,sum,na.rm=TRUE)
    }
    if(!is.null(result_geneRP_up) && is.null(result_geneRP_down)) {
      result_geneRP <-result_geneRP_up
      result_geneRP$sumRP <- result_geneRP$sumRP_up
    }
    if(is.null(result_geneRP_up) && !is.null(result_geneRP_down)) {
      result_geneRP <-result_geneRP_down
      result_geneRP$sumRP <- -result_geneRP$sumRP_down
    }
    if(is.null(result_geneRP_up) && is.null(result_geneRP_down)){
      return(NULL)
    }else{
    result_geneRP <- result_geneRP %>% dplyr::arrange(-sumRP)
    result_geneRP$RP_rank <- rownames(result_geneRP) %>% as.numeric()
    return(result_geneRP)
    }
  })
  regulatory_potential <- reactive({
    if(input$Species != "not selected"){
    data <- RNAseqDEG_anno()
    result_geneRP <- RP()
    merge_data <- integrate_ChIP_RNA(
      result_geneRP = result_geneRP,
      result_geneDiff = data,lfc_threshold = log(input$DEG_fc,2),padj_threshold = input$DEG_fdr
    )
    return(merge_data)
    }
  })
 
  output$ks_plot <- renderPlot({    
    if(!is.null(input$peak_distance) && !is.null(RNAseqDEG()) && 
       !is.null(regulatory_potential()) && !is.na(input$DEG_fdr) && 
       !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      regulatory_potential()
    }
  })
  output$DEG_fc <- renderUI({
    numericInput("DEG_fc","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr <- renderUI({
    numericInput("DEG_fdr","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  output$peak_distance <- renderUI({
    sliderInput("peak_distance","Regulatory range (distance (kb) from TSS)",
                min=0,max=200,value=100,step = 5)
  })
  output$RNAseqGroup <- renderUI({
    if(input$Species != "not selected" &&!is.null(mmAnno_up()) && !is.null(mmAnno_down()) && !is.null(input$peak_distance)){
      if(!is.null(RP_all_table())){
      selectInput("RNAseqGroup","Group (RNAseq-Epigenome)",
                  unique(RP_all_table()$Group),
                  multiple = FALSE)
      }
    }
  })
  pari_RNAseq_boxplot <- reactive({
    RNA <- RNAseqDEG_anno()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    data <- merge(RNA,RP(), by="gene_id",all=T)
    data$group <- "Others"
    data$group[data$sumRP > input$DEG_fdr] <- "RP > 1"
    data$group[data$sumRP < -input$DEG_fdr] <- "RP < -1"
    data$group <- factor(data$group,levels=c("Others","RP > 1","RP < -1"),ordered=TRUE)
    
    collist <- unique(data$group)
    if (length(collist) >= 3){
      stat.test <- data %>% tukey_hsd(log2FoldChange ~ group)
      stat.test <- stat.test %>% add_significance("p.adj")
      stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      col <- c("gray","#F8766D","#00BFC4")
    }else{
      group1 <- dplyr::filter(data, group == collist[1])
      group2 <- dplyr::filter(data, group == collist[2])
      if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
      stat.test <- data %>% t_test(log2FoldChange ~ group)
      stat.test <- stat.test %>% add_significance()
      stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      col <-c("gray","#F8766D")
      }else stat.test <- NULL
    }
    if(!is.null(stat.test)){
    p <- try(ggpubr::ggboxplot(data, x = "group", y = "log2FoldChange",
                           fill = "group", scales = "free", add = "jitter",
                           xlab = FALSE, ylab = "log2FoldChange")+theme_bw(base_size = 15)+
      xlab(NULL)+ylab("RNAseq log2FoldChange")+scale_fill_manual(values = col) + stat_pvalue_manual(stat.test,hide.ns = T,size = 5))
    }
    return(p)
  })
  output$RNAseq_boxplot_error <- renderText({
    if(!is.null(input$peak_distance) && !is.null(RNAseqDEG()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      RNA <- RNAseqDEG_anno()
      RNA <- dplyr::filter(RNA, !is.na(gene_id))
      data <- merge(RNA,RP(), by="gene_id",all=T)
      data$group <- "Others"
      data$group[data$sumRP > input$DEG_fdr] <- "RP > 1"
      data$group[data$sumRP < -input$DEG_fdr] <- "RP < -1"
      data$group <- factor(data$group,levels=c("Others","RP > 1","RP < -1"),ordered=TRUE)
      collist <- unique(data$group)
      if (length(collist) < 3){
      group1 <- dplyr::filter(data, group == collist[1])
      group2 <- dplyr::filter(data, group == collist[2])
      if(length(rownames(group1)) <= 1 || length(rownames(group2)) <= 1){
        print("boxplot: There are few genes with |RP| > 1")
      }
      }
    }
  })
  RNAseq_popu <- reactive({
    RNA <- RNAseqDEG_anno()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$group[RNA$log2FoldChange > log2(input$DEG_fc) & RNA$padj < input$DEG_fdr] <- "up"
    RNA$group[RNA$log2FoldChange < -log2(input$DEG_fc) & RNA$padj < input$DEG_fdr] <- "down"
    data <- merge(RNA,RP(), by="gene_id",all=T)
    if(length(which(unique(data$group) == "up")) == 1){
    up <- data %>% dplyr::filter(group == "up")
    up_total <- length(rownames(up))
    if(length(which(up$sumRP > 0))  > 0){
    up_epiUp <- length(rownames(dplyr::filter(up, sumRP > 0)))
    }else up_epiUp <- 0
    if(length(which(up$sumRP < 0))  > 0){
    up_epiDown <- length(rownames(dplyr::filter(up, sumRP < 0)))
    }else up_epiDown <- 0
    up_per <- c(up_total-(up_epiUp + up_epiDown), up_epiDown, up_epiUp)
    }else up_per <- c(0,0,0)
    if(length(which(unique(data$group) == "down")) == 1){
    down <- data %>% dplyr::filter(group == "down")
    down_total <- length(rownames(down))
    if(length(which(down$sumRP > 0))  > 0){
    down_epiUp <- length(rownames(dplyr::filter(down, sumRP > 0)))
    }else down_epiUp <- 0
    if(length(which(down$sumRP < 0))  > 0){
    down_epiDown <- length(rownames(dplyr::filter(down, sumRP < 0)))
    }else down_epiDown <- 0
    down_per <- c(down_total-(down_epiUp + down_epiDown), down_epiDown, down_epiUp)
    }else down_per <- c(0,0,0)
    x <- data.frame(RNA = c("up", "up","up","down","down","down"), 
                    Regulatory_potential = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"),
                                                    paste0("down:NS (",down_per[1],")"),paste0("down:down (",down_per[2],")"),paste0("down:up (",down_per[3],")")),
                                                  levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"),
                                                             paste0("down:NS (",down_per[1],")"),paste0("down:down (",down_per[2],")"),paste0("down:up (",down_per[3],")"))),
                    num = c(up_per,down_per),
                    col = c("F8766D","#00BFC4","grey","F8766D","#00BFC4","grey"))
    p <- ggplot(x, aes(x = RNA, y = num, fill = Regulatory_potential)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#00BFC4","#F8766D","grey","#00BFC4","#F8766D"))
    return(p)
  })
  output$int_boxplot <- renderPlot({
    if(!is.null(input$peak_distance) && !is.null(RNAseqDEG()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      pari_RNAseq_boxplot()
    }
  })
  
  output$int_bar <- renderPlot({
    if(!is.null(input$peak_distance) && !is.null(RNAseqDEG()) && 
       !is.na(input$DEG_fdr) && !is.na(input$DEG_fc) && !is.null(bw_count()) && input$Species != "not selected"){
      gridExtra::grid.arrange(RNAseq_popu(), ChIPseq_popu(), ncol = 1)
    }
  })
  ChIPseq_popu <- reactive({
    RNA <- RNAseqDEG_anno()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$gene_category <- "NS"
    RNA$gene_category[RNA$log2FoldChange > log2(input$DEG_fc) & RNA$padj < input$DEG_fdr] <- "up"
    RNA$gene_category[RNA$log2FoldChange < -log2(input$DEG_fc) & RNA$padj < input$DEG_fdr] <- "down"
    data <- merge(RNA,RP(), by="gene_id")
    data$epigenome_category <- "up"
    data$epigenome_category[data$sumRP < 0] <- "down"
    if(length(which(unique(data$epigenome_category) == "up")) == 1){
      up <- data %>% dplyr::filter(epigenome_category == "up")
      up_total <- length(rownames(up))
      if(length(which(unique(up$gene_category) == "up")) == 1){
        up_RNAUp <- length(rownames(dplyr::filter(up, gene_category == "up")))
      }else up_RNAUp <- 0
      if(length(which(unique(up$gene_category) == "down")) == 1){
      up_RNADown <- length(rownames(dplyr::filter(up, gene_category == "down")))
      }else up_RNADown <- 0
      up_per <- c(up_total-(up_RNAUp + up_RNADown), up_RNADown, up_RNAUp)
    }else{
      up_per <- c(0,0,0)
    } 
    if(length(which(unique(data$epigenome_category) == "down")) == 1){
    down <- data %>% dplyr::filter(epigenome_category == "down")
    down_total <- length(rownames(down))
    if(length(which(unique(down$gene_category) == "up")) == 1){
    down_RNAUp <- length(rownames(dplyr::filter(down, gene_category == "up")))
    }else down_RNAUp <- 0
    if(length(which(unique(down$gene_category) == "down")) == 1){
    down_RNADown <- length(rownames(dplyr::filter(down, gene_category == "down")))
    }else down_RNADown <- 0
    down_per <- c(down_total-(down_RNAUp + down_RNADown), down_RNADown, down_RNAUp)
    }else{
      down_per <- c(0,0,0)
    } 
    x <- data.frame(Regulatory_potential = c("up", "up","up","down","down","down"), 
                    RNA = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"),
                                   paste0("down:NS (",down_per[1],")"),paste0("down:down (",down_per[2],")"),paste0("down:up (",down_per[3],")")),
                                 levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"),
                                            paste0("down:NS (",down_per[1],")"),paste0("down:down (",down_per[2],")"),paste0("down:up (",down_per[3],")"))),
                    num = c(up_per,down_per),
                    col = c("#F8766D","#00BFC4","grey","#F8766D","#00BFC4","grey"))
    p <- ggplot(x, aes(x = Regulatory_potential, y = num, fill = RNA)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#00BFC4","#F8766D","grey","#00BFC4","#F8766D"))
    return(p)
  })
  RP_all_table <- reactive({
    target_result <- regulatory_potential()$data
    target_result$epigenome_category <- "up"
    target_result$epigenome_category[target_result$sumRP < 0] <- "down"
    table <- NULL
    if(!is.null(mmAnno_up()) && !is.null(mmAnno_down())) {
    table <- data.frame(Symbol = target_result$Symbol,
                        Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                        RNA_log2FC = -target_result$log2FoldChange,
                        RNA_padj = target_result$padj,
                        regulatory_potential = target_result$sumRP,
                        withUpPeakN = target_result$withPeakN_up,
                        withDownPeakN = target_result$withPeakN_down,
                        gene_id = target_result$gene_id)
    }else{
      if(!is.null(mmAnno_up())){
        table <- data.frame(Symbol = target_result$Symbol,
                            Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                            RNA_log2FC = -target_result$log2FoldChange,
                            RNA_padj = target_result$padj,
                            regulatory_potential = target_result$sumRP,
                            withUpPeakN = target_result$withPeakN_up,
                            gene_id = target_result$gene_id)
      }
      if(!is.null(mmAnno_down())){
        table <- data.frame(Symbol = target_result$Symbol,
                            Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                            RNA_log2FC = -target_result$log2FoldChange,
                            RNA_padj = target_result$padj,
                            regulatory_potential = target_result$sumRP,
                            withDownPeakN = target_result$withPeakN_down,
                            gene_id = target_result$gene_id)
      }
    }
    return(table)
  })
  
  RP_selected_table <- reactive({
    table <- RP_all_table() %>% dplyr::filter(Group == input$RNAseqGroup)
    return(table)
  })
  
  output$RP_table <- renderDT({
    if(!is.null(input$RNAseqGroup) && !is.null(input$peak_distance)){
      RP_selected_table() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$download_pairintbox = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling_boxplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pari_RNAseq_boxplot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_pairintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu(), ChIPseq_popu(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_pairKSplot = downloadHandler(
    filename = function(){
      paste0("KSplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(regulatory_potential())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_RP_table = downloadHandler(
    filename = function() {
      paste0("RP_summary_table.txt")
    },
    content = function(file){write.table(RP_all_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_selected_RP_table = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup,"_table.txt")
    },
    content = function(file){write.table(RP_selected_table(), file, row.names = F, sep = "\t", quote = F)}
  )
 
  int_selected_bed <- reactive({
    if(!is.null(RP_selected_table())){
      gene <- RP_selected_table()$gene_id
      type <- gsub(".+\\-","", input$RNAseqGroup)
      if(type == "up") {
        peak <- subset(mmAnno_up(), gene_id %in% gene)
        up_peak2 <- as.data.frame(peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      if(type == "down") {
        peak <- subset(mmAnno_down(), gene_id %in% gene)
        down_peak2 <- as.data.frame(peak)
        down_peak2$Row.names <- paste0(down_peak2$seqnames,":",down_peak2$start,"-",down_peak2$end)
        down_peak2 <- down_peak2 %>% distinct(Row.names, .keep_all = T)
        down_peak3 <- with(down_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(down_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(down_peak3)
      }
      return(y)
    }
  })
  output$download_selected_int_peak = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup,"_table.bed")
    },
    content = function(file){write.table(int_selected_bed(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  observeEvent(input$pair_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_pair_RP", open="KS_panel")
  }))
  
  
  
  ##Integrative trackplot-------
  observeEvent(input$RP_table_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_pair_RP", open="int_Trackplot_panel")
  }))
  output$int_igv_uprange <- renderUI({
    if(!is.null(int_goi_gene_position()) && !is.null(int_goi_promoter_position())){
      y <- int_goi_promoter_position()
      gene_position <- int_goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim <- renderUI({
    numericInput("int_igv_ylim","peak range:", value = 2, min = 0)
  })
  
  int_goi_promoter_position<- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table()[input$RP_table_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(mmAnno_up())) {
        up_peak <- subset(mmAnno_up(), gene_id %in% gene)
        mcols(up_peak) <- DataFrame(Group = "up")
      }
      if(!is.null(mmAnno_down())) {
        down_peak <- subset(mmAnno_down(), gene_id %in% gene)
        mcols(down_peak) <- DataFrame(Group = "up")
        y <- as.data.frame(down_peak)
      }
      if(!is.null(mmAnno_up()) && !is.null(mmAnno_down())) {
        peak <- c(up_peak,down_peak)
        y <- as.data.frame(peak)
      }else{
        if(!is.null(mmAnno_up())) y <- as.data.frame(up_peak)
        if(!is.null(mmAnno_down())) y <- as.data.frame(down_peak)
      }
      return(y)
    }
  })
  
  int_goi_gene_position <- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      gene <- RP_selected_table()[input$RP_table_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  output$int_trackplot_additional <- renderUI({
    fileInput("int_trackplot_additional1",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  int_track_additional_files <-reactive({
    if(!is.null(input$int_trackplot_additional1)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$int_trackplot_additional1[, 1])){
        file <- input$int_trackplot_additional1[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$int_trackplot_additional1[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  int_data_track <- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      return(data_trac(y=int_goi_promoter_position(),gene_position=int_goi_gene_position(),
                       gen=ref(),txdb=txdb(),org=org1(),filetype=input$data_file_type,
                       bw_files=bws(),bam_files=bam(),
                       track_additional_files=int_track_additional_files()))
    }
  })
  
  int_highlight_trackplot <- reactive({
    if(!is.null(input$RP_table_rows_selected) && !is.null(input$int_igv_uprange)){
      library(Gviz)
      y <- int_goi_promoter_position()
      gene_position <- int_goi_gene_position()
      chr <- gene_position$seqnames
      df <- int_data_track()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange[2] && y[i,]$end > input$int_igv_uprange[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          if(y[i,]$Group == "up") {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == "down") {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  int_goi_trackplot <- reactive({
    if(!is.null(input$RP_table_rows_selected) &&
       !is.null(int_goi_promoter_position()) && 
       !is.null(int_goi_gene_position()) && 
       !is.null(gtrack()) &&
       !is.null(input$int_igv_uprange)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot())){
        plot<- plotTracks(list(gtrack(), int_highlight_trackplot()),
                          from = input$int_igv_uprange[1], 
                          to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                          type="hist")
      }else{
        df <- int_data_track()
        df[["gtrack"]] <- gtrack()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange[1], 
                          to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_rows_selected)){
        int_goi_trackplot()
      }
    })
  })
  output$download_pair_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot())){
          plotTracks(list(gtrack(), int_highlight_trackplot()),
                            from = input$int_igv_uprange[1], 
                            to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                            type="hist")
        }else{
          df <- int_data_track()
          df[["gtrack"]] <- gtrack()
          plotTracks(df,
                            from = input$int_igv_uprange[1], 
                            to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                            type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #Integrative functional analysis--------
  int_Hallmark_set <- reactive({
    if(!is.null(input$intGeneset)){
      return(GeneList_for_enrichment(Species = input$Species, Gene_set = input$intGeneset, org = org1()))
    }
  })
  
  output$intGroup <- renderUI({
    if(!is.null(RP_all_table())){
      selectInput("intGroup","Group (RNAseq-Epigenome)",unique(RP_all_table()$Group),multiple = T)
    }
  })
  output$intGeneset <- renderUI({
    selectInput('intGeneset', 'Gene Set', gene_set_list)
  })
  
  selected_int_group <- reactive({
    group <- input$intGroup
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      table <- RP_all_table() %>% dplyr::filter(Group == name)
      head(table)
      df2 <- data.frame(ENTREZID = table$gene_id, Group = table$Group)
      df <- rbind(df,df2)
    }
    return(df)
  })
  int_enrich <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group(), Species = input$Species, org = org1(),
                                   H_t2g = int_Hallmark_set(),Gene_set = input$intGeneset))
  })
  int_enrich_list <- reactive({
    return(enrich_gene_list(data = selected_int_group(),
                            Gene_set = input$intGeneset, org = org1(), H_t2g = int_Hallmark_set()))
  })
  int_enrich_plot <- reactive({
    return(enrich_genelist(data = selected_int_group(),
                           enrich_gene_list = int_enrich_list()))
  })
  output$int_enrichment1 <- renderPlot({
    dotplot_for_output(data = int_enrich(),
                       plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                       Species = input$Species)
  })
  int_enrich_table <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich()), H_t2g = int_Hallmark_set(), Gene_set = input$intGeneset))
  })
  output$int_enrichment_result <- DT::renderDataTable({
    int_enrich_table()
  })
  output$download_pair_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich(),
                           plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                           Species = input$Species)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_pair_int_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset,".txt")
    },
    content = function(file){write.table(int_enrich_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  ##Venn diagram -----------
  output$Spe_motif_venn <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT_venn <- renderText({
    if(input$Species_venn == "not selected") print("Please select 'Species'")
  })
  
  txdb_venn <- reactive({
    if(input$Species_venn != "not selected"){
      return(txdb_function(Species = input$Species_venn))
    }
  })
  org_venn <- reactive({
    return(org(Species = input$Species_venn))
  })
  gene_position_venn <- reactive({
    if(input$Species_venn != "not selected"){
      return(genes(txdb_venn()))
    }
  })
  #venn-------
  Venn_peak_call_files <- reactive({
    if(is.null(input$peak_call_file_venn1)){
      if(input$goButton_venn > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$peak_call_file_venn1[, 1])){
        file <- input$peak_call_file_venn1[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$peak_call_file_venn1[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  
  bws_venn <- reactive({
    if(is.null(input$file_venn1)){
      if(input$goButton_venn > 0 ){
        df<-list()
        df[["A_1"]] <- "data/bigwig/A_1.BigWig"
        df[["A_2"]] <- "data/bigwig/A_2.BigWig"
        df[["B_1"]] <- "data/bigwig/B_1.BigWig"
        df[["B_2"]] <- "data/bigwig/B_2.BigWig"
        return(df)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$file_venn1[, 1])){
        file <- input$file_venn1[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$file_venn1[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  venn_overlap <- reactive({
    withProgress(message = "Preparing intersection, takes a few minutes",{
    ol <- findOverlapsOfPeaks(Venn_peak_call_files())
    names(ol$peaklist) <- gsub("///","-",names(ol$peaklist))
    return(ol)
    })
  })
  make_venn <-reactive({
    df<-list()
    for(name in names(Venn_peak_call_files())){
      data <- as.data.frame(Venn_peak_call_files()[[name]])
      data$strand <- as.character(data$strand)
      df[[name]] <- data
    }
    return(df)
  })

  output$venn <- renderPlot({
    withProgress(message = "Venn diagram",{
      if(!is.null(Venn_peak_call_files())){
        
ggVennPeaks(make_venn(),label_size = 5, alpha = .2)
      }
    })
  })
  
  output$download_vennplot = downloadHandler(
    filename = function(){
      paste0("Venn_diagram",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 7
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(ggVennPeaks(make_venn(),label_size = 5, alpha = .2))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$select_file2 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("intersect_select", "Select intersect",
                  c("not_selected",names(venn_overlap()$peaklist)),
                  selected = "not_selected",multiple = F)
    }
  })
  
  selected_grange <- reactive({
    if(input$intersect_select != "not_selected"){
      return(venn_overlap()$peaklist[[input$intersect_select]])
    }
  })
  
  output$venn_peak_distribution <- renderPlot({
    withProgress(message = "Peak distribution",{
      if(!is.null(input$intersect_select)){
      if(input$intersect_select != "not_selected" && !is.null(selected_grange()) &&
         input$Species_venn != "not selected"){
        genomicElementDistribution(selected_grange(), 
                                   TxDb = txdb_venn())
      }
      }
    })
  })
  
  selected_annoData_table <- reactive({
    withProgress(message = "Preparing annotation",{
      overlaps.anno <- annotatePeak(genomicElementDistribution(selected_grange(), 
                                                               TxDb = txdb_venn(),
                                                               plot = F)$peaks,
                                    TxDb = txdb_venn()) %>% as.data.frame()
      my.symbols <- overlaps.anno$geneId
      gene_IDs<-id_convert(my.symbols,input$Species_venn,type="ENTREZID")
      colnames(gene_IDs) <- c("geneId","NearestGene")
      data <- merge(gene_IDs,overlaps.anno,by="geneId")
      data <- data[,2:11] %>% distinct(peakNames, .keep_all = T)
      return(data)
    })
  })
  output$selected_intersect_annotation <- DT::renderDT({
    if(input$intersect_select != "not_selected" && !is.null(selected_grange()) &&
       input$Species_venn != "not selected" ){
      selected_annoData_table() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  
  output$download_venn_peak_distribution = downloadHandler(
    filename = function(){
      paste0("Peak_distribution-",input$intersect_select,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4.5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        genomicElementDistribution(selected_grange(), 
                                   TxDb = txdb_venn())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_selected_intersect_annotation_table <- downloadHandler(
    filename = function() {
      paste0("Intersect_table-",input$intersect_select,".txt")
    },
    content = function(file){write.table(as.data.frame(selected_annoData_table()), 
                                         file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_selected_intersect_annotation_table_bed <- downloadHandler(
    filename = function() {
      paste0("Intersect_table-",input$intersect_select,".bed")
    },
    content = function(file){write.table(as.data.frame(selected_annoData_table()[,2:4]), 
                                         file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  
  ###Venn Peak pattern comparison--------
  output$peak_pattern_venn_heat_range <- renderUI({
    numericInput("peak_venn_range","Intensity range",value=10,min = 0.1,step=2)
  })
  
  peak_venn_grange <- reactive({
    return(selected_grange())
  })
  
  peak_venn_cvglists <-reactive({
    venn2 <- peak_venn_grange()
    feature.recentered <- reCenterPeaks(venn2, width=4000)
    files <- bws_venn()
    cvglists <- sapply(files, import,which=feature.recentered,as="RleList")
    names(cvglists) <- gsub(".+\\/","",gsub("\\..+$", "", names(files)))
    return(cvglists)
  })
  
  peak_venn_sig <- reactive({
    venn2 <- peak_venn_grange()
    cvglists <- peak_venn_cvglists()
    feature.center <- reCenterPeaks(venn2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000,)
    return(sig)
  })
  
  peak_venn_feature_center <- reactive({
    venn2 <- peak_venn_grange()
    feature.center <- reCenterPeaks(venn2, width=1)
    return(feature.center)
  })
  peak_venn_alinedHeatmap <- reactive({
    if(!is.null(input$peak_venn_range)){
      range <- c()
      for(n in length(names(bws_venn()))) {range <- c(range, input$peak_venn_range)}
      heatmap <- featureAlignedHeatmap(peak_venn_sig(), peak_venn_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = brewer.pal(n = 9, name = 'OrRd'))
      return(heatmap)
    }
  })
  output$peak_pattern_venn_heatmap <- renderPlot({
    if(!is.null(input$intersect_select) && 
       input$intersect_select != "not_selected" && !is.null(bws_venn())){
      withProgress(message = "feature aligned heatmap",{
        plot_grid(peak_venn_alinedHeatmap())
      })
    }
  })
  output$peak_pattern_venn_line <- renderPlot({
    if(!is.null(input$intersect_select) && input$intersect_select != "not_selected" && !is.null(bws_venn())){
      withProgress(message = "feature aligned distribution",{
        featureAlignedDistribution(peak_venn_sig(), peak_venn_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
      })
    }
  })
  output$download_peak_pattern_venn_heatmap = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_venn_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_venn_alinedHeatmap()))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_venn_line = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_venn_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        featureAlignedDistribution(peak_venn_sig(), peak_venn_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##venn_trackplot----------
  ref_venn <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species_venn))
    return(ref)
  })
  
  
  
  observeEvent(input$selected_intersect_annotation_rows_selected, ({
    updateCollapse(session,id =  "int_result_collapse_panel", open="Trackplot_venn_panel")
  }))
  observeEvent(input$selected_intersect_annotation_rows_selected,({
    if(input$Species_venn != "not selected"){
    if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
      y <- goi_promoter_position_venn()
      gene_position <- goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
      updateSliderInput(session,"igv_venn_uprange","Range:",value = c(start_position,end_position),
                        step = 100, min = start_position - 10000, max = end_position + 10000)
    }
    }
  }))
  
  output$igv_venn_uprange <- renderUI({
    if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
      y <- goi_promoter_position_venn()
      gene_position <- goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("igv_venn_uprange","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_venn_ylim <- renderUI({
    numericInput("igv_venn_ylim","peak range:", value = 2, min = 0)
  })
  
  gtrack_venn <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=0.8)
      return(gtrack)
    })
  })
  
  goi_promoter_position_venn<- reactive({
    if(!is.null(input$selected_intersect_annotation_rows_selected)){
      library(Gviz)
      label_data <- selected_annoData_table()[input$selected_intersect_annotation_rows_selected,]
      return(label_data)
    }
  })
  
  goi_gene_position_venn <- reactive({
    if(!is.null(goi_promoter_position_venn())){
      label_data <- goi_promoter_position_venn()$NearestGene
      gene_IDs<- id_convert(label_data,input$Species_venn,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position_venn(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  
  data_track_venn <-reactive({
    return(data_trac(y=goi_promoter_position_venn(),gene_position=goi_gene_position_venn(),
                     gen=ref_venn(),txdb=txdb_venn(),org=org_venn(),
                     bw_files=bws_venn(),track_additional_files=NULL))
  })
  
  highlight_trackplot_venn <- reactive({
    library(Gviz)
    y <- goi_promoter_position_venn()
    gene_position <- goi_gene_position_venn()
    chr <- gene_position$seqnames
    df <- data_track_venn()
    if(y$start < input$igv_venn_uprange[2] && y$end > input$igv_venn_uprange[1]){
      withProgress(message = "Highlighting the selected region",{
        ht <- HighlightTrack(trackList = df,
                             start = y$start, width = y$width,
                             chromosome = chr, alpha=0.5)
      })
    }else ht <- NULL
    return(ht)
  })
  goi_trackplot_venn <- reactive({
    withProgress(message = "Trackplot",{
      if(!is.null(input$selected_intersect_annotation_rows_selected) &&
         !is.null(input$igv_venn_uprange)){
        library(Gviz)
        if(!is.null(highlight_trackplot_venn())){
          plot<- plotTracks(list(gtrack_venn(), highlight_trackplot_venn()),
                            from = input$igv_venn_uprange[1], 
                            to = input$igv_venn_uprange[2],ylim=c(0,input$igv_venn_ylim),
                            type="hist")
        }else{
          df <- data_track_venn()
          df[["gtrack"]] <- gtrack_venn()
          plot<- plotTracks(df,
                            from = input$igv_venn_uprange[1], 
                            to = input$igv_venn_uprange[2],ylim=c(0,input$igv_venn_ylim),
                            type="hist") 
        }
        return(plot)
      }
    })
  })
  output$trackplot_venn_goi <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      goi_trackplot_venn()
    })
  })
  
  output$download_venn_trackplot = downloadHandler(
    filename = "Trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        plotTracks(list(gtrack_venn(), highlight_trackplot_venn()),
                   from = goi_gene_position_venn()$start-input$igv_venn_uprange, 
                   to = goi_gene_position_venn()$end+input$igv_venn_downrange,ylim=c(0,input$igv_venn_ylim),
                   type="hist")
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##Venn motif -----------
  output$homer_showCategory_venn <- renderUI({
    sliderInput("homer_showCategory_venn","Most significant motifs", value=5, min=1,max=20)
  })
  output$homer_size_venn <- renderUI({
    radioButtons('homer_size_venn','Type of the region for motif finding',
                 c('given (exact size)'="given",
                   'custom size'="custom"
                 ),selected = "custom")
  })
  output$homer_bg_venn <- renderUI({
      radioButtons('homer_bg_venn','Background sequence',
                   c('random'="random",
                     'bed files'="peakcalling"
                   ),selected = "random")
  })
  output$homer_bg2_venn <- renderUI({
    if(!is.null(input$homer_bg_venn)){
      if(input$homer_bg_venn == "peakcalling"){
        fileInput('homer_bg2_venn',
                  'Select bed files',
                  accept = c("bed","narrowPeak"),
                  multiple = TRUE,
                  width = "80%")
        }}
  })
  Venn_peak_call_files_homer <- reactive({
    if(is.null(input$homer_bg2_venn)){
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$homer_bg2_venn[, 1])){
        file <- input$homer_bg2_venn[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$homer_bg2_venn[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      files3 <- soGGi:::runConsensusRegions(GRangesList(files2), "none")
      return(files3)
    }
  })
  
  output$homer_size2_venn <- renderUI({
    if(!is.null(input$homer_size_venn)){
    if(input$homer_size_venn == "custom"){
      numericInput('homer_size2_venn','Size of the region for motif finding',value=200, step=100)}}
  })
  output$homer_unknown_venn <- renderUI({
    selectInput("homer_unknown_venn","Type of enrichment analysis",c("known motif","known and de novo motifs"), selected = "known motif")
  })
  observeEvent(input$homer_unknown_venn,({
    updateActionButton(session,"motifButton_venn", "Start")
  }))
  
  
  output$venn_whichGroup1 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      selectInput("venn_whichGroup1", "Select intersects", choices = c(names(venn_overlap()$peaklist)),multiple = TRUE)
    }
  })
  preMotif_list_venn <- reactive({
    df <- list()
    for(name in input$venn_whichGroup1){
      df[[name]] <- venn_overlap()$peaklist[[name]]
    }
    return(df)
  })

  enrich_motif_venn <- reactive({
    if(input$homer_size_venn == "given") size <- "given"
    if(input$homer_size_venn == "custom") size <- input$homer_size2_venn
    if(input$motifButton_venn > 0 && !is.null(preMotif_list_venn()) 
       && input$Species_venn != "not selected" && !is.null(input$homer_unknown_venn)){
      return(findMotif(df= preMotif_list_venn(),Species = input$Species_venn,size=size,back = input$homer_bg_venn,
                       motif=input$homer_unknown_venn, other_data = Venn_peak_call_files_homer(),type="Other"))
    }
  })
  
  output$motif_venn_plot <- renderPlot({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && 
       !is.null(input$homer_unknown_venn) && input$Species_venn != "not selected"){
      homer_Motifplot(df = enrich_motif_venn(),showCategory = input$homer_showCategory_venn)
    }
  })
  
  output$download_motif_venn_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = enrich_motif_venn(),showCategory = input$homer_showCategory_venn)
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  output$download_motif_venn_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )

  motif_table_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      return(known_motif(enrich_motif_venn()))
    }
  })
  denovo_motif_table_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      return(denovo_motif(enrich_motif_venn()))
    }
  })
  
  output$denovo_motif_venn_result <- DT::renderDT({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && input$Species_venn != "not selected"){
      denovo_motif_table_venn()
    }
  })
  
  output$motif_venn_result <- DT::renderDT({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn()) && input$Species_venn != "not selected"){
      motif_table_venn()
    }
  })
  
  
  output$download_denovo_motif_venn_table = downloadHandler(
    filename = function() {
      paste0("denovo_motif_table",".txt")
    },
    content = function(file){write.table(denovo_motif_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_homer_report_venn = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif_venn()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  #Venn functional analysis------
  output$Gene_set_venn <- renderUI({
    selectInput('Gene_set_venn', 'Gene Set', gene_set_list)
  })
  
  output$venn_whichGroup2 <- renderUI({
    if(!is.null(Venn_peak_call_files())){
      if(!is.null(venn_overlap())){
        selectInput("venn_whichGroup2", "Select intersects", choices = c(names(venn_overlap()$peaklist)),multiple = TRUE)
      }
    }
  })
  
  Hallmark_set_venn <- reactive({
    if(!is.null(input$Gene_set_venn)){
      return(GeneList_for_enrichment(Species = input$Species_venn, Gene_set = input$Gene_set_venn, org = org_venn()))
    }
  })
  
  selected_grange_venn_list <- reactive({
    if(!is.null(input$venn_whichGroup2)){
      Glist <- GRangesList()
      for(name in input$venn_whichGroup2){
        data <- venn_overlap()$peaklist[[name]]
        Glist[[name]] <- data
      }
      return(Glist)
    }
  })
  
  enrichment_enricher_venn <- reactive({
    data <- selected_grange_venn_list()
    H_t2g <- Hallmark_set_venn()
    df <- list()
    source <- ref_for_GREAT(input$Species_venn)
    for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species_venn),source)
      df[[name]] <- res
    }
    return(df)
  })
  
  enrichment_1_1_venn <- reactive({
    df <- enrichment_enricher_venn()
    if(!is.null(input$Gene_set_venn) && input$Species_venn != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                          "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
        }else group1 <- NULL
        group1$Group <- name
        data <- rbind(data, group1)
      }
      if(length(data$id) != 0) {
        return(data)
      }else return(NULL)
    }else{return(NULL)}
  })
  
  pair_enrich1_H_venn <- reactive({
    if(!is.null(input$Gene_set_venn) && input$Species_venn != "not selected"){
      df <- enrichment_enricher_venn()
      data3 <- selected_grange_venn_list()
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
          if (length(group1$p_adjust_hyper) > 5){
            group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
            group1 <- group1[1:5,]
          }
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                             "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      if(length(data$Description) != 0){
        data$Description <- gsub("_", " ", data$Description)
        data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
        data$x <- gsub("-","", data$x)
        data <- dplyr::arrange(data, x)
        idx <- order(data[["x"]], decreasing = FALSE)
        data$Description <- factor(data$Description,
                                   levels=rev(unique(data$Description[idx])))
        p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                        geom_point() +
                        scale_color_continuous(low="red", high="blue",
                                               guide=guide_colorbar(reverse=TRUE)) +
                        scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                        scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      }else p1 <- NULL
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment_venn <- renderPlot({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){
      dotplot_for_output(data = selected_grange_venn_list(),
                         plot_genelist = pair_enrich1_H_venn(), Gene_set = input$Gene_set_venn, 
                         Species = input$Species_venn)
    }
  })
  
  
  output$download_venn_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H_venn()
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_venn_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set_venn,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1_venn()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  pair_enrich_table_venn <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1_venn()), 
                            H_t2g = Hallmark_set_venn(), Gene_set = input$Gene_set_venn))
  })
  
  output$venn_enrichment_result <- DT::renderDataTable({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){
      as.data.frame(enrichment_1_1_venn())
    }
  })
  
  
  
  output$whichGeneSet_venn <- renderUI({
    if(!is.null(input$intersection_venn_fun)){
      group <- input$intersection_venn_fun
      set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1_venn()), Group == group)$id)
      selectInput('Pathway_list_venn', 'Pathway list', set_list)
    }
  })
  
  region_gene_associate_venn <- reactive({
    if(!is.null(input$Pathway_list_venn)){
      set_list <- input$Pathway_list_venn
      df <- enrichment_enricher_venn()
      data3 <- selected_grange_venn_list()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_venn_associations <- DT::renderDT({
    if(!is.null(input$Pathway_list_venn) && 
       !is.null(input$intersection_venn_fun) && input$Species_venn != "not selected"){
      region_gene_associate_venn()
    }
  })
  output$whichGenes_venn <- renderUI({
    if(!is.null(input$venn_whichGroup2) && input$Species_venn != "not selected"){ 
      selectInput('intersection_venn_fun', 'Select intersection', names(enrichment_enricher_venn()))
    }
  })
  
  region_gene_associate_plot_venn <- reactive({
    set_list <- input$Pathway_list_venn
    df <- enrichment_enricher_venn()
    name <- input$intersection_venn_fun
    if(!is.null(df[[name]])) {
      res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
    }else res <- NULL
    return(res)
  })
  output$region_gene_associations_venn_plot <- renderPlot({
    if(!is.null(input$Pathway_list_venn) && 
       !is.null(input$intersection_venn_fun) && input$Species_venn != "not selected"){
      region_gene_associate_plot_venn()
    }
  })
  output$download_venn_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$intersection_venn_fun,"-",input$Pathway_list_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$venn_pdf_width
        set_list <- input$Pathway_list_venn
        df <- enrichment_enricher_venn()
        name <- input$intersection_venn_fun
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_venn_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations-",input$intersection_venn_fun,"-",input$Pathway_list_venn,".txt")
    },
    content = function(file){write.table(region_gene_associate_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  #Venn regulatory potential------
  output$venn_select_RNA <- renderUI({
    if(!is.null(Venn_peak_call_files())){
    selectInput("intersect_RNA", "Select intersect",
                c("not_selected",names(venn_overlap()$peaklist)),
                selected = "not_selected",multiple = F)
    }
  })

  peak_venn_grange_RNA  <- reactive({
    if(!is.null(input$intersect_RNA)){
    if(input$intersect_RNA != "not_selected"){
      return(venn_overlap()$peaklist[[input$intersect_RNA]])
    }}
  })
  
  
  output$vennRNAseqresult <- renderUI({
    fileInput("venn_DEG_result",
              "Select RNA-seq DEG result file",
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  RNAseqDEG_venn <- reactive({
    return(RNAseqDEGimport(tmp=input$venn_DEG_result$datapath,
                           exampleButton=input$goButton_venn))
  })
  
  output$venn_DEG_result <- renderDataTable({
    if(!is.null(RNAseqDEG_venn())){
      RNAseqDEG_venn() 
    }
  })
  
  RNAseqDEG_anno_venn <- reactive({
    return(RNAseqDEG_ann(RNAdata=RNAseqDEG_venn(),Species=input$Species_venn))
  })
  
  mmAnno_venn <- reactive({
    return(mmAnno(peak=peak_venn_grange_RNA(),
                  txdb=txdb_venn(),peak_distance=input$peak_distance_venn))
  })
  
  RP_venn <- reactive({
    return(RP_f(mmAnno=mmAnno_venn(),txdb=txdb_venn()))
  })
  regulatory_potential_venn <- reactive({
    return(regulatory_potential_f(species=input$Species_venn,data=RNAseqDEG_anno_venn(),
                                  result_geneRP= RP_venn(),DEG_fc=input$DEG_fc_venn,
                                  DEG_fdr=input$DEG_fdr_venn))
  })
  
  output$ks_plot_venn <- renderPlot({    
    if(!is.null(input$peak_distance_venn) && !is.null(RNAseqDEG_venn()) && 
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && 
       input$Species_venn != "not selected"  && !is.null(mmAnno_venn())){
      regulatory_potential_venn()
    }
  })
  output$DEG_fc_venn <- renderUI({
    numericInput("DEG_fc_venn","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr_venn <- renderUI({
    numericInput("DEG_fdr_venn","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  output$peak_distance_venn <- renderUI({
    sliderInput("peak_distance_venn","Regulatory range (distance (kb) from TSS)",
                min=0,max=200,value=100,step = 5)
  })
  
  output$RNAseqGroup_venn <- renderUI({
    if(input$Species_venn != "not selected" && !is.null(mmAnno_venn()) && !is.null(input$peak_distance_venn)){
      if(!is.null(RP_all_table_venn())){
        selectInput("RNAseqGroup_venn","Group (RNAseq-Epigenome)",
                    unique(RP_all_table_venn()$Group),
                    multiple = FALSE)
      }
    }
  })
  
  RNAseq_boxplot_venn <- reactive({
    RNA <- RNAseqDEG_anno_venn()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    data <- merge(RNA,RP_venn(), by="gene_id",all=T)
    data$group <- "Others"
    data$group[data$sumRP > 1] <- "RP > 1"
    data$group <- factor(data$group,levels=c("Others","RP > 1"),ordered=TRUE)
    
    collist <- unique(data$group)
    group1 <- dplyr::filter(data, group == collist[1])
    group2 <- dplyr::filter(data, group == collist[2])
    if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
      stat.test <- data %>% t_test(log2FoldChange ~ group)
      stat.test <- stat.test %>% add_significance()
      stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      col <-c("gray","#F8766D")
    }else stat.test <- NULL
    if(!is.null(stat.test)){
    p <- try(ggpubr::ggboxplot(data, x = "group", y = "log2FoldChange",
                           fill = "group", scales = "free", add = "jitter",
                           xlab = FALSE, ylab = "log2FoldChange")+theme_bw(base_size = 15)+
      xlab(NULL)+scale_fill_manual(values = col) + stat_pvalue_manual(stat.test,hide.ns = T, size = 5))
    }
    return(p)
  })
  output$RNAseq_boxplot_venn_error <- renderText({
    if(!is.null(input$peak_distance_venn) && !is.null(RNAseqDEG_venn()) && 
       input$Species_venn != "not selected" && !is.null(mmAnno_venn())){
      RNA <- RNAseqDEG_anno_venn()
      RNA <- dplyr::filter(RNA, !is.na(gene_id))
      data <- merge(RNA,RP_venn(), by="gene_id",all=T)
      data$group <- "Others"
      data$group[data$sumRP > 1] <- "RP > 1"
      data$group <- factor(data$group,levels=c("Others","RP > 1"),ordered=TRUE)
      
      collist <- unique(data$group)
      group1 <- dplyr::filter(data, group == collist[1])
      group2 <- dplyr::filter(data, group == collist[2])
      if(length(rownames(group1)) <= 1 || length(rownames(group2)) <= 1){
        print("boxplot: There are few genes with RP > 1")
    }
    }
  })
  
  RNAseq_popu_venn <- reactive({
    RNA <- RNAseqDEG_anno_venn()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$group[RNA$log2FoldChange > log2(input$DEG_fc_venn) & RNA$padj < input$DEG_fdr_venn] <- "up"
    RNA$group[RNA$log2FoldChange < -log2(input$DEG_fc_venn) & RNA$padj < input$DEG_fdr_venn] <- "down"
    data <- merge(RNA,RP_venn(), by="gene_id",all=T)
    if(length(which(unique(data$group) == "up")) == 1){
      up <- data %>% dplyr::filter(group == "up")
      up_total <- length(rownames(up))
      if(length(which(up$sumRP > 0))  > 0){
        up_epiUp <- length(rownames(dplyr::filter(up, sumRP > 0)))
      }else up_epiUp <- 0
      up_per <- c(up_total-(up_epiUp), up_epiUp)
    }else up_per <- c(0,0)
    if(length(which(unique(data$group) == "down")) == 1){
      down <- data %>% dplyr::filter(group == "down")
      down_total <- length(rownames(down))
      if(length(which(down$sumRP > 0))  > 0){
        down_epiUp <- length(rownames(dplyr::filter(down, sumRP > 0)))
      }else down_epiUp <- 0
      down_per <- c(down_total-(down_epiUp), down_epiUp)
    }else down_per <- c(0,0)
    x <- data.frame(RNA = c("up", "up","down","down"), 
                    Regulatory_potential = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:up (",up_per[2],")"),
                                                    paste0("down:NS (",down_per[1],")"),paste0("down:up (",down_per[2],")")),
                                                  levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:up (",up_per[2],")"),
                                                             paste0("down:NS (",down_per[1],")"),paste0("down:up (",down_per[2],")"))),
                    num = c(up_per,down_per),
                    col = c("F8766D","grey","F8766D","grey"))
    p <- ggplot(x, aes(x = RNA, y = num, fill = Regulatory_potential)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#F8766D","grey","#F8766D"))
    return(p)
  })
  output$int_box_venn <- renderPlot({
    withProgress(message = "Boxplot",{
      if(!is.null(input$peak_distance_venn) && !is.null(RNAseqDEG_venn()) && 
         input$Species_venn != "not selected" && !is.null(mmAnno_venn())){
        RNAseq_boxplot_venn()
      }
    })
  })
  output$int_bar_venn <- renderPlot({
    if(!is.null(input$peak_distance_venn) && !is.null(RNAseqDEG_venn()) && 
       !is.na(input$DEG_fdr_venn) && !is.na(input$DEG_fc_venn) && 
       input$Species_venn != "not selected" && !is.null(mmAnno_venn())){
      gridExtra::grid.arrange(RNAseq_popu_venn(), ChIPseq_popu_venn(), ncol = 1)
    }
  })
  ChIPseq_popu_venn <- reactive({
    RNA <- RNAseqDEG_anno_venn()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$gene_category <- "NS"
    RNA$gene_category[RNA$log2FoldChange > log2(input$DEG_fc_venn) & RNA$padj < input$DEG_fdr_venn] <- "up"
    RNA$gene_category[RNA$log2FoldChange < -log2(input$DEG_fc_venn) & RNA$padj < input$DEG_fdr_venn] <- "down"
    data <- merge(RNA,RP_venn(), by="gene_id")
    data$epigenome_category <- "NS"
    data$epigenome_category[data$sumRP > 0] <- "up"
    if(length(which(unique(data$epigenome_category) == "up")) == 1){
      up <- data %>% dplyr::filter(epigenome_category == "up")
      up_total <- length(rownames(up))
      if(length(which(unique(up$gene_category) == "up")) == 1){
        up_RNAUp <- length(rownames(dplyr::filter(up, gene_category == "up")))
      }else up_RNAUp <- 0
      if(length(which(unique(up$gene_category) == "down")) == 1){
        up_RNADown <- length(rownames(dplyr::filter(up, gene_category == "down")))
      }else up_RNADown <- 0
      up_per <- c(up_total-(up_RNAUp + up_RNADown), up_RNADown, up_RNAUp)
    }else{
      up_per <- c(0,0,0)
    } 
    x <- data.frame(Regulatory_potential = c("up", "up","up"), 
                    RNA = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")")),
                                 levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"))),
                    num = up_per,
                    col = c("#F8766D","#00BFC4","grey"))
    p <- ggplot(x, aes(x = Regulatory_potential, y = num, fill = RNA)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#00BFC4","#F8766D"))
    return(p)
  })
  RP_all_table_venn <- reactive({
    target_result <- regulatory_potential_venn()$data
    target_result$epigenome_category <- "up"
    table <- NULL
    if(!is.null(mmAnno_venn())) {
      table <- data.frame(Symbol = target_result$Symbol,
                          Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                          RNA_log2FC = -target_result$log2FoldChange,
                          RNA_padj = target_result$padj,
                          regulatory_potential = target_result$sumRP,
                          withUpPeakN = target_result$withPeakN,
                          gene_id = target_result$gene_id)
    }
    return(table)
  })
  
  RP_selected_table_venn <- reactive({
    table <- RP_all_table_venn() %>% dplyr::filter(Group == input$RNAseqGroup_venn)
    return(table)
  })
  
  output$RP_table_venn <- renderDT({
    if(!is.null(input$RNAseqGroup_venn) && 
       !is.null(input$peak_distance_venn && !is.null(mmAnno_venn()))){
      RP_selected_table_venn() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  
  output$download_vennintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu_venn(), ChIPseq_popu_venn(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_vennintbox = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential boxplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(RNAseq_boxplot_venn())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_vennKSplot = downloadHandler(
    filename = function(){
      paste0("KSplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(regulatory_potential_venn())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_RP_venn_table = downloadHandler(
    filename = function() {
      paste0("RP_summary_table.txt")
    },
    content = function(file){write.table(RP_all_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_selected_RP_venn_table = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup_venn,"_table.txt")
    },
    content = function(file){write.table(RP_selected_tabl_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  int_selected_bed_venn <- reactive({
    if(!is.null(RP_selected_table_venn())){
      gene <- RP_selected_table_venn()$gene_id
      y <- NULL
      if(!is.null(mmAnno_venn())) {
        up_peak <- subset(mmAnno_venn(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  output$download_selected_int_peak_venn = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup_venn,"_table.bed")
    },
    content = function(file){write.table(int_selected_bed_venn(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  observeEvent(input$venn_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_venn_RP", open="KS_panel")
  }))
  #Venn Integrative trackplot-------
  observeEvent(input$RP_table_venn_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_venn_RP", open="int_Trackplot_panel")
  }))
  output$int_igv_uprange_venn <- renderUI({
    if(!is.null(int_goi_gene_position_venn()) && !is.null(int_goi_promoter_position_venn())){
      y <- int_goi_promoter_position_venn()
      gene_position <- int_goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange_venn","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim_venn <- renderUI({
    numericInput("int_igv_ylim_venn","peak range:", value = 2, min = 0)
  })
  
  int_goi_promoter_position_venn<- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table_venn()[input$RP_table_venn_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(mmAnno_venn())) {
        up_peak <- subset(mmAnno_venn(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  
  int_goi_gene_position_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      gene <- RP_selected_table_venn()[input$RP_table_venn_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position_venn(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  output$int_trackplot_additional_venn <- renderUI({
    fileInput("int_trackplot_additional1_venn",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  int_track_additional_files_venn <-reactive({
    if(!is.null(input$int_trackplot_additional1_venn)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$int_trackplot_additional1_venn[, 1])){
        file <- input$int_trackplot_additional1_venn[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$int_trackplot_additional1_venn[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  int_data_track_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected)){
      return(data_trac(y=int_goi_promoter_position_venn(),gene_position=int_goi_gene_position_venn(),
                       gen=ref_venn(),txdb=txdb_venn(),org=org_venn(),
                       bw_files=bws_venn(),
                       track_additional_files=int_track_additional_files_venn()))
    }
  })
  
  int_highlight_trackplot_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected) && !is.null(input$int_igv_uprange_venn)){
      library(Gviz)
      y <- int_goi_promoter_position_venn()
      gene_position <- int_goi_gene_position_venn()
      chr <- gene_position$seqnames
      df <- int_data_track_venn()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange_venn[2] && y[i,]$end > input$int_igv_uprange_venn[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          if(y[i,]$Group == "up") {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == "down") {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  int_goi_trackplot_venn <- reactive({
    if(!is.null(input$RP_table_venn_rows_selected) &&
       !is.null(int_goi_promoter_position_venn()) && 
       !is.null(int_goi_gene_position_venn()) && 
       !is.null(gtrack_venn()) &&
       !is.null(input$int_igv_uprange_venn)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot_venn())){
        plot<- plotTracks(list(gtrack_venn(), int_highlight_trackplot_venn()),
                          from = input$int_igv_uprange_venn[1], 
                          to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                          type="hist")
      }else{
        df <- int_data_track_venn()
        df[["gtrack"]] <- gtrack_venn()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange_venn[1], 
                          to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi_venn <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_venn_rows_selected)){
        int_goi_trackplot_venn()
      }
    })
  })
  output$download_venn_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot_venn())){
          plot<- plotTracks(list(gtrack_venn(), int_highlight_trackplot_venn()),
                            from = input$int_igv_uprange_venn[1], 
                            to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                            type="hist")
        }else{
          df <- int_data_track_venn()
          df[["gtrack"]] <- gtrack_venn()
          plot<- plotTracks(df,
                            from = input$int_igv_uprange_venn[1], 
                            to = input$int_igv_uprange_venn[2],ylim=c(0,input$int_igv_ylim_venn),
                            type="hist")
        }
        ptirn(plot)
        dev.off()
        incProgress(1)
      })
    }
  )
  #Venn integrative functional enrichment-------
  int_Hallmark_set_venn <- reactive({
    if(!is.null(input$intGeneset_venn)){
      return(GeneList_for_enrichment(Species = input$Species_venn, Gene_set = input$intGeneset_venn, org = org_venn()))
    }
  })
  
  output$intGroup_venn <- renderUI({
    if(!is.null(RP_all_table_venn())){
      selectInput("intGroup_venn","Group (RNAseq-Epigenome)",unique(RP_all_table_venn()$Group),multiple = T)
    }
  })
  output$intGeneset_venn <- renderUI({
    selectInput('intGeneset_venn', 'Gene Set', gene_set_list)
  })
  
  selected_int_group_venn <- reactive({
    group <- input$intGroup_venn
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      table <- RP_all_table_venn() %>% dplyr::filter(Group == name)
      head(table)
      df2 <- data.frame(ENTREZID = table$gene_id, Group = table$Group)
      df <- rbind(df,df2)
    }
    return(df)
  })
  int_enrich_venn <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group_venn(), Species = input$Species_venn, org = org_venn(),
                                   H_t2g = int_Hallmark_set_venn(),Gene_set = input$intGeneset_venn))
  })
  int_enrich_list_venn <- reactive({
    return(enrich_gene_list(data = selected_int_group_venn(),
                            Gene_set = input$intGeneset_venn, org = org_venn(), H_t2g = int_Hallmark_set_venn()))
  })
  int_enrich_plot_venn <- reactive({
    return(enrich_genelist(data = selected_int_group_venn(),
                           enrich_gene_list = int_enrich_list_venn()))
  })
  output$int_enrichment1_venn <- renderPlot({
    dotplot_for_output(data = int_enrich_venn(),
                       plot_genelist = int_enrich_plot_venn(), Gene_set = input$intGeneset_venn, 
                       Species = input$Species_venn)
  })
  int_enrich_table_venn <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich_venn()), H_t2g = int_Hallmark_set_venn(), Gene_set = input$intGeneset_venn))
  })
  output$int_enrichment_result_venn <- DT::renderDataTable({
    int_enrich_table_venn()
  })
  output$download_venn_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset_venn,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich_venn(),
                           plot_genelist = int_enrich_plot_venn(), Gene_set = input$intGeneset_venn, 
                           Species = input$Species_venn)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_venn_int_enrichment_venn_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset_venn,".txt")
    },
    content = function(file){write.table(int_enrich_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ##Clustering--------
  # input data ------------------------------------------------------------------------------
  org1_clustering <- reactive({
    return(org(Species = input$Species_clustering))
  })
  txdb_clustering <- reactive({
    if(input$Species_clustering != "not selected"){
      return(txdb_function(Species = input$Species_clustering))
    }
  })
  promoter_region_clustering <- reactive({
    if(input$Genomic_region_clustering == "Promoter"){
      if(input$Species_clustering != "not selected"){
        return(promoter_clustering(txdb_clustering(),upstream = input$upstream_clustering, downstream = input$downstream_clustering))
      }
    }else return(promoter_clustering(upstream = input$upstream_clustering, downstream = input$downstream_clustering,
                          input_type = "Genome-wide",files =peak_call_files_clustering()))
  })
  gene_position_clustering <- reactive({
    if(input$Species_clustering != "not selected"){
      return(genes(txdb_clustering()))
    }
  })
  
  bws_clustering <- reactive({
    if(input$data_file_type_clustering == "Row1"){
      if(is.null(input$file1_clustering)){
        if(input$goButton_clustering > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file1_clustering[, 1])){
          file <- input$file1_clustering[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$file1_clustering[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
        return(files)
      }
    }
  })
  bam_clustering <- reactive({
    if(input$data_file_type_clustering == "Row2"){
      if(is.null(input$file_bam_clustering)){
        if(input$goButton_clustering > 0 ){
          df<-list()
          df[["A_1"]] <- "data/bigwig/A_1.BigWig"
          df[["A_2"]] <- "data/bigwig/A_2.BigWig"
          df[["B_1"]] <- "data/bigwig/B_1.BigWig"
          df[["B_2"]] <- "data/bigwig/B_2.BigWig"
          return(df)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$file_bam_clustering[, 1])){
          file <- input$file_bam_clustering[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$file_bam_clustering[nr,]$name))
          files <- c(files,file)
        }
        names(files)<-name
        return(files)
      }
    }
  })
  
  peak_call_files_clustering <- reactive({
    if(input$Genomic_region_clustering == "Genome-wide"){
      if(is.null(input$peak_call_file1_clustering)){
        if(input$goButton_clustering > 0 ){
          df <- list()
          df[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
          df[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
          df[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
          df[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
          name <- c("Ctrl_1","Ctrl_2","Sen_1","Sen_2")
          files2 <- lapply(df, GetGRanges, simple = TRUE)
          names(files2)<-name
          return(files2)
        }
        return(NULL)
      }else{
        files<-c()
        name<-c()
        for(nr in 1:length(input$peak_call_file1_clustering[, 1])){
          file <- input$peak_call_file1_clustering[[nr, 'datapath']]
          name <- c(name, gsub("\\..+$", "", input$peak_call_file1_clustering[nr,]$name))
          files <- c(files,file)
        }
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        names(files2)<-name
        return(files2)
      }
    }
  })
  
  bw_files_clustering <- reactive({
    return(BigWigFileList(bws_clustering()))
  })
  
  output$input_bw_files_clustering <- DT::renderDT({
    if(input$data_file_type_clustering == "Row1"){
      uploaded_files = names(bws_clustering())
      data.frame(uploaded_files = uploaded_files)
    }else{
      uploaded_files = names(bam_clustering())
      as.data.frame(uploaded_files) 
    }
  })
  
  output$input_peak_call_files_clustering <- DT::renderDataTable({
    if(input$Genomic_region_clustering == "Genome-wide"){
      uploaded_files = names(peak_call_files_clustering())
      as.data.frame(uploaded_files)
    }
  })
  genelist_clustering <- reactive({
    withProgress(message = "Importing a gene list file, please wait",{
      tmp <- input$genelist_file1_clustering$datapath
      if(is.null(input$genelist_file1_clustering) && input$goButton_clustering > 0 )  tmp = "data/RNAseq.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
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
  })
  bw_count_clustering2 <- reactive({
    if(input$data_file_type_clustering == "Row1" && !is.null(bw_files_clustering())){
      if(input$Genomic_region_clustering == "Promoter"){
        if(input$Species_clustering != "not selected"){
          count <-Bigwig2count(bw = bw_files_clustering(),promoter_region_clustering(),
                               Species = input$Species_clustering,input_type =input$Genomic_region_clustering)
          return(count) 
        }
      }else{
        if(!is.null(peak_call_files_clustering())){
          return(Bigwig2count(bw = bw_files_clustering(),promoter_region_clustering(),
                              input_type =input$Genomic_region_clustering))
        }
      }
    }
    if(input$Species_clustering != "not selected" && input$data_file_type_clustering == "Row2" && 
       !is.null(bam_clustering())){
      withProgress(message = "converting bam to gene count data",{
        switch(input$Pair_or_single_clustering,
               "Paired-end" = seq <- TRUE,
               "Single-end" = seq <- FALSE)
        regionsToCount <- data.frame(GeneID = paste("ID", seqnames(promoter_region_clustering()), 
                                                    start(promoter_region_clustering()), end(promoter_region_clustering()), sep = "_"), 
                                     Chr = seqnames(promoter_region_clustering()), 
                                     Start = start(promoter_region_clustering()), End = end(promoter_region_clustering()), 
                                     Strand = strand(promoter_region_clustering()))
        a <- as.data.frame(promoter_region_clustering())
        Row.name <- paste0(a$seqnames,":",a$start,"-",a$end)
        if(input$Genomic_region_clustering == "Promoter"){
          count <- featureCounts(bam_clustering(),annot.ext = regionsToCount,isPairedEnd = seq,
                                 countMultiMappingReads =F,maxFragLength = 100)$counts
          rownames(count) <- Row.name
          colnames(count) <- names(bam_clustering())
          return(count)
        }else{
          if(!is.null(peak_call_files_clustering())){
            count <- featureCounts(bam_clustering(),annot.ext = regionsToCount,isPairedEnd = seq,
                                   countMultiMappingReads =F,maxFragLength = 100)$counts
            rownames(count) <- Row.name
            colnames(count) <- names(bam_clustering())
            return(count)
          }
        }
      })
    }
  })
  bw_count_clustering <- reactive({
    if(input$data_file_type_clustering == "Row1" && !is.null(bw_files_clustering())){
      if(input$Genomic_region_clustering == "Promoter"){
        if(input$Species_clustering != "not selected"){
          count <-bw_count_clustering2()
          if(!is.null(genelist_clustering())){
            colnum <- length(colnames(count))
            count <- merge(count,genelist_clustering(),by=0)
            rownames(count) <- count$Row.names
            count <- count[,2:(1+colnum)]
          }
          return(count) 
        }
      }else{
          return(bw_count_clustering2())
      }
    }
    if(input$Species_clustering != "not selected" && input$data_file_type_clustering == "Row2" && 
       !is.null(bam_clustering())){
      return(bw_count_clustering2())
    }
  })
  output$raw_count_table_clustering <- DT::renderDataTable({
    bw_count_clustering()
  })
  
  output$download_raw_count_clustering_table <- downloadHandler(
    filename = function() {
      if (input$Genomic_region_clustering=='Promoter'){
        paste0("bigwig2count-promoter(", -input$upstream,"-",input$downstream,").txt")
      }else{
        paste0("bigwig2count-genomeWide",".txt")
      }},
    content = function(file){write.table(bw_count_clustering(), file, row.names = T, sep = "\t", quote = F)}
  )
  observeEvent(bws_clustering(), ({
    updateCollapse(session,id =  "input_collapse_panel_clustering", open="raw_count_panel")
  }))
  observeEvent(bam_clustering(), ({
    updateCollapse(session,id =  "input_collapse_panel_clustering", open="raw_count_panel")
  }))
  observeEvent(peak_call_files_clustering(), ({
    updateCollapse(session,id =  "input_collapse_panel_clustering", open="peak_call_files_panel")
  }))
  
  # clustering PCA ------------------------------------------------------------------------------
  output$clustering_pca_error <- renderText({
    p <- length(colnames(bw_count_clustering()))
    if(p < 3){
        print("PCA: The sample number must be > 2")
    }else return(NULL)
  })
   output$download_clustering_PCA = downloadHandler(
    filename = "PCA-MDS-dendrogram.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(PCAplot(data = bw_count_clustering(),plot=TRUE))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA_clustering <- renderPlot({
    withProgress(message = "preparing result table",{
      if(is.null(bw_count_clustering())){
        return(NULL)
      }else{
        p <- length(colnames(bw_count_clustering()))
        if(p > 2){
          print(PCAplot(data = bw_count_clustering(),plot=TRUE))
        }
      }
    })
  })
  
  output$clustering_PCA_data <- DT::renderDataTable({
    PCAplot(data = bw_count_clustering(),plot=FALSE)
  })
  
  output$download_clustering_PCA_table = downloadHandler(
    filename ="PCA_table.txt",
    content = function(file){write.table(PCAplot(data = bw_count_clustering(),plot=FALSE), 
                                          file, row.names = T, sep = "\t", quote = F)}
  )
  #Clustering umap-------
  output$clustering_umap_n <- renderUI({
    sliderInput("clustering_n_neighbors", "n_neighbors", min = 2,
                max=100, step = 1,
                value = 15)
  })
  
  clustering_umap_plot <- reactive({
    data <- bw_count_clustering()
    if(is.null(input$clustering_n_neighbors)){
      return(NULL)
    }else{
      if(is.null(data)){
        return(NULL)
      }else{
        p<- try(umap_plot(data = data, n_neighbors = input$clustering_n_neighbors))
        return(p)
      }
    }
  })
  
  output$clustering_umap_error <- renderText({
    p <- clustering_umap_plot()
    if(length(p) == 1){
      if (class(p) == "try-error") {
        print("umap: number of neighbors must be smaller than number of items")
      }
    }else return(NULL)
  })
  
  
  output$clustering_umap <- renderPlot({
    p <- clustering_umap_plot()
    withProgress(message = "umap",{
      if(length(p) == 1){
        if (class(p) == "try-error") {
          return(NULL)
        }
      }else{
        print(p)
      }
    })
  })
  
  output$download_clustering_umap = downloadHandler(
    filename = function(){
      paste0(input$clustering_n_neighbors,"_neighbors_umap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 4.7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(clustering_umap_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  #Clustering correlation plot-----
  corrplot <-reactive({
    data <- bw_count_clustering()
    mcor <- cor(data)
    p <- ggcorrplot(corr = mcor, hc.order = TRUE, method = "square",
                    colors = c("#4b61ba", "white", "red"), lab = TRUE)
    p <- p + scale_fill_gradient2(limit = c(0,1), low = "#4b61ba", high =  "red", mid = "white", midpoint = 0.5)
    return(p)
  })
  output$correlationplot <- renderPlot({
    if(!is.null(bw_count_clustering())){
      corrplot()
    }
  })
  
  output$download_clustering_corrplot = downloadHandler(
    filename = function(){
      paste0("Corrplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(corrplot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #Clustering kmeans--------
  bw_count_clustering_anno <- reactive({
    if(!is.null(bw_count_clustering()) && !is.null(txdb_clustering())){
      withProgress(message = "preparing annotation",{
        data <- peak_kmeans_grange()
        data2 <- annotatePeak(data, TxDb= txdb_clustering())
        data <- as.data.frame(as.GRanges(data2))
        Row.name <- paste0(data$seqnames,":",data$start,"-",data$end)
        data$locus <- Row.name
        my.symbols <- data$geneId
        gene_IDs<-id_convert(my.symbols,input$Species_clustering,type="ENTREZID")
        colnames(gene_IDs) <- c("geneId","NearestGene")
        data <- merge(data, gene_IDs, by="geneId")
        data <- data %>% distinct(locus, .keep_all = T)
        rownames(data)<-data$locus
        return(data)
      })
    }
  })
  
  bw_count_clusterin_anno2 <- reactive({
    if(input$Genomic_region_clustering == "Genome-wide"){
      anno <- bw_count_clustering_anno()
      anno <- dplyr::arrange(anno, locus)
      data <- data.frame(NearestGene = anno$NearestGene,
                         locus = anno$locus)
      clustering_kmeans_pattern_extract <- clustering_kmeans_pattern_extract()
      clustering_kmeans_pattern_extract$locus <- rownames(clustering_kmeans_pattern_extract())
      clustering_kmeans_pattern_extract <- dplyr::arrange(clustering_kmeans_pattern_extract, locus)
      data2 <- merge(data,clustering_kmeans_pattern_extract,by="locus")
      rownames(data2) <- data2$locus
      data2 <- data2[, - which(colnames(data2) == "locus")]
      return(data2)
    }
  })
  
  output$clustering_kmeans_num <- renderUI({
    if(is.null(bw_count_clustering())){
      return(NULL)
    }else{
      withProgress(message = "Preparing kmeans clustering",{
        sliderInput("clustering_kmeans_number", "k-means number", min = 1,
                    max=20, step = 1,
                    value = 2)
      })
    }
  })
  
  output$kmeans_cv <- renderUI({
    sliderInput("kmeans_cv", "Most variable loci:", min = 0,
                max=10000, step = 1000,
                value = 2000)
  })
  
  
  bw_count_clustering_cutoff <- reactive({
    data <- bw_count_clustering()
    if(is.null(data) || is.null(input$kmeans_cv)){
      return(NULL)
    }else{
      data2 <- data[order(apply(data,1,mad), decreasing = T)[1:input$kmeans_cv],]
      return(data2)
    }
  })
  
  clustering_data_z <- reactive({
    data <- bw_count_clustering_cutoff()
    if(is.null(data)){
      return(NULL)
    }else{
      data.z <- genescale(data, axis = 1, method = "Z")
      data.z <- na.omit(data.z)
      return(data.z)
    }
  })
  
  clustering_kmeans <- reactive({
    data.z <- clustering_data_z()
    if(is.null(data.z)){
      return(NULL)
    }else{
      withProgress(message = "k-means clustering",{
        ht <- Heatmap(data.z, name = "z-score",
                      column_order = colnames(data.z),
                      clustering_method_columns = 'ward.D2',
                      row_km= input$clustering_kmeans_number, cluster_row_slices = F, row_km_repeats = 100,
                      show_row_names = F,column_names_side = "top")
        ht <- draw(ht)
        return(ht)
      })
    }
  })
  
  clustering_kmeans_cluster <- reactive({
    ht <- clustering_kmeans()
    data.z <- clustering_data_z()
    data <- bw_count_clustering_cutoff()
    if(is.null(ht) || is.null(data.z)){
      return(NULL)
    }else{
      r.dend <- row_dend(ht)
      rcl.list <- row_order(ht)
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      if(!is.null(input$clustering_kmeans_number)){
        if(length(lapply(rcl.list, function(x) length(x))) != input$clustering_kmeans_number){
          return(NULL)
        }else{
          for (i in 1:length(row_order(ht))){ if (i == 1) {
            clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
            out <- cbind(clu, paste("cluster", i, sep=""))
            colnames(out) <- c("GeneID", "Cluster")} else {
              clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
              clu <- cbind(clu, paste("cluster", i, sep=""))
              out <- rbind(out, clu)}}
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  output$clustering_select_kmean <- renderUI({
    withProgress(message = "preparing kmeans clustering",{
    clusters <- clustering_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("clustering_select_kmean", "cluster_list", choices = c(unique(clusters$Cluster)),multiple = T)
    } 
    })
  })
  
  clustering_kmeans_pattern_extract <- reactive({
    clusters <- clustering_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$clustering_select_kmean)){
        return(NULL)
      }else{
        cluster_name <- input$clustering_select_kmean
        clusterCount <- dplyr::filter(clusters, Cluster %in% cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  output$clustering_kmeans_extract_table <- renderDT({
    if(!is.null(input$clustering_select_kmean)){
      if(input$Genomic_region_clustering == "Genome-wide"){
    if(input$Species_clustering == "not selected"){
      clustering_kmeans_pattern_extract() %>%
        datatable(
          selection = "single")
    }else{
      if(!is.null(peak_kmeans_grange()) ){
      bw_count_clusterin_anno2() %>%
          datatable(
            selection = "single")
      }
    } 
      }else{
        clustering_kmeans_pattern_extract() %>%
          datatable(
            selection = "single")
      }
    }
  })
  
  output$clustering_kmeans_heatmap <- renderPlot({
    withProgress(message = "plot heatmap",{
    ht <- clustering_kmeans()
    if(is.null(ht)){
      return(NULL)
    }else{
      print(ht)
    }
    })
  })
  
  output$download_clustering_kmeans_extract_count_bed = downloadHandler(
    filename = function() {
      paste0("kmeans_selected_table",".bed")
    },
    content = function(file){
        write.table(as.data.frame(peak_kmeans_grange()), file, row.names = F, col.names = F,sep = "\t", quote = F)
      }
  )
  output$download_clustering_kmeans_heatmap = downloadHandler(
    filename = function() {
      paste(input$clustering_kmeans_number,"kmeans_heatmap.pdf",sep = "_")
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(clustering_kmeans())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_clustering_kmeans_extract_count = downloadHandler(
    filename = function() {
      paste0("kmeans_selected_table",".txt")
    },
    content = function(file){
      if(input$Species_clustering == "not selected"){
      write.table(clustering_kmeans_pattern_extract(), file, row.names = T, sep = "\t", quote = F)
      }else write.table(bw_count_clusterin_anno2(), file, row.names = T, sep = "\t", quote = F)
        }
  )
  #clustering kmeans peak--------
  output$peak_pattern_kmeans_additional <- renderUI({
    fileInput("peak_pattern_kmeans_add",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  kmeans_additional <-reactive({
    if(!is.null(input$peak_pattern_kmeans_add)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$peak_pattern_kmeans_add[, 1])){
        file <- input$peak_pattern_kmeans_add[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$peak_pattern_kmeans_add[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  output$peak_pattern_kmeans_heat_range <- renderUI({
    numericInput("peak_kmeans_range","Intensity range",value=10,min = 0.1,step=2)
  })
  peak_kmeans_grange <- reactive({
    if(!is.null(clustering_kmeans_pattern_extract())){
      if(input$Genomic_region_clustering == "Genome-wide"){
        kmeans <- range_changer(clustering_kmeans_pattern_extract())
        kmeans <- with(kmeans, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
      }else{
        kmeans <- symbol2gene_id(clustering_kmeans_pattern_extract(),org1_clustering()) %>% distinct(gene_id, .keep_all = T)
        kmeans <- subset(promoter_region_clustering(), gene_id %in% kmeans$gene_id) 
      }
      return(kmeans)
    }else return(NULL)
  })
  peak_kmeans_cvglists <-reactive({
    kmeans2 <- peak_kmeans_grange()
    feature.recentered <- reCenterPeaks(kmeans2, width=4000)
    files <- bws_clustering()
    if(!is.null(kmeans_additional())) files <- c(files, kmeans_additional())
    cvglists <- sapply(files, import,which=feature.recentered,as="RleList")
    names(cvglists) <- gsub(".+\\/","",gsub("\\..+$", "", names(files)))
    return(cvglists)
  })
  
  peak_kmeans_sig <- reactive({
    kmeans2 <- peak_kmeans_grange()
    cvglists <- peak_kmeans_cvglists()
    feature.center <- reCenterPeaks(kmeans2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000)
    return(sig)
  })
  
  peak_kmeans_feature_center <- reactive({
    kmeans2 <- peak_kmeans_grange()
    feature.center <- reCenterPeaks(kmeans2, width=1)
    return(feature.center)
  })
  peak_kmeans_alinedHeatmap <- reactive({
    if(!is.null(input$peak_kmeans_range)){
      range <- c()
      for(n in length(names(bws_clustering()))) {range <- c(range, input$peak_kmeans_range)}
      heatmap <- featureAlignedHeatmap(peak_kmeans_sig(), peak_kmeans_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = brewer.pal(n = 9, name = 'OrRd'))
      return(heatmap)
    }
  })
  output$peak_pattern_kmeans_heatmap <- renderPlot({
    withProgress(message = "feature aligned heatmap",{
      if(!is.null(clustering_kmeans_pattern_extract())){
        plot_grid(peak_kmeans_alinedHeatmap())
      }
    })
  })
  output$peak_pattern_kmeans_line <- renderPlot({
    withProgress(message = "feature aligned distribution",{
      if(!is.null(clustering_kmeans_pattern_extract())){
        featureAlignedDistribution(peak_kmeans_sig(), peak_kmeans_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
      }
    })
  })
  
  output$download_peak_pattern_kmeans_heatmap = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_heatmap",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(peak_kmeans_alinedHeatmap()))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_peak_pattern_kmeans_line = downloadHandler(
    filename = function(){
      paste0("download_peak_pattern_alignedDistibution",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        featureAlignedDistribution(peak_kmeans_sig(), peak_kmeans_feature_center(),
                                   upstream=2000, downstream=2000,
                                   type="l")
        dev.off()
        incProgress(1)
      })
    }
  )
  #Clustering trackplot-------
  ref_clustering <- reactive({
    ref <- gsub(".+\\(","",gsub(")", "", input$Species_clustering))
    return(ref)
  })
  

  observeEvent(input$clustering_kmeans_extract_table_rows_selected, ({
    updateCollapse(session,id =  "clustering_kmeans_collapse_panel", open="kmeans_track_panel")
  }))
  observeEvent(input$clustering_kmeans_extract_table_rows_selected,({
    if(input$Species_clustering != "not selected"){
    if(!is.null(goi_gene_position_clustering()) && !is.null(goi_promoter_position_clustering())){
      y <- goi_promoter_position_clustering()
      gene_position <- goi_gene_position_clustering()
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
      updateSliderInput(session,"igv_uprange_clustering","Range:",
                        value = c(start_position,end_position),
                        step = 100, min = start_position - 10000, max = end_position + 10000)
    }
    }
  }))
  
  output$igv_uprange_clustering <- renderUI({
    if(input$Species_clustering != "not selected"){
    if(!is.null(goi_gene_position_clustering()) && !is.null(goi_promoter_position_clustering())){
      y <- goi_promoter_position_clustering()
      gene_position <- goi_gene_position_clustering()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("igv_uprange_clustering","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
    }
  })
  output$igv_ylim_clustering <- renderUI({
    numericInput("igv_ylim_clustering","peak range:", value = 2, min = 0)
  })
  
  gtrack_clustering <- reactive({
    withProgress(message = "Preparing track",{
      library(Gviz)
      gtrack <- Gviz::GenomeAxisTrack(cex=1)
      return(gtrack)
    })
  })
  
  goi_promoter_position_clustering<- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      library(Gviz)
      if(input$Genomic_region_clustering == "Promoter"){
        label_data <- rownames(clustering_kmeans_pattern_extract()[input$clustering_kmeans_extract_table_rows_selected,])
        gene_IDs<- id_convert(label_data,input$Species_clustering,type="SYMBOL_single")
        y <- as.data.frame(subset(promoter_region_clustering(), gene_id %in% gene_IDs))
      }else{
        data <- range_changer(bw_count_clusterin_anno2()[input$clustering_kmeans_extract_table_rows_selected,])
        y <- as.data.frame(with(data, GRanges(seqnames = chr, 
                                 ranges = IRanges(start,end))))
      }
      return(y)
    }
  })
  
  goi_gene_position_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      if(input$Genomic_region_clustering == "Promoter"){
        label_data <- rownames(clustering_kmeans_pattern_extract()[input$clustering_kmeans_extract_table_rows_selected,])
      }else{
        label_data <- bw_count_clusterin_anno2()[input$clustering_kmeans_extract_table_rows_selected,]$NearestGene
      }
      gene_IDs<- id_convert(label_data,input$Species_clustering,type="SYMBOL_single")
      gene_position <- as.data.frame(subset(gene_position_clustering(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  output$trackplot_additional_clustering <- renderUI({
    fileInput("trackplot_additional1_clustering",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  track_additional_files_clustering <-reactive({
    if(!is.null(input$trackplot_additional1_clustering)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$trackplot_additional1_clustering[, 1])){
        file <- input$trackplot_additional1_clustering[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$trackplot_additional1_clustering[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  data_track_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      return(data_trac(y=goi_promoter_position_clustering(),gene_position=goi_gene_position_clustering(),
                       gen=ref_clustering(),txdb=txdb_clustering(),org=org1_clustering(),filetype=input$data_file_type_clustering,
                       bw_files=bws_clustering(),bam_files=bam_clustering(),
                       track_additional_files=track_additional_files_clustering()))
    }
  })
  
  highlight_trackplot_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected)){
      library(Gviz)
      y <- goi_promoter_position_clustering()
      gene_position <- goi_gene_position_clustering()
      chr <- gene_position$seqnames
      df <- data_track_clustering()
      if(y$start < input$igv_uprange_clustering[2] && y$end > input$igv_uprange_clustering[1]){
        withProgress(message = "Highlighting the selected region",{
          ht <- HighlightTrack(trackList = df,
                               start = y$start, width = y$width,
                               chromosome = chr, alpha=0.5)
        })
      }else ht <- NULL
      return(ht)
    }
  })
  goi_trackplot_clustering <- reactive({
    if(!is.null(input$clustering_kmeans_extract_table_rows_selected) &&
       !is.null(goi_promoter_position_clustering()) && 
       !is.null(goi_gene_position_clustering()) && 
       !is.null(gtrack_clustering()) &&
       !is.null(input$igv_uprange_clustering)){
      library(Gviz)
      if(!is.null(highlight_trackplot_clustering())){
        plot<- plotTracks(list(gtrack_clustering(), highlight_trackplot_clustering()),
                          from = input$igv_uprange_clustering[1], 
                          to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                          type="hist",cex.main = 1.25)
      }else{
        df <- data_track_clustering()
        df[["gtrack"]] <- gtrack_clustering()
        plot<- plotTracks(df,
                          from = input$igv_uprange_clustering[1], 
                          to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                          type="hist",cex.main = 1.25)
      }
      return(plot)
    }
  })
  output$trackplot_goi_clustering <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$clustering_kmeans_extract_table_rows_selected)  &&
         input$Species_clustering != "not selected"){
        goi_trackplot_clustering()
      }
    })
  })
  
  output$download_clustering_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(highlight_trackplot_clustering())){
          plotTracks(list(gtrack_clustering(), highlight_trackplot_clustering()),
                     from = input$igv_uprange_clustering[1], 
                     to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                     type="hist")
        }else{
          df <- data_track_clustering()
          df[["gtrack"]] <- gtrack_clustering()
          plotTracks(df,
                     from = input$igv_uprange_clustering[1], 
                     to = input$igv_uprange_clustering[2],ylim=c(0,input$igv_ylim_clustering),
                     type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  #Clustering regulatory potential------
  output$warning_kmeans_RNA <- renderText({
    if(is.null(clustering_kmeans_cluster())) print("Please press the 'k-means clustering panel' in advance")
  })
  
  output$clustering_select_kmean_RNA <- renderUI({
    withProgress(message = "preparing kmeans clustering",{
      clusters <- clustering_kmeans_cluster()
      if(is.null(clusters)){
        return(NULL)
      }else{
        selectInput("clustering_select_kmean_RNA", "k-means cluster_list", choices = c(unique(clusters$Cluster)),multiple = T)
      } 
    })
  })
  
  clustering_kmeans_pattern_extract_RNA <- reactive({
    clusters <- clustering_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$clustering_select_kmean_RNA)){
        return(NULL)
      }else{
        cluster_name <- input$clustering_select_kmean_RNA
        clusterCount <- dplyr::filter(clusters, Cluster %in% cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  peak_kmeans_grange_RNA <- reactive({
    if(!is.null(clustering_kmeans_pattern_extract_RNA())){
      if(input$Genomic_region_clustering == "Genome-wide"){
        kmeans <- range_changer(clustering_kmeans_pattern_extract_RNA())
        kmeans <- with(kmeans, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
      }else{
        kmeans <- symbol2gene_id(clustering_kmeans_pattern_extract_RNA(),org1_clustering()) %>% distinct(gene_id, .keep_all = T)
        kmeans <- subset(promoter_region_clustering(), gene_id %in% kmeans$gene_id) 
      }
      return(kmeans)
    }else return(NULL)
  })
  
  
  output$clusteringRNAseqresult <- renderUI({
    fileInput("clustering_DEG_result",
              "Select RNA-seq DEG result file",
              accept = c("txt","csv","xlsx"),
              multiple = TRUE,
              width = "80%")
  })
  RNAseqDEG_clustering <- reactive({
    return(RNAseqDEGimport(tmp=input$clustering_DEG_result$datapath,
                           exampleButton=input$goButton_clustering))
  })
  
  output$clustering_DEG_result <- renderDataTable({
    if(!is.null(RNAseqDEG_clustering())){
      RNAseqDEG_clustering() 
    }
  })
  
  RNAseqDEG_anno_clustering <- reactive({
    return(RNAseqDEG_ann(RNAdata=RNAseqDEG_clustering(),Species=input$Species_clustering))
  })
  
  mmAnno_clustering <- reactive({
    return(mmAnno(peak=peak_kmeans_grange_RNA(),genomic_region=input$Genomic_region_clustering,
                  txdb=txdb_clustering(),peak_distance=input$peak_distance_clustering))
  })
  
  RP_clustering <- reactive({
    return(RP_f(mmAnno=mmAnno_clustering(),txdb=txdb_clustering()))
  })
  regulatory_potential_clustering <- reactive({
      return(regulatory_potential_f(species=input$Species_clustering,data=RNAseqDEG_anno_clustering(),
                                  result_geneRP= RP_clustering(),DEG_fc=input$DEG_fc_clustering,
                                  DEG_fdr=input$DEG_fdr_clustering))
  })
  
  output$ks_plot_clustering <- renderPlot({    
    if(!is.null(input$peak_distance_clustering) && !is.null(RNAseqDEG_clustering()) && 
       !is.na(input$DEG_fdr_clustering) && !is.na(input$DEG_fc_clustering) && !is.null(bw_count_clustering()) && 
       input$Species_clustering != "not selected"  && !is.null(mmAnno_clustering())){
      regulatory_potential_clustering()
    }
  })
  output$DEG_fc_clustering <- renderUI({
    numericInput("DEG_fc_clustering","Fold change cutoff for RNA-seq data",
                 min=0,max=NA,value=1.5,step = 0.5)
  })
  output$DEG_fdr_clustering <- renderUI({
    numericInput("DEG_fdr_clustering","FDR cutoff for RNA-seq data",
                 min=0,max=1, value=0.05,step = 0.001)
  })
  output$peak_distance_clustering <- renderUI({
    sliderInput("peak_distance_clustering","Regulatory range (distance (kb) from TSS)",
                min=0,max=200,value=100,step = 5)
  })
  
  output$RNAseqGroup_clustering <- renderUI({
    if(input$Species_clustering != "not selected" &&!is.null(mmAnno_clustering()) && !is.null(input$peak_distance_clustering)){
      if(!is.null(RP_all_table_clustering())){
        selectInput("RNAseqGroup_clustering","Group (RNAseq-Epigenome)",
                    unique(RP_all_table_clustering()$Group),
                    multiple = FALSE)
      }
    }
  })
  
  RNAseq_boxplot_clustering <- reactive({
    RNA <- RNAseqDEG_anno_clustering()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    data <- merge(RNA,RP_clustering(), by="gene_id",all=T)
    data$group <- "Others"
    data$group[data$sumRP > 1] <- "RP > 1"
    data$group <- factor(data$group,levels=c("Others","RP > 1"),ordered=TRUE)
    
    collist <- unique(data$group)
    group1 <- dplyr::filter(data, group == collist[1])
    group2 <- dplyr::filter(data, group == collist[2])
    if(length(rownames(group1)) >1 && length(rownames(group2)) >1){
      stat.test <- data %>% t_test(log2FoldChange ~ group)
      stat.test <- stat.test %>% add_significance()
      stat.test <- stat.test %>% add_xy_position(scales = "free", step.increase = 0.2)
      col <-c("gray","#F8766D")
    }else stat.test <- NULL
    if(!is.null(stat.test)){
    p <- try(ggpubr::ggboxplot(data, x = "group", y = "log2FoldChange",
                           fill = "group", scales = "free", add = "jitter",
                           xlab = FALSE, ylab = "log2FoldChange")+theme_bw(base_size = 15)+
      xlab(NULL)+scale_fill_manual(values = col) + stat_pvalue_manual(stat.test,hide.ns = T, size = 5))
    }
    return(p)
  })
  output$RNAseq_boxplot_clustering_error <- renderText({
    if(!is.null(input$peak_distance_clustering) && !is.null(RNAseqDEG_clustering()) && 
       !is.na(input$DEG_fdr_clustering) && !is.na(input$DEG_fc_clustering) && 
       !is.null(bw_count_clustering()) && input$Species_clustering != "not selected" && 
       !is.null(mmAnno_clustering())){
      RNA <- RNAseqDEG_anno_clustering()
      RNA <- dplyr::filter(RNA, !is.na(gene_id))
      data <- merge(RNA,RP_clustering(), by="gene_id",all=T)
      data$group <- "Others"
      data$group[data$sumRP > 1] <- "RP > 1"
      data$group <- factor(data$group,levels=c("Others","RP > 1"),ordered=TRUE)
      
      collist <- unique(data$group)
      group1 <- dplyr::filter(data, group == collist[1])
      group2 <- dplyr::filter(data, group == collist[2])
      if(length(rownames(group1)) <= 1 || length(rownames(group2)) <= 1){
        print("boxplot: There are few genes with RP > 1")
      }
    }
  })
  
  RNAseq_popu_clustering <- reactive({
    RNA <- RNAseqDEG_anno_clustering()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$group[RNA$log2FoldChange > log2(input$DEG_fc_clustering) & RNA$padj < input$DEG_fdr_clustering] <- "up"
    RNA$group[RNA$log2FoldChange < -log2(input$DEG_fc_clustering) & RNA$padj < input$DEG_fdr_clustering] <- "down"
    data <- merge(RNA,RP_clustering(), by="gene_id",all=T)
    if(length(which(unique(data$group) == "up")) == 1){
      up <- data %>% dplyr::filter(group == "up")
      up_total <- length(rownames(up))
      if(length(which(up$sumRP > 0))  > 0){
        up_epiUp <- length(rownames(dplyr::filter(up, sumRP > 0)))
      }else up_epiUp <- 0
      up_per <- c(up_total-(up_epiUp), up_epiUp)
    }else up_per <- c(0,0)
    if(length(which(unique(data$group) == "down")) == 1){
      down <- data %>% dplyr::filter(group == "down")
      down_total <- length(rownames(down))
      if(length(which(down$sumRP > 0))  > 0){
        down_epiUp <- length(rownames(dplyr::filter(down, sumRP > 0)))
      }else down_epiUp <- 0
      down_per <- c(down_total-(down_epiUp), down_epiUp)
    }else down_per <- c(0,0)
    x <- data.frame(RNA = c("up", "up","down","down"), 
                    Regulatory_potential = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:up (",up_per[2],")"),
                                                    paste0("down:NS (",down_per[1],")"),paste0("down:up (",down_per[2],")")),
                                                  levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:up (",up_per[2],")"),
                                                             paste0("down:NS (",down_per[1],")"),paste0("down:up (",down_per[2],")"))),
                    num = c(up_per,down_per),
                    col = c("F8766D","grey","F8766D","grey"))
    p <- ggplot(x, aes(x = RNA, y = num, fill = Regulatory_potential)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#F8766D","grey","#F8766D"))
    return(p)
  })
  output$int_box_clustering <- renderPlot({
    withProgress(message = "Boxplot",{
    if(!is.null(input$peak_distance_clustering) && !is.null(RNAseqDEG_clustering()) && 
       !is.na(input$DEG_fdr_clustering) && !is.na(input$DEG_fc_clustering) && 
       !is.null(bw_count_clustering()) && input$Species_clustering != "not selected" && 
       !is.null(mmAnno_clustering())){
      RNAseq_boxplot_clustering()
    }
    })
  })
  output$int_bar_clustering <- renderPlot({
    if(!is.null(input$peak_distance_clustering) && !is.null(RNAseqDEG_clustering()) && 
       !is.na(input$DEG_fdr_clustering) && !is.na(input$DEG_fc_clustering) && 
       !is.null(bw_count_clustering()) && input$Species_clustering != "not selected" && 
       !is.null(mmAnno_clustering())){
      gridExtra::grid.arrange(RNAseq_popu_clustering(), ChIPseq_popu_clustering(), ncol = 1)
    }
  })
  ChIPseq_popu_clustering <- reactive({
    RNA <- RNAseqDEG_anno_clustering()
    RNA <- dplyr::filter(RNA, !is.na(gene_id))
    RNA$gene_category <- "NS"
    RNA$gene_category[RNA$log2FoldChange > log2(input$DEG_fc_clustering) & RNA$padj < input$DEG_fdr_clustering] <- "up"
    RNA$gene_category[RNA$log2FoldChange < -log2(input$DEG_fc_clustering) & RNA$padj < input$DEG_fdr_clustering] <- "down"
    data <- merge(RNA,RP_clustering(), by="gene_id")
    data$epigenome_category <- "NS"
    data$epigenome_category[data$sumRP > 0] <- "up"
    if(length(which(unique(data$epigenome_category) == "up")) == 1){
      up <- data %>% dplyr::filter(epigenome_category == "up")
      up_total <- length(rownames(up))
      if(length(which(unique(up$gene_category) == "up")) == 1){
        up_RNAUp <- length(rownames(dplyr::filter(up, gene_category == "up")))
      }else up_RNAUp <- 0
      if(length(which(unique(up$gene_category) == "down")) == 1){
        up_RNADown <- length(rownames(dplyr::filter(up, gene_category == "down")))
      }else up_RNADown <- 0
      up_per <- c(up_total-(up_RNAUp + up_RNADown), up_RNADown, up_RNAUp)
    }else{
      up_per <- c(0,0,0)
    } 
    x <- data.frame(Regulatory_potential = c("up", "up","up"), 
                    RNA = factor(c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")")),
                                 levels = c(paste0("up:NS (",up_per[1],")"),paste0("up:down (",up_per[2],")"),paste0("up:up (",up_per[3],")"))),
                    num = up_per,
                    col = c("#F8766D","#00BFC4","grey"))
    p <- ggplot(x, aes(x = Regulatory_potential, y = num, fill = RNA)) + geom_bar(stat = "identity") +
      theme_bw(base_size = 15) + coord_flip() + scale_fill_manual(values = c("grey","#00BFC4","#F8766D"))
    return(p)
  })
  RP_all_table_clustering <- reactive({
    target_result <- regulatory_potential_clustering()$data
    target_result$epigenome_category <- "up"
    table <- NULL
    if(!is.null(mmAnno_clustering())) {
      table <- data.frame(Symbol = target_result$Symbol,
                          Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                          RNA_log2FC = -target_result$log2FoldChange,
                          RNA_padj = target_result$padj,
                          regulatory_potential = target_result$sumRP,
                          withUpPeakN = target_result$withPeakN,
                          gene_id = target_result$gene_id)
    }
    return(table)
  })
  
  RP_selected_table_clustering <- reactive({
    table <- RP_all_table_clustering() %>% dplyr::filter(Group == input$RNAseqGroup_clustering)
    return(table)
  })
  
  output$RP_table_clustering <- renderDT({
    if(!is.null(input$RNAseqGroup_clustering) && 
       !is.null(input$peak_distance_clustering && !is.null(mmAnno_clustering()))){
      RP_selected_table_clustering() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  
  output$download_clusteringintbar = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential profiling",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        gridExtra::grid.arrange(RNAseq_popu_clustering(), ChIPseq_popu_clustering(), ncol = 1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_clusteringintbox = downloadHandler(
    filename = function(){
      paste0("RNA-regulatory_potential boxplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(RNAseq_boxplot_clustering())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_clusteringKSplot = downloadHandler(
    filename = function(){
      paste0("KSplot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(regulatory_potential_clustering())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_RP_clustering_table = downloadHandler(
    filename = function() {
      paste0("RP_summary_table.txt")
    },
    content = function(file){write.table(RP_all_table_clustering(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_selected_RP_clustering_table = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup_clustering,"_table.txt")
    },
    content = function(file){write.table(RP_selected_tabl_clustering(), file, row.names = F, sep = "\t", quote = F)}
  )
  int_selected_bed_clustering <- reactive({
    if(!is.null(RP_selected_table_clustering())){
      gene <- RP_selected_table_clustering()$gene_id
      y <- NULL
      if(!is.null(mmAnno_clustering())) {
        up_peak <- subset(mmAnno_clustering(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  output$download_selected_int_peak_clustering = downloadHandler(
    filename = function() {
      paste0("RP_",input$RNAseqGroup_clustering,"_table.bed")
    },
    content = function(file){write.table(int_selected_bed_clustering(), file, row.names = F, col.names = F,sep = "\t", quote = F)}
  )
  observeEvent(input$clustering_DEG_result, ({
    updateCollapse(session,id =  "input_collapse_clustering_RP", open="KS_panel")
  }))
  
  #Clustering Integrative trackplot-------
  observeEvent(input$RP_table_clustering_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_clustering_RP", open="int_Trackplot_panel")
  }))
  output$int_igv_uprange_clustering <- renderUI({
    if(!is.null(int_goi_gene_position_clustering()) && !is.null(int_goi_promoter_position_clustering())){
      y <- int_goi_promoter_position_clustering()
      gene_position <- int_goi_gene_position_clustering()
      start_position <- min(c(y$start,gene_position$start))-1000
      end_position <- max(c(y$end,gene_position$end))+1000
      sliderInput("int_igv_uprange_clustering","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim_clustering <- renderUI({
    numericInput("int_igv_ylim_clustering","peak range:", value = 2, min = 0)
  })
  
  int_goi_promoter_position_clustering<- reactive({
    if(!is.null(input$RP_table_clustering_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table_clustering()[input$RP_table_clustering_rows_selected,]$gene_id
      y <- NULL
      if(!is.null(mmAnno_clustering())) {
        up_peak <- subset(mmAnno_clustering(), gene_id %in% gene)
        up_peak2 <- as.data.frame(up_peak)
        up_peak2$Row.names <- paste0(up_peak2$seqnames,":",up_peak2$start,"-",up_peak2$end)
        up_peak2 <- up_peak2 %>% distinct(Row.names, .keep_all = T)
        up_peak3 <- with(up_peak2,GRanges(seqnames,IRanges(start,end)))
        mcols(up_peak3) <- DataFrame(Group = "up")
        y <- as.data.frame(up_peak3)
      }
      return(y)
    }
  })
  
  int_goi_gene_position_clustering <- reactive({
    if(!is.null(input$RP_table_clustering_rows_selected)){
      gene <- RP_selected_table_clustering()[input$RP_table_clustering_rows_selected,]$gene_id
      gene_position <- as.data.frame(subset(gene_position_clustering(), gene_id %in% gene))
      return(gene_position)
    }
  })
  
  output$int_trackplot_additional_clustering <- renderUI({
    fileInput("int_trackplot_additional1_clustering",
              "Select additional bigwig files",
              accept = c("bw","BigWig"),
              multiple = TRUE,
              width = "80%")
  })
  int_track_additional_files_clustering <-reactive({
    if(!is.null(input$int_trackplot_additional1_clustering)){
      files<-c()
      name<-c()
      for(nr in 1:length(input$int_trackplot_additional1_clustering[, 1])){
        file <- input$int_trackplot_additional1_clustering[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$int_trackplot_additional1_clustering[nr,]$name))
        files <- c(files,file)
      }
      names(files)<-name
      return(files)
    }
  })
  
  int_data_track_clustering <- reactive({
    if(!is.null(input$RP_table_clustering_rows_selected)){
      return(data_trac(y=int_goi_promoter_position_clustering(),gene_position=int_goi_gene_position_clustering(),
                       gen=ref_clustering(),txdb=txdb_clustering(),org=org1_clustering(),filetype=input$data_file_type_clustering,
                       bw_files=bws_clustering(),bam_files=bam_clustering(),
                       track_additional_files=int_track_additional_files_clustering()))
    }
  })
  
  int_highlight_trackplot_clustering <- reactive({
    if(!is.null(input$RP_table_clustering_rows_selected) && !is.null(input$int_igv_uprange_clustering)){
      library(Gviz)
      y <- int_goi_promoter_position_clustering()
      gene_position <- int_goi_gene_position_clustering()
      chr <- gene_position$seqnames
      df <- int_data_track_clustering()
      start <-c()
      width <- c()
      col <- c()
      fill <- c()
      for(i in 1:length(rownames(y))){
        if(y[i,]$start < input$int_igv_uprange_clustering[2] && y[i,]$end > input$int_igv_uprange_clustering[1]){
          start <- c(start, y[i,]$start)
          width <- c(width, y[i,]$width)
          if(y[i,]$Group == "up") {
            col <- c(col,"#E41A1C")
            fill <- c(fill,"#FBB4AE")
          }
          if(y[i,]$Group == "down") {
            col <- c(col,"#377EB8")
            fill <- c(fill,"#B3CDE3")
          }
        }
      }
      if(length(start) != 0){
        ht <- HighlightTrack(trackList = df,
                             start = start, width = width,col = col,fill = fill,
                             chromosome = chr, alpha=0.5) 
      }else ht <- NULL
      return(ht)
    }
  })
  int_goi_trackplot_clustering <- reactive({
    if(!is.null(input$RP_table_clustering_rows_selected) &&
       !is.null(int_goi_promoter_position_clustering()) && 
       !is.null(int_goi_gene_position_clustering()) && 
       !is.null(gtrack_clustering()) &&
       !is.null(input$int_igv_uprange_clustering)){
      library(Gviz)
      if(!is.null(int_highlight_trackplot_clustering())){
        plot<- plotTracks(list(gtrack_clustering(), int_highlight_trackplot_clustering()),
                          from = input$int_igv_uprange_clustering[1], 
                          to = input$int_igv_uprange_clustering[2],ylim=c(0,input$int_igv_ylim_clustering),
                          type="hist")
      }else{
        df <- int_data_track_clustering()
        df[["gtrack"]] <- gtrack_clustering()
        plot<- plotTracks(df,
                          from = input$int_igv_uprange_clustering[1], 
                          to = input$int_igv_uprange_clustering[2],ylim=c(0,input$int_igv_ylim_clustering),
                          type="hist")
      }
      return(plot)
    }
  })
  output$int_trackplot_goi_clustering <- renderPlot({
    withProgress(message = "Preparing trackplot",{
      if(!is.null(input$RP_table_clustering_rows_selected)){
        int_goi_trackplot_clustering()
      }
    })
  })
  output$download_clustering_int_trackplot = downloadHandler(
    filename = "trackplot.pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        if(!is.null(int_highlight_trackplot_clustering())){
          plotTracks(list(gtrack_clustering(), int_highlight_trackplot_clustering()),
                     from = input$int_igv_uprange_clustering[1], 
                     to = input$int_igv_uprange_clustering[2],ylim=c(0,input$int_igv_ylim_clustering),
                     type="hist")
        }else{
          df <- int_data_track_clustering()
          df[["gtrack"]] <- gtrack_clustering()
          plotTracks(df,
                     from = input$int_igv_uprange_clustering[1], 
                     to = input$int_igv_uprange_clustering[2],ylim=c(0,input$int_igv_ylim_clustering),
                     type="hist")
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #Clustering integrative functional enrichment-------
  int_Hallmark_set_clustering <- reactive({
    if(!is.null(input$intGeneset_clustering)){
      return(GeneList_for_enrichment(Species = input$Species_clustering, Gene_set = input$intGeneset_clustering, org = org1_clustering()))
    }
  })
  
  output$intGroup_clustering <- renderUI({
    if(!is.null(RP_all_table_clustering())){
      selectInput("intGroup_clustering","Group (RNAseq-Epigenome)",unique(RP_all_table_clustering()$Group),multiple = T)
    }
  })
  output$intGeneset_clustering <- renderUI({
    selectInput('intGeneset_clustering', 'Gene Set', gene_set_list)
  })
  
  selected_int_group_clustering <- reactive({
    group <- input$intGroup_clustering
    df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
    colnames(df) <- c("ENTREZID","Group")
    for(name in group){
      table <- RP_all_table_clustering() %>% dplyr::filter(Group == name)
      head(table)
      df2 <- data.frame(ENTREZID = table$gene_id, Group = table$Group)
      df <- rbind(df,df2)
    }
    return(df)
  })
  int_enrich_clustering <- reactive({
    return(enrich_viewer_forMulti2(data3 = selected_int_group_clustering(), Species = input$Species_clustering, org = org1_clustering(),
                                   H_t2g = int_Hallmark_set_clustering(),Gene_set = input$intGeneset_clustering))
  })
  int_enrich_list_clustering <- reactive({
    return(enrich_gene_list(data = selected_int_group_clustering(),
                            Gene_set = input$intGeneset_clustering, org = org1_clustering(), H_t2g = int_Hallmark_set_clustering()))
  })
  int_enrich_plot_clustering <- reactive({
    return(enrich_genelist(data = selected_int_group_clustering(),
                           enrich_gene_list = int_enrich_list_clustering()))
  })
  output$int_enrichment1_clustering <- renderPlot({
    dotplot_for_output(data = int_enrich_clustering(),
                       plot_genelist = int_enrich_plot_clustering(), Gene_set = input$intGeneset_clustering, 
                       Species = input$Species_clustering)
  })
  int_enrich_table_clustering <- reactive({
    return(enrich_for_table(data = as.data.frame(int_enrich_clustering()), H_t2g = int_Hallmark_set_clustering(), Gene_set = input$intGeneset_clustering))
  })
  output$int_enrichment_result_clustering <- DT::renderDataTable({
    int_enrich_table_clustering()
  })
  output$download_clustering_int_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$intGeneset_clustering,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$clustering_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$clustering_pdf_height
        if(input$clustering_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$clustering_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = int_enrich(),
                           plot_genelist = int_enrich_plot(), Gene_set = input$intGeneset, 
                           Species = input$Species)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_clustering_int_enrichment_clustering_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$intGeneset_clustering,".txt")
    },
    content = function(file){write.table(int_enrich_table_clustering(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ## Enrichment viewer------
  output$Spe_motif_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  output$Spe_GREAT_enrich <- renderText({
    if(input$Species_enrich == "not selected") print("Please select 'Species'")
  })
  
  
  txdb_enrich <- reactive({
    if(input$Species_enrich != "not selected"){
      return(txdb_function(Species = input$Species_enrich))
    }
  })
  org_enrich <- reactive({
    return(org(Species = input$Species_enrich))
  })
  
  Enrich_peak_call_files <- reactive({
    if(is.null(input$enrich_data_file)){
      if(input$goButton_enrich > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$enrich_data_file[, 1])){
        file <- input$enrich_data_file[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$enrich_data_file[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  
  output$enrichment_input <- DT::renderDataTable({
    uploaded_files = names(Enrich_peak_call_files())
    as.data.frame(uploaded_files)
  })
  
  #Enrichment viewer functional analysis------
  output$Gene_set_enrich <- renderUI({
    selectInput('Gene_set_enrich', 'Gene Set', gene_set_list)
  })
  
  Hallmark_set_enrich <- reactive({
    if(!is.null(input$Gene_set_enrich)){
      return(GeneList_for_enrichment(Species = input$Species_enrich, Gene_set = input$Gene_set_enrich, org = org_enrich()))
    }
  })
  
  enrichment_enricher_enrich <- reactive({
    data <- Enrich_peak_call_files()
    H_t2g <- Hallmark_set_enrich()
    df <- list()
    source <- ref_for_GREAT(input$Species_enrich)
    for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g,input$Species_enrich),source)
      df[[name]] <- res
    }
    return(df)
  })
  
  enrichment_1_1_enrich <- reactive({
    df <- enrichment_enricher_enrich()
    if(!is.null(input$Gene_set_enrich) && input$Species_enrich != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                          "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
        }else group1 <- NULL
        group1$Group <- name
        data <- rbind(data, group1)
      }
      if(length(data$id) != 0) {
        return(data)
      }else return(NULL)
    }else{return(NULL)}
  })
  
  pair_enrich1_H_enrich <- reactive({
    if(!is.null(input$Gene_set_enrich) && input$Species_enrich != "not selected"){
      df <- enrichment_enricher_enrich()
      data3 <- Enrich_peak_call_files()
      data <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getEnrichmentTable(df[[name]]))
          group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
          if (length(group1$p_adjust_hyper) > input$enrich_showCategory){
            group1 <- group1[sort(group1$p_adjust_hyper, decreasing = F, index=T)$ix,]
            group1 <- group1[1:input$enrich_showCategory,]
          }
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      colnames(data) <-  colnames(data) <- c("Description", "genome_fraction", "observed_region_hits", "fold_enrichment", "p_value", "p_adjust", "mean_tss_dist", 
                                             "observed_gene_hits", "gene_set_size","fold_enrichment_hyper","p_value_hyper","p_adjust_hyper","Group")
      if(length(data$Description) != 0){
        data$Description <- gsub("_", " ", data$Description)
        data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "p_adjust_hyper"))))))
        data$x <- gsub(":","", data$x)
        data <- dplyr::arrange(data, x)
        idx <- order(data[["x"]], decreasing = FALSE)
        data$Description <- factor(data$Description,
                                   levels=rev(unique(data$Description[idx])))
        p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                        geom_point() +
                        scale_color_continuous(low="red", high="blue",
                                               guide=guide_colorbar(reverse=TRUE)) +
                        scale_size(range=c(1, 6))+ theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) + 
                        scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      }else p1 <- NULL
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment_enrich <- renderPlot({
    if(input$Species_enrich != "not selected" && !is.null(Enrich_peak_call_files())){
      dotplot_for_output(data = Enrich_peak_call_files(),
                         plot_genelist = pair_enrich1_H_enrich(), Gene_set = input$Gene_set_enrich, 
                         Species = input$Species_enrich)
    }
  })
  
  
  output$download_enrich_enrichment = downloadHandler(
    filename = function(){
      paste0("Enrichment-",input$Gene_set_enrich,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H_enrich()
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_enrich_enrichment_table = downloadHandler(
    filename = function() {
      paste0("Enrichment_table-",input$Gene_set_enrich,".txt")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1_enrich()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  pair_enrich_table_enrich <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1_enrich()), 
                            H_t2g = Hallmark_set_enrich(), Gene_set = input$Gene_set_enrich))
  })
  
  output$enrich_enrichment_result <- DT::renderDataTable({
    if(input$Species_enrich != "not selected" && !is.null(Enrich_peak_call_files())){
      as.data.frame(enrichment_1_1_enrich())
    }
  })
  
  
  
  output$whichGeneSet_enrich <- renderUI({
    if(!is.null(input$intersection_enrich_fun)){
      group <- input$intersection_enrich_fun
      set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1_enrich()), Group == group)$id)
      selectInput('Pathway_list_enrich', 'Pathway list', set_list)
    }
  })
  
  region_gene_associate_enrich <- reactive({
    if(!is.null(input$Pathway_list_enrich)){
      set_list <- input$Pathway_list_enrich
      df <- enrichment_enricher_enrich()
      data3 <- Enrich_peak_call_files()
      data <- data.frame(matrix(rep(NA, 9), nrow=1))[numeric(0), ]
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(rGREAT::getRegionGeneAssociations(df[[name]], term_id = set_list))
          group1$Group <- paste(name, "\n(",length(as.data.frame(data3[[name]])$start),")",sep = "")
        }else group1 <- NULL
        data <- rbind(data, group1)
      }
      data <- apply(data,2,as.character)
      return(data)
    }
  })
  output$region_gene_enrich_associations <- DT::renderDT({
    if(!is.null(input$Pathway_list_enrich) && 
       !is.null(input$intersection_enrich_fun) && input$Species_enrich != "not selected"){
      region_gene_associate_enrich()
    }
  })
  output$whichGenes_enrich <- renderUI({
    if(input$Species_enrich != "not selected"){ 
      selectInput('intersection_enrich_fun', 'Select intersection', names(enrichment_enricher_enrich()))
    }
  })
  
  region_gene_associate_plot_enrich <- reactive({
    set_list <- input$Pathway_list_enrich
    df <- enrichment_enricher_enrich()
    name <- input$intersection_enrich_fun
    if(!is.null(df[[name]])) {
      res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
    }else res <- NULL
    return(res)
  })
  output$region_gene_associations_enrich_plot <- renderPlot({
    if(!is.null(input$Pathway_list_enrich) && 
       !is.null(input$intersection_enrich_fun) && input$Species_enrich != "not selected"){
      region_gene_associate_plot_enrich()
    }
  })
  output$download_enrich_region_gene_associations_plot = downloadHandler(
    filename = function(){
      paste0("region_gene_associations-",input$intersection_enrich_fun,"-",input$Pathway_list_enrich,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$enrich_pdf_width
        set_list <- input$Pathway_list_enrich
        df <- enrichment_enricher_enrich()
        name <- input$intersection_enrich_fun
        pdf(file, height = pdf_height, width = pdf_width)
        rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_enrich_region_gene_associations = downloadHandler(
    filename = function() {
      paste0("region_gene_associations-",input$intersection_enrich_fun,"-",input$Pathway_list_enrich,".txt")
    },
    content = function(file){write.table(region_gene_associate_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  ##Enrich motif -----------
  output$homer_size_enrich <- renderUI({
    radioButtons('homer_size_enrich','Type of the region for motif finding',
                 c('given (exact size)'="given",
                   'custom size'="custom"
                 ),selected = "custom")
  })
  output$homer_size2_enrich <- renderUI({
    if(!is.null(input$homer_size_enrich)){
    if(input$homer_size_enrich == "custom"){
      numericInput('homer_size2_enrich','Size of the region for motif finding',value=200, step=100)
      }}
  })
  output$homer_bg_enrich <- renderUI({
    radioButtons('homer_bg_enrich','Background sequence',
                 c('random'="random",
                   'bed files'="peakcalling"
                 ),selected = "random")
  })
  output$homer_bg2_enrich <- renderUI({
    if(!is.null(input$homer_bg_enrich)){
      if(input$homer_bg_enrich == "peakcalling"){
        fileInput('homer_bg2_enrich',
                  'Select bed files',
                  accept = c("bed","narrowPeak"),
                  multiple = TRUE,
                  width = "80%")
      }}
  })
  Enrich_peak_call_files_homer <- reactive({
    if(is.null(input$homer_bg2_enrich)){
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$homer_bg2_enrich[, 1])){
        file <- input$homer_bg2_enrich[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$homer_bg2_enrich[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      files3 <- soGGi:::runConsensusRegions(GRangesList(files2), "none")
      return(files3)
    }
  })
  output$homer_unknown_enrich <- renderUI({
    selectInput("homer_unknown_enrich","Type of enrichment analysis",c("known motif","known and de novo motifs"), selected = "known motif")
  })
  observeEvent(input$homer_unknown_enrich,({
    updateActionButton(session,"motifButton_enrich", "Start")
  }))

  enrich_motif_enrich <- reactive({
    if(input$homer_size_enrich == "given") size <- "given"
    if(input$homer_size_enrich == "custom") size <- input$homer_size2_enrich
    if(input$motifButton_enrich > 0 && !is.null(Enrich_peak_call_files()) 
       && input$Species_enrich != "not selected" && !is.null(input$homer_unknown_enrich)){
      return(findMotif(df= Enrich_peak_call_files(), Species = input$Species_enrich,size=size,back = input$homer_bg_enrich,
                       motif=input$homer_unknown_enrich, other_data = Enrich_peak_call_files_homer(),type="Other"))
    }
  })

  output$motif_enrich_plot <- renderPlot({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && 
       !is.null(input$homer_unknown_enrich) && input$Species_enrich != "not selected"){
      homer_Motifplot(df = enrich_motif_enrich(),showCategory = input$enrich_showCategory)
    }
  })
  
  output$download_motif_enrich_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- homer_Motifplot(df = enrich_motif_enrich(),showCategory = input$enrich_showCategory)
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  

  output$download_motif_enrich_table = downloadHandler(
    filename = function() {
      paste0("known_motif_table",".txt")
    },
    content = function(file){write.table(motif_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  

  denovo_motif_table_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      return(denovo_motif(enrich_motif_enrich()))
    }
  })
  
  output$denovo_motif_enrich_result <- DT::renderDT({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && input$Species_enrich != "not selected"){
      denovo_motif_table_enrich()
    }
  })
  
  motif_table_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      return(known_motif(enrich_motif_enrich()))
    }
  })
  
  output$motif_enrich_result <- DT::renderDT({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich()) && input$Species_enrich != "not selected"){
      motif_table_enrich()
    }
  })
  
  output$download_denovo_motif_enrich_table = downloadHandler(
    filename = function() {
      paste0("denovo_motif_table",".txt")
    },
    content = function(file){write.table(denovo_motif_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_homer_report_enrich = downloadHandler(
    filename = function() {
      paste0("HOMER_report",".zip")
    },
    content = function(fname){
      fs <- c()
      path_list <- enrich_motif_enrich()
      base_dir <- gsub("\\/.+$", "", path_list[[names(path_list)[1]]])
      for(name in names(path_list)){
        files <-list.files(path_list[[name]],pattern = "*.*")
        for(i in 1:length(files)){
          data <- paste0(path_list[[name]],"/",files[[i]])
          fs <- c(fs, data)
        }
      }
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip"
  )
  
  ##More-------------
  #Bed tool----------
  output$bed_file1_warning <- renderText({
    if(is.null(bed_peak_call_files1())) print("Please upload bed files (A)")
  })
  output$bed_file2_warning <- renderText({
    if(input$data_file_type_bed != "type1" && 
       is.null(bed_peak_call_files2())) print("Please upload bed files (B)")
  })
  
  bed_peak_call_files1 <- reactive({
    if(is.null(input$bed_data_file)){
      if(input$goButton_bed > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$bed_data_file[, 1])){
        file <- input$bed_data_file[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$bed_data_file[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  bed_peak_call_files2 <- reactive({
    if(is.null(input$bed_data_file2)){
      if(input$goButton_bed > 0 ){
        files <- list()
        files[["A_1"]] <- "data/peakcall/A_1_peaks.narrowPeak"
        files[["A_2"]] <- "data/peakcall/A_2_peaks.narrowPeak"
        files[["B_1"]] <- "data/peakcall/B_1_peaks.narrowPeak"
        files[["B_2"]] <- "data/peakcall/B_2_peaks.narrowPeak"
        files2 <- lapply(files, GetGRanges, simple = TRUE)
        return(files2)
      }
      return(NULL)
    }else{
      files<-c()
      name<-c()
      for(nr in 1:length(input$bed_data_file2[, 1])){
        file <- input$bed_data_file2[[nr, 'datapath']]
        name <- c(name, gsub("\\..+$", "", input$bed_data_file2[nr,]$name))
        files <- c(files,file)
      }
      files2 <- lapply(files, GetGRanges, simple = TRUE)
      names(files2)<-name
      return(files2)
    }
  })
  bed_input<-reactive({
    if(input$data_file_type_bed == "type1"){
      uploaded_files = names(bed_peak_call_files1())
    }else if(input$data_file_type_bed == "type2"){
      uploaded_files <- c()
      for(i in 1:length(names(bed_peak_call_files1()))){
        uploaded_files <- c(uploaded_files, paste0(names(bed_peak_call_files1()[i])," - ", 
                                                   names(bed_peak_call_files2()[i])))
      }
    }else {
      if(input$intersect_bed == "default") mode <- "_intersect_"
      if(input$intersect_bed == "wa") mode <- "_intersect-wa_"
      if(input$intersect_bed == "exclude") mode <- "_intersect-v_"
      uploaded_files <- c()
      for(i in 1:length(names(bed_peak_call_files1()))){
        uploaded_files <- c(uploaded_files, paste0(names(bed_peak_call_files1()[i]),mode, 
                                                   names(bed_peak_call_files2()[i])))
      }
    }
    return(as.data.frame(uploaded_files))
  })
  
  output$bed_input <- DT::renderDataTable({
    bed_input()
  })
  
  output$bed_merge_dist <- renderUI({
    numericInput("bed_merge_dist","Gap size (bp)",value = 1000,min = 1)
  })
  
  interval_merge <- reactive({
    if(input$data_file_type_bed == "type1"){
    files <- bed_peak_call_files1()
    df <- list()
    for(name in names(files)){
      file <- files[[name]]
      df[[name]] <- merge_bed(file, max_dist = input$bed_merge_dist)
    }
    }else if(input$data_file_type_bed == "type2"){
      files <- bed_peak_call_files1()
      files2 <- bed_peak_call_files2()
      df <- list()
      for(i in 1:length(names(files))){
        name <- paste0(names(bed_peak_call_files1()[i])," - ", 
                              names(bed_peak_call_files2()[i]))
        df[[name]] <- subtract_bed(files[[i]],files2[[i]])
      }
    }else{
      if(input$intersect_bed == "default") mode <- "_intersect_"
      if(input$intersect_bed == "wa") mode <- "_intersect-wa_"
      if(input$intersect_bed == "exclude") mode <- "_intersect-v_"
      files <- bed_peak_call_files1()
      files2 <- bed_peak_call_files2()
      df <- list()
      for(i in 1:length(names(files))){
        name <- paste0(names(bed_peak_call_files1()[i]),mode, 
                       names(bed_peak_call_files2()[i]))
        df[[name]] <- intersect_bed(files[[i]],files2[[i]],mode=input$intersect_bed)
      }
    }
    return(df)
  })
  
  output$download_bed <- downloadHandler(
    filename = function() {
      paste("bed_files", "zip", sep=".")
    },
    content = function(fname) {
      withProgress(message = "Preparing download, please wait",{
        fs <- c()
        files <- interval_merge()
        for(name in names(files)){
          file_name <- paste0(name, ".bed")
          fs <- c(fs, file_name)
          write.table(as.data.frame(files[[name]]), file_name,quote = F, row.names = F, col.names = F,sep = "\t")
        }
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  
})