popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=10000000*1024^2)
  # pair-wise ------------------------------------------------------------------------------
  output$Spe <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  org1 <- reactive({
    return(org(Species = input$Species))
  })
  txdb <- reactive({
    if(input$Species != "not selected"){
    return(txdb_function(Species = input$Species))
    }
  })
  promoter_region <- reactive({
    if(input$Species != "not selected"){
      if(input$Genomic_region == "Promoter"){
        return(promoter(txdb(),upstream = input$upstream, downstream = input$downstream))
      }else return(promoter(txdb(),upstream = input$upstream, downstream = input$downstream,
                            input_type = "Genome-wide",files =peak_call_files()))
    }
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
          df[["Ctrl_1"]] <- "data/bigwig/Growing-p53_1.BigWig"
          df[["Ctrl_2"]] <- "data/bigwig/Growing-p53_2.BigWig"
          df[["Sen_1"]] <- "data/bigwig/OIS-p53_1.BigWig"
          df[["Sen_2"]] <- "data/bigwig/OIS-p53_2.BigWig"
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
  
  peak_call_files <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      if(input$data_file_type == "Row1"){
        if(is.null(input$peak_call_file1)){
          if(input$goButton > 0 ){
            df <- list()
            df[["Ctrl_1"]] <- "data/peakcall/Growing-p53_1_peaks.narrowPeak"
            df[["Ctrl_2"]] <- "data/peakcall/Growing-p53_2_peaks.narrowPeak"
            df[["Sen_1"]] <- "data/peakcall/OIS-p53_1_peaks.narrowPeak"
            df[["Sen_2"]] <- "data/peakcall/OIS-p53_2_peaks.narrowPeak"
            name <- c("Ctrl_1","Ctrl_2","Sen_1","Sen_2")
            files2 <- lapply(df, ChIPQC:::GetGRanges, simple = TRUE)
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
          files2 <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
          names(files2)<-name
          print(files2)
          return(files2)
        }
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
    }
  })
  
  output$input_peak_call_files <- DT::renderDataTable({
    if(input$Genomic_region == "Genome-wide"){
    if(input$data_file_type == "Row1"){
      uploaded_files = names(peak_call_files())
      as.data.frame(uploaded_files)
    }
    }
  })
  
  bw_count <- reactive({
    if(input$Species != "not selected" && input$data_file_type == "Row1" && 
       !is.null(bw_files())){
      if(input$Genomic_region == "Promoter"){
        return(Bigwig2count(bw = bw_files(),promoter_region(),
                            Species = input$Species,input_type =input$Genomic_region))
      }else{
        if(!is.null(peak_call_files())){
          return(Bigwig2count(bw = bw_files(),promoter_region(),
                              Species = input$Species,input_type =input$Genomic_region))
        }
      }
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
  
  deg_result <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      count <- bw_count()
      collist <- gsub("\\_.+$", "", colnames(count))
      dds <- dds()
      contrast <- c("con", unique(collist))
      res <- results(dds,  contrast = contrast)
      if(input$FDR_method == "IHW") {
        ihw_res <- ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
        res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
      }
      if(input$FDR_method == "Qvalue") {
        res <- results(dds,  contrast = contrast)
        qvalue <- qvalue::qvalue(res$pvalue)
        res$padj <- qvalue$qvalues
      }
      res <- as.data.frame(res)
      return(res)
    }
  })
  
  deg_norm_count <- reactive({
    if(is.null(bw_count())){
      return(NULL)
    }else{
      count <- bw_count()
      collist <- gsub("\\_.+$", "", colnames(count))
      group <- data.frame(con = factor(collist))
      dds <- dds()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
      return(normalized_counts)
    }
  })
  
  observeEvent(bws(), ({
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
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
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
    gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                    keytype = "ENTREZID",
                                    columns = c("SYMBOL", "ENTREZID"))
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
    if(input$Genomic_region == "Promoter"){
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
    if(is.null(bw_files()) || input$Species == "not selected"){
      return(NULL)
    }else{
      print(PCAplot(data = deg_norm_count(),plot=TRUE))
    }
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
    sliderInput("xrange","X_axis range:",min = -100,
                max=100, step = 1,
                value = c(-5, 5))
  })
  output$volcano_y <- renderUI({
    sliderInput("yrange","Y_axis range:",min = 0, max= 300, step = 1,
                value = 10)
  })
  
  pair_volcano <- reactive({
    if(!is.null(input$xrange)){
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
          Color <- c("blue","green","darkgray","red")
          for(name in label_data){
            data$color[data$Row.names == name] <- "GOI"
          }
        }else{
          Color <- c("blue","darkgray","red")
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
    if(!is.null(input$xrange)){
      if(is.null(bw_count())){
        return(NULL)
      }else{
        
        print(pair_volcano())
      }
    }
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
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(volcano_heatmap())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  ##trackplot----------
  ref <- reactive({
    switch (input$Species,
            "Mus musculus (mm10)" = ref <- "mm10",
            "Homo sapiens (hg19)" = ref <- "hg19",
            "Homo sapiens (hg38)" = ref <- "hg38")
    return(ref)
  })
  
  
  
  observeEvent(input$DEG_result_rows_selected, ({
    updateCollapse(session,id =  "input_collapse_pair_DEG", open="Trackplot_panel")
  }))
  observeEvent(input$DEG_result_rows_selected,({
    if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
    y <- goi_promoter_position()
    gene_position <- goi_gene_position()
    start_position <- min(c(y$start,gene_position$start))
    end_position <- max(c(y$end,gene_position$end))
    updateSliderInput(session,"igv_uprange","Range:",
                       value = c(start_position,end_position),
                      step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  }))
  
  output$igv_uprange <- renderUI({
    if(!is.null(goi_gene_position()) && !is.null(goi_promoter_position())){
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
    sliderInput("igv_uprange","Range:",value = c(start_position,end_position),
                step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_ylim <- renderUI({
    numericInput("igv_ylim","peak range:", value = 10, min = 0)
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
        gene_IDs<-AnnotationDbi::select(org1(),keys = label_data,
                                        keytype = "SYMBOL",
                                        columns = "ENTREZID")
        y <- as.data.frame(subset(promoter_region(), gene_id %in% gene_IDs))
      }else{
        label_data <- rownames(deg_result2()[input$DEG_result_rows_selected,])
        print(label_data)
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
      gene_IDs<-AnnotationDbi::select(org1(),keys = label_data,
                                      keytype = "SYMBOL",
                                      columns = "ENTREZID")
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
        print(name)
        files <- c(files,file)
      }
      names(files)<-name
      print(input$trackplot_additional1[[1, 'datapath']])
      print(files)
      return(files)
    }
  })
  
  data_track <- reactive({
    if(!is.null(input$DEG_result_rows_selected)){
      library(Gviz)
      y <- goi_promoter_position()
      gene_position <- goi_gene_position()
      chr <- gene_position$seqnames
      gen <- ref()
      txdb <- txdb()
      grtrack <- GeneRegionTrack(txdb,
                                 chromosome = chr, name = "UCSC known genes",
                                 transcriptAnnotation = "tx_name",genome = gen,
                                 background.title = "grey",cex = 1.25)
      bw_files <- bws()
      cond1 <- unique(gsub("\\_.+$", "", names(bw_files)))[1]
      cond2 <- unique(gsub("\\_.+$", "", names(bw_files)))[2]
      if(!is.null(track_additional_files())) bw_files <- c(bw_files, track_additional_files())
      df <- list()
      for(name in names(bw_files)){
        cond <- gsub("\\_.+$", "", name)
        if(cond == cond1) col = "black"
        if(cond == cond2) col = "darkred"
        if(cond != cond1 && cond != cond2) col = "darkblue"
        df[[name]] <- DataTrack(range = bw_files[[name]], type = "l",genome = gen,
                                name = gsub("\\..+$", "", name), window = -1,
                                chromosome = chr, background.title = col,cex = 1.25,
                                col.histogram = "darkblue", 
                                fill.histogram = "darkblue",)
      }
      df[["grtrack"]] <- grtrack
      return(df)
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
                             chromosome = chr)
      })
      }else ht <- NULL
      return(ht)
    }
  })
  goi_trackplot <- reactive({
    if(!is.null(input$DEG_result_rows_selected) &&
       !is.null(goi_promoter_position()) && 
       !is.null(goi_gene_position()) && 
       !is.null(gtrack()) &&
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
      if(!is.null(input$DEG_result_rows_selected)){
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
        plotTracks(list(gtrack(), highlight_trackplot()),
                   from = goi_gene_position()$start-input$igv_uprange, 
                   to = goi_gene_position()$end+input$igv_downrange,ylim=c(0,input$igv_ylim),
                   type="hist")
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
        switch (input$Species,
                "Homo sapiens (hg38)" = source <- "TxDb.Hsapiens.UCSC.hg38.knownGene",
                "Homo sapiens (hg19)" = source <- "TxDb.Hsapiens.UCSC.hg19.knownGene",
                "Mus musculus (mm10)" = source <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
        )
        data_degcount_up2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_up()))
        data_degcount_up3 <- with(data_degcount_up2, GRanges(seqnames = seqnames, 
                                                             ranges = IRanges(start,end),
                                                             ENTREZID = geneId))
        res_up = rGREAT::great(gr = data_degcount_up3,gene_sets = gene_list_for_enrichment_genome(H_t2g),source)
        data_degcount_down2 <- dplyr::filter(data3, locus %in% rownames(data_degcount_down()))
        data_degcount_down3 <- with(data_degcount_down2, GRanges(seqnames = seqnames, 
                                                                 ranges = IRanges(start,end),
                                                                 ENTREZID = geneId))
        res_down = rGREAT::great(data_degcount_down3,gene_list_for_enrichment_genome(H_t2g), source)
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
                            scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
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
            if(!is.null(df[[name]])) {
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
              print(data)
              idx <- order(data[["x"]], decreasing = FALSE)
              data$Description <- factor(data$Description,
                                         levels=rev(unique(data$Description[idx])))
              p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                              geom_point() +
                              scale_color_continuous(low="red", high="blue",
                                                     guide=guide_colorbar(reverse=TRUE)) +
                              scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
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
          pdf_width <- 12
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
    if(input$Genomic_region == "Genome-wide" && !is.null(input$Up_or_Down)){
      group <- input$Up_or_Down
    set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1()), Group == group)$id)
    print(set_list)
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
      if(!is.null(df[[name]])) {
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
       !is.null(input$Up_or_Down)){
    set_list <- input$Pathway_list
    df <- enrichment_enricher()
    name <- input$Up_or_Down
    if(!is.null(df[[name]])) {
      res <- rGREAT::plotRegionGeneAssociations(df[[name]], term_id = set_list)
    }else res <- NULL
    return(res)
    }
  })
  output$region_gene_associations_plot <- renderPlot({
    region_gene_associate_plot()
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
  preMotif_list <- reactive({
    df <- list()
    df[["Up"]] <- data_degcount_up()
    df[["Down"]] <- data_degcount_down()
    return(df)
  })
  pwms_motif <- reactive({
    if(input$motifButton > 0 && !is.null(preMotif_list())){
    return(pwms(Species = input$Species))
    }
  })
  
  enrich_motif <- reactive({
    if(input$motifButton > 0 && !is.null(pwms_motif()) && 
       !is.null(preMotif_list())){
      if(input$Genomic_region == "Genome-wide"){
      return(MotifAnalysis(df= preMotif_list(), anno_data = deg_result_anno2(),
                           Species = input$Species, pwms = pwms_motif(), consensus = promoter_region()))
      }else return(MotifAnalysis(df= data_degcount2(), anno_data = promoter_region(),
                                 Species = input$Species, pwms = pwms_motif(),
                                 type = "Promoter"))
    }
  })
  
  output$motif_plot <- renderPlot({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      Motifplot(df2 = enrich_motif(),padj = 0.05)
    }
  })
  output$download_motif_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- Motifplot(df2 = enrich_motif(),padj = 0.05)
        if(input$pair_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_motif_table = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_table",".txt")
    },
    content = function(file){write.table(motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  motif_table <- reactive({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      df2 <- enrich_motif()
      df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
      for(name in names(df2)){
        res <- df2[[name]]
        res <- dplyr::filter(res, X1 > -log10(0.05))
        res <- res %>% dplyr::arrange(-X1.1)
        df <- rbind(df, res)
      }
      colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                        "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                        "Group")
      df$padj <- 10^(-df$negLog10Padj)
      return(df)
    }
  })

  output$motif_result <- DT::renderDT({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      motif_table()
    }
  })
  
  promoter_motif_region <- reactive({
    target_motif <- motif_table()[input$motif_result_rows_selected,]
    if(!is.null(input$motif_result_rows_selected)){
      if(input$Genomic_region == "Genome-wide"){
      res <- MotifRegion(df= preMotif_list(), anno_data = deg_result_anno2(),
                         target_motif = target_motif,Species = input$Species)  
      }else{
        res <- MotifRegion(df= data_degcount2(), anno_data = promoter_region(),
                           target_motif = target_motif,Species = input$Species,
                           type = "Promoter")
        res <- res[,-1]
      }
      return(res)
    }
  })
  
  output$promoter_motif_region_table <- renderDataTable({
    target_motif <- motif_table()[input$motif_result_rows_selected,]
    if(!is.null(input$motif_result_rows_selected)){
      promoter_motif_region()
    }
  })
  
  output$download_promoter_motif_region = downloadHandler(
    filename = function() {
      paste0("motif_region",".txt")
    },
    content = function(file){write.table(promoter_motif_region(), file, row.names = F, sep = "\t", quote = F)}
  )
  observeEvent(motif_table()[input$motif_result_rows_selected,], ({
    updateCollapse(session,id =  "Promoter_motif_collapse_panel", open="Promoter_motif_region_panel")
  }))
  
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
        files[["data/peakcall/Growing-p53_1_peaks.narrowPeak"]] <- "data/peakcall/Growing-p53_1_peaks.narrowPeak"
        files[["data/peakcall/Growing-p53_2_peaks.narrowPeak"]] <- "data/peakcall/Growing-p53_2_peaks.narrowPeak"
        files[["data/peakcall/OIS-p53_1_peaks.narrowPeak"]] <- "data/peakcall/OIS-p53_1_peaks.narrowPeak"
        files[["data/peakcall/OIS-p53_2_peaks.narrowPeak"]] <- "data/peakcall/OIS-p53_2_peaks.narrowPeak"
      }
      print(files)
      Glist <- files2GRangelist(files)
      Glist[["Consensus_region"]] <- promoter_region()
      return(Glist)
    }
    }
  })
  output$input_peak_distribution <- renderPlot({
    withProgress(message = "Plot peak distribution for peak call files",{
    if(!is.null(input_peak_list()) && !is.null(txdb()) && input$Genomic_region == "Genome-wide"){
      print(input_peak_list())
      genomicElementDistribution(input_peak_list(), 
                                 TxDb = txdb()) 
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
      if(!is.null(deg_peak_list()) && input$Genomic_region == "Genome-wide"){
        genomicElementDistribution(deg_peak_list(), 
                                   TxDb = txdb()) 
      }
    })
  })

  
  output$download_input_peak_distribution = downloadHandler(
    filename = function() {
      paste0("input_peak_distribution",".pdf")
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
        genomicElementDistribution(input_peak_list(), 
                                   TxDb = txdb()) 
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_deg_peak_distribution = downloadHandler(
    filename = function() {
      paste0("DAR_peak_distribution",".pdf")
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
        genomicElementDistribution(deg_peak_list(), 
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
      print(name)
      files <- c(files,file)
    }
    names(files)<-name
    print(input$peak_pattern_up_add[[1, 'datapath']])
    print(files)
    return(files)
    }
  })
  output$peak_pattern_up_heat_range <- renderUI({
    numericInput("peak_up_range","Heatmap y-axis range",value=10,step=2)
  })
peak_up_grange <- reactive({
  if(input$Genomic_region == "Genome-wide"){
  up <- range_changer(data_degcount_up())
  up2 <- with(up, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
  }else{
    up <- symbol2gene_id(data_degcount_up(),org1()) %>% distinct(gene_id, .keep_all = T)
    up <- subset(promoter_region(), gene_id %in% up$gene_id) 
  }
  return(up)
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
    feature.recentered <- reCenterPeaks(up2, width=4000)
    cvglists <- peak_up_cvglists()
    feature.center <- reCenterPeaks(up2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000)
    return(sig)
  })
  
  peak_up_feature_center <- reactive({
    up2 <- peak_up_grange()
    feature.center <- reCenterPeaks(up2, width=1)
    sig <- peak_up_sig()
    keep <- rowSums(sig[[2]]) > 0
    feature.center <- feature.center[keep]
    return(feature.center)
  })
  peak_up_alinedHeatmap <- reactive({
    if(!is.null(input$peak_up_range)){
    range <- c()
    for(n in length(names(bws()))) {range <- c(range, input$peak_up_range)}
    heatmap <- featureAlignedHeatmap(peak_up_sig(), peak_up_feature_center(),
                                     upstream=2000, downstream=2000,
                                     upper.extreme=range,color = c("white","red"))
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
    numericInput("peak_down_range","Heatmap y-axis range",value=10,step=2)
  })
  peak_down_grange <- reactive({
    if(input$Genomic_region == "Genome-wide"){
      down <- range_changer(data_degcount_down())
      down2 <- with(down, GRanges(seqnames = chr,ranges = IRanges(start=start,end=end)))
    }else{
      down <- symbol2gene_id(data_degcount_down(),org1()) %>% distinct(gene_id, .keep_all = T)
      down <- subset(promoter_region(), gene_id %in% down$gene_id)
    }
    return(down)
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
    feature.recentered <- reCenterPeaks(down2, width=4000)
    cvglists <- peak_down_cvglists()
    feature.center <- reCenterPeaks(down2, width=1)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000)
    return(sig)
  })
  
  peak_down_feature_center <- reactive({
    down2 <- peak_down_grange()
    feature.center <- reCenterPeaks(down2, width=1)
    sig <- peak_down_sig()
    keep <- rowSums(sig[[2]]) > 0
    feature.center <- feature.center[keep]
    return(feature.center)
  })
  peak_down_alinedHeatmap <- reactive({
    if(!is.null(input$peak_down_range)){
      range <- c()
      for(n in length(names(bws()))) {range <- c(range, input$peak_down_range)}
      heatmap <- featureAlignedHeatmap(peak_down_sig(), peak_down_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = c("white","red"))
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
    RNAseqDEG()
  })
  
  RNAseqDEG_anno <- reactive({
    RNAdata <- RNAseqDEG()
    RNAdata$log2FoldChange <- -RNAdata$log2FoldChange
    if(str_detect(rownames(RNAdata)[1], "ENS")){
      my.symbols <- gsub("\\..*","", rownames(RNAdata))
      gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                      keytype = "ENSEMBL",
                                      columns = c("ENSEMBL","SYMBOL","ENTREZID"))
      colnames(gene_IDs) <- c("EnsemblID","Symbol","gene_id")
      RNAdata$EnsemblID <- gsub("\\..*","", rownames(RNAdata))
      gene_IDs <- gene_IDs %>% distinct(EnsemblID, .keep_all = T)
    }else{
      my.symbols <- rownames(RNAdata)
      gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                      keytype = "SYMBOL",
                                      columns = c("SYMBOL", "ENTREZID"))
      colnames(gene_IDs) <- c("Symbol", "gene_id")
      gene_IDs <- gene_IDs %>% distinct(Symbol, .keep_all = T)
      RNAdata$Symbol <- rownames(RNAdata) 
    }
    data <- merge(RNAdata, gene_IDs, by="Symbol")
    return(data)
  })
  
  mmAnno_up <- reactive({
    up_peak <- peak_up_grange()
    if(input$Genomic_region == "Promoter"){ 
      up_peak <- up_peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      up_peak <- with(up_peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
    }
    range <- input$peak_distance * 1000
    mcols(up_peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(up_peak))))
    mmAnno_up <- mm_geneScan(up_peak, txdb(),upstream = range,downstream = range)
    return(mmAnno_up)
  })
  mmAnno_down <- reactive({
    down_peak <- peak_down_grange()
    if(input$Genomic_region == "Promoter"){ 
      down_peak <- down_peak  %>% as.data.frame() %>% distinct(start, .keep_all = T)
      down_peak <- with(down_peak, GRanges(seqnames = seqnames,ranges = IRanges(start=start,end=end)))
    }
    range <- input$peak_distance * 1000
    mcols(down_peak) <- DataFrame(feature_id = paste0("peak_", seq_len(length(down_peak))))
    mmAnno_down <- mm_geneScan(down_peak, txdb(),upstream = range,downstream = range)
    return(mmAnno_down)
  })
  
  regulatory_potential <- reactive({
    mmAnno_up <- mmAnno_up()
    print(mmAnno_up)
    result_geneRP_up <- calcRP_TFHit(mmAnno = mmAnno_up,Txdb = txdb())
    mmAnno_down <- mmAnno_down()
    result_geneRP_down <- calcRP_TFHit(mmAnno = mmAnno_down,Txdb = txdb())
    colnames(result_geneRP_up)[2:4] <- paste0(colnames(result_geneRP_up)[2:4], "_up")
    colnames(result_geneRP_down)[2:4] <- paste0(colnames(result_geneRP_down)[2:4], "_down")
    result_geneRP <- merge(result_geneRP_up,result_geneRP_down,by="gene_id",all =T)
    result_geneRP$sumRP <- apply(data.frame(result_geneRP$sumRP_up,
                                            -result_geneRP$sumRP_down),1,sum,na.rm=TRUE)
    result_geneRP <- result_geneRP %>% dplyr::arrange(-sumRP)
    result_geneRP$RP_rank <- rownames(result_geneRP) %>% as.numeric()
    merge_data <- integrate_ChIP_RNA(
      result_geneRP = result_geneRP,
      result_geneDiff = RNAseqDEG_anno(),lfc_threshold = log2(input$DEG_fc),padj_threshold = input$DEG_fdr
    )
    return(merge_data)
  })
  
  output$ks_plot <- renderPlot({
    if(!is.null(input$peak_distance)){
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
                 min=0,max=1000,value=100,step = 10)
  })
  output$RNAseqGroup <- renderUI({
    if(!is.null(RP_all_table())){
      selectInput("RNAseqGroup","Group (RNAseq-Epigenome)",
                  unique(RP_all_table()$Group),
                  multiple = FALSE) 
    }
  })
  RP_all_table <- reactive({
    target_result <- regulatory_potential()$data
    target_result$epigenome_category <- "up"
    target_result$epigenome_category[target_result$sumRP < 0] <- "down"
    table <- data.frame(Symbol = target_result$Symbol,
                        Group = paste0(target_result$gene_category,"-",target_result$epigenome_category),
                        RNA_log2FC = -target_result$log2FoldChange,
                        RNA_padj = target_result$padj,
                        regulatory_potential = target_result$sumRP,
                        withUpPeakN = target_result$withPeakN_up,
                        withDownPeakN = target_result$withPeakN_down,
                        gene_id = target_result$gene_id)
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
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        regulatory_potential()
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
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
      sliderInput("int_igv_uprange","Range:",value = c(start_position,end_position),
                  step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$int_igv_ylim <- renderUI({
    numericInput("int_igv_ylim","peak range:", value = 10, min = 0)
  })
  
  int_goi_promoter_position<- reactive({
    if(!is.null(input$RP_table_rows_selected)){
      library(Gviz)
      gene <- RP_selected_table()[input$RP_table_rows_selected,]$gene_id
      up_peak <- subset(mmAnno_up(), gene_id %in% gene)
      down_peak <- subset(mmAnno_down(), gene_id %in% gene)
      mcols(up_peak) <- DataFrame(Group = "up")
      mcols(down_peak) <- DataFrame(Group = "down")
      peak <- c(up_peak,down_peak)
      y <- as.data.frame(peak)
      print(y)
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
        print(name)
        files <- c(files,file)
      }
      names(files)<-name
      print(input$int_trackplot_additional1[[1, 'datapath']])
      print(files)
      return(files)
    }
  })
  
  int_data_track <- reactive({
      if(!is.null(input$RP_table_rows_selected)){
        library(Gviz)
        y <- int_goi_promoter_position()
        gene_position <- int_goi_gene_position()
        chr <- gene_position$seqnames
        gen <- ref()
        txdb <- txdb()
        grtrack <- GeneRegionTrack(txdb,
                                   chromosome = chr, name = "UCSC known genes",
                                   transcriptAnnotation = "tx_name",genome = gen,
                                   background.title = "grey",cex = 1.25)
        bw_files <- bws()
        cond1 <- unique(gsub("\\_.+$", "", names(bw_files)))[1]
        cond2 <- unique(gsub("\\_.+$", "", names(bw_files)))[2]
        if(!is.null(int_track_additional_files())) bw_files <- c(bw_files, int_track_additional_files())
        df <- list()
        for(name in names(bw_files)){
          cond <- gsub("\\_.+$", "", name)
          if(cond == cond1) col = "black"
          if(cond == cond2) col = "darkred"
          if(cond != cond1 && cond != cond2) col = "darkblue"
          df[[name]] <- DataTrack(range = bw_files[[name]], type = "l",genome = gen,
                                  name = gsub("\\..+$", "", name), window = -1,
                                  chromosome = chr, background.title = col,cex = 1.25,
                                  col.histogram = "darkblue", 
                                  fill.histogram = "darkblue",)
        }
        df[["grtrack"]] <- grtrack
        return(df)
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
                             chromosome = chr) 
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
        plotTracks(list(gtrack(), int_highlight_trackplot()),
                   from = input$int_igv_uprange[1], 
                   to = input$int_igv_uprange[2],ylim=c(0,input$int_igv_ylim),
                   type="hist")
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
            files[["Growing-p53_1"]] <- "data/peakcall/Growing-p53_1_peaks.narrowPeak"
            files[["Growing-p53_2"]] <- "data/peakcall/Growing-p53_2_peaks.narrowPeak"
            files2 <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
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
          files2 <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
          names(files2)<-name
          print(files2)
          return(files2)
        }
  })
  
  bws_venn <- reactive({
      if(is.null(input$file_venn1)){
        if(input$goButton_venn > 0 ){
          df<-list()
          df[["Ctrl_1"]] <- "data/bigwig/Growing-p53_1.BigWig"
          df[["Ctrl_2"]] <- "data/bigwig/Growing-p53_2.BigWig"
          df[["Sen_1"]] <- "data/bigwig/OIS-p53_1.BigWig"
          df[["Sen_2"]] <- "data/bigwig/OIS-p53_2.BigWig"
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
    ol <- findOverlapsOfPeaks(Venn_peak_call_files())
    names(ol$peaklist) <- gsub("///",":",names(ol$peaklist))
    return(ol)
  })
  
  output$venn <- renderPlot({
    withProgress(message = "Venn diagram, takes a few minutes",{
    if(!is.null(Venn_peak_call_files())){
      makeVennDiagram(venn_overlap()) 
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
        makeVennDiagram(venn_overlap()) 
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
    if(input$intersect_select != "not_selected" && 
       !is.null(bws_venn()) && input$Species_venn != "not selected"){
    Glist <- GRangesList()
    Glist[[input$intersect_select]] <- selected_grange()
    genomicElementDistribution(Glist, 
                               TxDb = txdb_venn()) 
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
    gene_IDs<-AnnotationDbi::select(org_venn(),keys = my.symbols,
                                    keytype = "ENTREZID",
                                    columns = c("SYMBOL", "ENTREZID"))
    colnames(gene_IDs) <- c("geneId","NearestGene")
    data <- merge(gene_IDs,overlaps.anno,by="geneId")
    data <- data[,2:11] %>% distinct(peakNames, .keep_all = T)
    return(data)
    })
  })
  output$selected_intersect_annotation <- DT::renderDT({
    if(input$intersect_select != "not_selected" && input$Species_venn != "not selected"){
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
        Glist <- GRangesList()
        Glist[[input$intersect_select]] <- selected_grange()
        pdf(file, height = pdf_height, width = pdf_width)
        genomicElementDistribution(Glist, 
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
    numericInput("peak_venn_range","Heatmap y-axis range",value=10,step=2)
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
    feature.recentered <- reCenterPeaks(venn2, width=4000)
    cvglists <- peak_venn_cvglists()
    feature.center <- reCenterPeaks(venn2, width=1)
    print(venn2)
    print(feature.center)
    print(cvglists)
    sig <- featureAlignedSignal(cvglists, feature.center,
                                upstream=2000, downstream=2000,)
    return(sig)
  })
  
  peak_venn_feature_center <- reactive({
    venn2 <- peak_venn_grange()
    feature.center <- reCenterPeaks(venn2, width=1)
    sig <- peak_venn_sig()
    keep <- rowSums(sig[[2]]) > 0
    feature.center <- feature.center[keep]
    return(feature.center)
  })
  peak_venn_alinedHeatmap <- reactive({
    if(!is.null(input$peak_venn_range)){
      range <- c()
      for(n in length(names(bws_venn()))) {range <- c(range, input$peak_venn_range)}
      heatmap <- featureAlignedHeatmap(peak_venn_sig(), peak_venn_feature_center(),
                                       upstream=2000, downstream=2000,
                                       upper.extreme=range,color = c("white","red"))
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
    if(!is.null(input$intersect_select) && input$intersect_select != "not_selected"){
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
    switch (input$Species_venn,
            "Mus musculus (mm10)" = ref <- "mm10",
            "Homo sapiens (hg19)" = ref <- "hg19",
            "Homo sapiens (hg38)" = ref <- "hg38")
    return(ref)
  })
  
  
  
  observeEvent(input$selected_intersect_annotation_rows_selected, ({
    updateCollapse(session,id =  "int_result_collapse_panel", open="Trackplot_venn_panel")
  }))
  observeEvent(input$selected_intersect_annotation_rows_selected,({
    if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
      y <- goi_promoter_position_venn()
      gene_position <- goi_gene_position_venn()
      start_position <- min(c(y$start,gene_position$start))
      end_position <- max(c(y$end,gene_position$end))
      updateSliderInput(session,"igv_venn_uprange","Range:",value = c(start_position,end_position),
                        step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  }))
  
  output$igv_venn_uprange <- renderUI({
    if(!is.null(goi_promoter_position_venn()) && !is.null(goi_gene_position_venn())){
    y <- goi_promoter_position_venn()
    gene_position <- goi_gene_position_venn()
    start_position <- min(c(y$start,gene_position$start))
    end_position <- max(c(y$end,gene_position$end))
    sliderInput("igv_venn_uprange","Range:",value = c(start_position,end_position),
                step = 100, min = start_position - 10000, max = end_position + 10000)
    }
  })
  output$igv_venn_ylim <- renderUI({
    numericInput("igv_venn_ylim","peak range:", value = 10, min = 0)
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
        print(label_data)
      return(label_data)
    }
  })
  
  goi_gene_position_venn <- reactive({
    if(!is.null(goi_promoter_position_venn())){
      label_data <- goi_promoter_position_venn()$NearestGene
      print(label_data)
      gene_IDs<-AnnotationDbi::select(org_venn(),keys = label_data,
                                      keytype = "SYMBOL",
                                      columns = "ENTREZID")
      gene_position <- as.data.frame(subset(gene_position_venn(), gene_id %in% gene_IDs))
      return(gene_position)
    }
  })
  
  data_track_venn <-reactive({
      library(Gviz)
      y <- goi_promoter_position_venn()
      gene_position <- goi_gene_position_venn()
      chr <- gene_position$seqnames
      gen <- ref_venn()
      txdb <- txdb_venn()
      grtrack <- GeneRegionTrack(txdb,
                                 chromosome = chr, name = "UCSC known genes",
                                 transcriptAnnotation = "tx_name",genome = gen,
                                 background.title = "grey",cex = 1.25)
      bw_files <- bws_venn()
      df <- list()
      for(name in names(bw_files)){
        cond <- gsub("\\_.+$", "", name)
        df[[name]] <- DataTrack(range = bw_files[[name]], type = "l",genome = gen,
                                name = gsub("\\..+$", "", name), window = -1,
                                chromosome = chr, background.title = "black",cex = 1.25,
                                col.histogram = "darkblue", 
                                fill.histogram = "darkblue",)
      }
      df[["grtrack"]] <- grtrack
      return(df)
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
                               chromosome = chr)
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
  pwms_motif_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(preMotif_list_venn())){
      return(pwms(Species = input$Species_venn))
    }
  })
  
  venn_consensus <- reactive({
    files <- Venn_peak_call_files()
    consensus <- soGGi:::runConsensusRegions(GRangesList(files), "none")
    return(consensus)
  })
  
  enrich_motif_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(pwms_motif_venn()) && 
       !is.null(preMotif_list_venn())){
        return(MotifAnalysis(df= preMotif_list_venn(), 
                             Species = input$Species_venn, pwms = pwms_motif_venn(),
                             consensus = venn_consensus()))
    }
  })
  
  output$motif_venn_plot <- renderPlot({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      Motifplot(df2 = enrich_motif_venn(),padj = 0.05)
    }
  })
  
  
  output$download_motif_venn_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- Motifplot(df2 = enrich_motif_venn(),padj = 0.05)
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
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
      paste0("motif_enrichment_table",".txt")
    },
    content = function(file){write.table(motif_table_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  motif_table_venn <- reactive({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      df2 <- enrich_motif_venn()
      df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
      for(name in names(df2)){
        res <- df2[[name]]
        res <- dplyr::filter(res, X1 > -log10(0.05))
        res <- res %>% dplyr::arrange(-X1.1)
        df <- rbind(df, res)
      }
      colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                        "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                        "Group")
      df$padj <- 10^(-df$negLog10Padj)
      return(df)
    }
  })
  
  output$motif_venn_result <- DT::renderDT({
    if(input$motifButton_venn > 0 && !is.null(enrich_motif_venn())){
      motif_table_venn()
    }
  })
  
  promoter_motif_region_venn <- reactive({
    target_motif <- motif_table_venn()[input$motif_venn_result_rows_selected,]
    if(input$motifButton_venn > 0 && !is.null(input$motif_venn_result_rows_selected)){
        res <- MotifRegion(df= preMotif_list_venn(),
                           target_motif = target_motif,Species = input$Species_venn)  
      return(res)
    }
  })
  
  output$promoter_motif_region_table <- renderDataTable({
    target_motif <- motif_table_venn()[input$motif_venn_result_rows_selected,]
    if(input$motifButton_venn > 0 && !is.null(input$motif_venn_result_rows_selected)){
      promoter_motif_region_venn()
    }
  })
  
  output$download_promoter_motif_venn_region = downloadHandler(
    filename = function() {
      paste0("motif_region",".txt")
    },
    content = function(file){write.table(promoter_motif_region_venn(), file, row.names = F, sep = "\t", quote = F)}
  )
  observeEvent(motif_table_venn()[input$motif_venn_result_rows_selected,], ({
    updateCollapse(session,id =  "Promoter_motif_venn_collapse_panel", open="Promoter_motif_region_venn_panel")
  }))
  
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
      print(input$venn_whichGroup2)
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
      switch (input$Species_venn,
              "Homo sapiens (hg38)" = source <- "TxDb.Hsapiens.UCSC.hg38.knownGene",
              "Homo sapiens (hg19)" = source <- "TxDb.Hsapiens.UCSC.hg19.knownGene",
              "Mus musculus (mm10)" = source <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
      )
      for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g),source)
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
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            print(data)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
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
          pdf_height <- 5
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 5
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
      print(set_list)
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
  
  ## Enrichment viewer------
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
        files[["Growing_1"]] <- "data/peakcall/Growing-p53_1_peaks.narrowPeak"
        files[["Growing_2"]] <- "data/peakcall/Growing-p53_2_peaks.narrowPeak"
        files[["OIS_1"]] <- "data/peakcall/OIS-p53_1_peaks.narrowPeak"
        files[["OIS_2"]] <- "data/peakcall/OIS-p53_2_peaks.narrowPeak"
        files2 <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
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
      files2 <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
      names(files2)<-name
      print(files2)
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
    switch (input$Species_enrich,
            "Homo sapiens (hg38)" = source <- "TxDb.Hsapiens.UCSC.hg38.knownGene",
            "Homo sapiens (hg19)" = source <- "TxDb.Hsapiens.UCSC.hg19.knownGene",
            "Mus musculus (mm10)" = source <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    )
    for(name in names(data)){
      res = rGREAT::great(gr = data[[name]],gene_sets = gene_list_for_enrichment_genome(H_t2g),source)
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
        print(data)
        idx <- order(data[["x"]], decreasing = FALSE)
        data$Description <- factor(data$Description,
                                   levels=rev(unique(data$Description[idx])))
        p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="p_adjust_hyper",size="fold_enrichment_hyper"))+
                        geom_point() +
                        scale_color_continuous(low="red", high="blue",
                                               guide=guide_colorbar(reverse=TRUE)) +
                        scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
                        scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
      }else p1 <- NULL
      p <- plot_grid(p1, nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment_enrich <- renderPlot({
    if(input$Species_enrich != "not selected"){
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
          pdf_width <- 5
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
    if(input$Species_enrich != "not selected"){
      as.data.frame(enrichment_1_1_enrich())
    }
  })
  
  
  
  output$whichGeneSet_enrich <- renderUI({
    if(!is.null(input$intersection_enrich_fun)){
      group <- input$intersection_enrich_fun
      set_list <- unique(dplyr::filter(as.data.frame(enrichment_1_1_enrich()), Group == group)$id)
      print(set_list)
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
  
  pwms_motif_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(Enrich_peak_call_files())){
      return(pwms(Species = input$Species_enrich))
    }
  })
  
  enrich_consensus <- reactive({
    files <- Enrich_peak_call_files()
    consensus <- soGGi:::runConsensusRegions(GRangesList(files), "none")
    return(consensus)
  })
  
  enrich_motif_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(pwms_motif_enrich()) && 
       !is.null(Enrich_peak_call_files())){
      return(MotifAnalysis(df= Enrich_peak_call_files(), 
                           Species = input$Species_enrich, pwms = pwms_motif_enrich()))
    }
  })
  
  output$motif_enrich_plot <- renderPlot({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      Motifplot(df2 = enrich_motif_enrich(),padj = 0.05)
    }
  })
  
  
  output$download_motif_enrich_plot = downloadHandler(
    filename = function() {
      paste0("motif_enrichment_plot",".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- Motifplot(df2 = enrich_motif_enrich(),padj = 0.05)
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6
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
      paste0("motif_enrichment_table",".txt")
    },
    content = function(file){write.table(motif_table_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  motif_table_enrich <- reactive({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      df2 <- enrich_motif_enrich()
      df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
      for(name in names(df2)){
        res <- df2[[name]]
        res <- dplyr::filter(res, X1 > -log10(0.05))
        res <- res %>% dplyr::arrange(-X1.1)
        df <- rbind(df, res)
      }
      colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                        "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                        "Group")
      df$padj <- 10^(-df$negLog10Padj)
      return(df)
    }
  })
  
  output$motif_enrich_result <- DT::renderDT({
    if(input$motifButton_enrich > 0 && !is.null(enrich_motif_enrich())){
      motif_table_enrich()
    }
  })
  
  promoter_motif_region_enrich <- reactive({
    target_motif <- motif_table_enrich()[input$motif_enrich_result_rows_selected,]
    if(input$motifButton_enrich > 0 && !is.null(input$motif_enrich_result_rows_selected)){
      res <- MotifRegion(df= Enrich_peak_call_files(),
                         target_motif = target_motif,Species = input$Species_enrich)  
      return(res)
    }
  })
  
  output$promoter_motif_region_table <- renderDataTable({
    target_motif <- motif_table_enrich()[input$motif_enrich_result_rows_selected,]
    if(input$motifButton_enrich > 0 && !is.null(input$motif_enrich_result_rows_selected)){
      promoter_motif_region_enrich()
    }
  })
  
  output$download_promoter_motif_enrich_region = downloadHandler(
    filename = function() {
      paste0("motif_region",".txt")
    },
    content = function(file){write.table(promoter_motif_region_enrich(), file, row.names = F, sep = "\t", quote = F)}
  )
  observeEvent(motif_table_enrich()[input$motif_enrich_result_rows_selected,], ({
    updateCollapse(session,id =  "Promoter_motif_enrich_collapse_panel", open="Promoter_motif_region_enrich_panel")
  }))
  

  
})